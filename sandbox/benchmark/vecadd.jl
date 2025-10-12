using BenchmarkTools
using KernelAbstractions
using CUDA
using StaticArrays

# compute a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5
@kernel function vector_add5_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4), @Const(a5), c2, c3, c4, c5)
    i = @index(Global)
    @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i]
end

function vector_add5_cpu(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, c2::S, c3::S, c4::S, c5::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    kernel! = vector_add5_kernel!(CPU(), 512)
    kernel!(a, a1, a2, a3, a4, a5, c2, c3, c4, c5; ndrange=length(a))
    
    return nothing
end

function vector_add5_gpu(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, c2::S, c3::S, c4::S, c5::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    kernel! = vector_add5_kernel!(CUDABackend(), 512)
    kernel!(a, a1, a2, a3, a4, a5, c2, c3, c4, c5; ndrange=length(a))
    return nothing
end


function broadcast_version(a, a1, a2, a3, a4, a5, c2, c3, c4, c5)
    @. a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5
    return nothing
end

function broadcast_inplace!(a, a1, a2, a3, a4, a5, c2, c3, c4, c5)
    @inbounds for i in eachindex(a)
        a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i]
    end
    return nothing
end

function comprehensive_benchmark()
    
    sizes = [10^3, 10^4, 10^5, 10^6, 10^7]
    c2, c3, c4, c5 = 0.1, 0.2, 0.3, 0.4
    
    println("="^60)
    println("Comprehensive Vector Addition Benchmark")
    println("="^60)
    
    for size in sizes
        println("\n" * "="^50)
        println("Array size: ", size)
        println("="^50)
        
        
        a_cpu = rand(Float64, size)
        a1_cpu = rand(Float64, size)
        a2_cpu = rand(Float64, size)
        a3_cpu = rand(Float64, size)
        a4_cpu = rand(Float64, size)
        a5_cpu = rand(Float64, size)
        
        
        println("CPU Results:")
        
        
        a_temp = copy(a_cpu)
        b1 = @benchmark broadcast_version($a_temp, $a1_cpu, $a2_cpu, $a3_cpu, $a4_cpu, $a5_cpu, $c2, $c3, $c4, $c5)
        println("  Broadcast (@.):     ", BenchmarkTools.prettytime(median(b1).time))
        
        
        a_temp = copy(a_cpu)
        b2 = @benchmark broadcast_inplace!($a_temp, $a1_cpu, $a2_cpu, $a3_cpu, $a4_cpu, $a5_cpu, $c2, $c3, $c4, $c5)
        println("  Inplace loop:       ", BenchmarkTools.prettytime(median(b2).time))
        
        
        a_temp = copy(a_cpu)
        b3 = @benchmark vector_add5_cpu($a_temp, $a1_cpu, $a2_cpu, $a3_cpu, $a4_cpu, $a5_cpu, $c2, $c3, $c4, $c5)
        println("  KA CPU kernel:      ", BenchmarkTools.prettytime(median(b3).time))
        
        
        result1 = copy(a_cpu); broadcast_version(result1, a1_cpu, a2_cpu, a3_cpu, a4_cpu, a5_cpu, c2, c3, c4, c5)
        result2 = copy(a_cpu); vector_add5_cpu(result2, a1_cpu, a2_cpu, a3_cpu, a4_cpu, a5_cpu, c2, c3, c4, c5)
        @assert maximum(abs.(result1 .- result2)) < 1e-12 "CPU results don't match!"
        
        
        if CUDA.functional() && size >= 10^4  
            println("\nGPU Results:")
            
            
            a_gpu = CuArray(a_cpu)
            a1_gpu = CuArray(a1_cpu)
            a2_gpu = CuArray(a2_cpu)
            a3_gpu = CuArray(a3_cpu)
            a4_gpu = CuArray(a4_cpu)
            a5_gpu = CuArray(a5_cpu)
            
            
            vector_add5_gpu(a_gpu, a1_gpu, a2_gpu, a3_gpu, a4_gpu, a5_gpu, c2, c3, c4, c5)
            CUDA.synchronize()
            
            a_temp_gpu = copy(a_gpu)
            b4 = @benchmark (CUDA.@sync @. $a_temp_gpu = $a1_gpu + $c2*$a2_gpu + $c3*$a3_gpu + $c4*$a4_gpu + $c5*$a5_gpu)
            println("  GPU Broadcast:     ", BenchmarkTools.prettytime(median(b4).time))
            
            a_temp_gpu = copy(a_gpu)
            b5 = @benchmark (CUDA.@sync vector_add5_gpu($a_temp_gpu, $a1_gpu, $a2_gpu, $a3_gpu, $a4_gpu, $a5_gpu, $c2, $c3, $c4, $c5))
            println("  KA GPU kernel:     ", BenchmarkTools.prettytime(median(b5).time))
            
            result_gpu1 = copy(a_gpu); CUDA.@sync @. result_gpu1 = a1_gpu + c2*a2_gpu + c3*a3_gpu + c4*a4_gpu + c5*a5_gpu
            result_gpu2 = copy(a_gpu); vector_add5_gpu(result_gpu2, a1_gpu, a2_gpu, a3_gpu, a4_gpu, a5_gpu, c2, c3, c4, c5); CUDA.synchronize()
            @assert maximum(abs.(Array(result_gpu1) .- Array(result_gpu2))) < 1e-12 "GPU results don't match!"
            
            cpu_result = Array(result_gpu1)
            @assert maximum(abs.(result1 .- cpu_result)) < 1e-12 "CPU-GPU results don't match!"
        else
            if !CUDA.functional()
                println("\nGPU: Not available")
            else
                println("\nGPU: Skipping small array (size < 10^4)")
            end
        end
    end
    
end

comprehensive_benchmark()
