module CUDAExt

using MicroMagnetic
using CUDA
using CUDA.CUSPARSE
using CUDA.CUFFT

CUDA.allowscalar(false)

#Backends with a raw in-place c2r transform (see the comments in
#src/micro/demag.jl): the real field overwrites the complex spectrum buffer
#(padded to 2*lenx = nx_fft+2 real rows per component), eliminating the
#separate h_pad allocation.  The raw cufft transform is unnormalized; 1/N is
#folded into the packed tensors by init_demag.
MicroMagnetic.inplace_inverse(::CuArray) = true
MicroMagnetic.inplace_inverse(::CuArray{Complex{T}}) where {T} = T === Float64 || T === Float32

mutable struct CUFFTInplacePlan
    handle::CUFFT.cufftHandle
    T::Type
    function CUFFTInplacePlan(handle::CUFFT.cufftHandle, T::Type)
        p = new(handle, T)
        finalizer(p) do pl
            CUFFT.cufftDestroy(pl.handle)
        end
        return p
    end
end

function MicroMagnetic.make_inplace_plan(M_pad::CuArray{Complex{T},4}, nx_fft::Int) where {T}
    lenx, ny, nz = size(M_pad)[1:3]
    plan = Ref{CUFFT.cufftHandle}(Cint(0))
    n = Cint[nz, ny, nx_fft]            #C order: the halved (first) dim goes last
    inembed = Cint[nz, ny, lenx]        #complex input dims
    onembed = Cint[nz, ny, 2 * lenx]    #padded real output dims
    ptype = T === Float64 ? CUFFT.CUFFT_Z2D : CUFFT.CUFFT_C2R
    CUFFT.cufftPlanMany(plan, Cint(3), n, inembed, Cint(1), Cint(lenx * ny * nz),
                        onembed, Cint(1), Cint(2 * lenx * ny * nz), ptype, Cint(3))
    return CUFFTInplacePlan(plan[], T)
end

function MicroMagnetic.inv_transform!(h_pad, M_pad::CuArray, h_plan::CUFFTInplacePlan,
                                      nx_fft::Int)
    CUFFT.cufftSetStream(h_plan.handle, CUDA.stream())
    ptr_real = reinterpret(CuPtr{Float64}, pointer(M_pad))
    if h_plan.T === Float64
        CUFFT.cufftExecZ2D(h_plan.handle, pointer(M_pad), ptr_real)
    else
        CUFFT.cufftExecC2R(h_plan.handle, pointer(M_pad), ptr_real)
    end
    return nothing
end

function set_cuda_backend()
    MicroMagnetic.all_backends[1] = CUDA.CUDABackend()
    MicroMagnetic.set_backend("cuda")
    MicroMagnetic.GPUSparseMatrixCSC[] = CuSparseMatrixCSC
    MicroMagnetic.GPUSparseMatrixCSR[] = CuSparseMatrixCSR
    MicroMagnetic.GPUMatrix[] = CuMatrix
    return nothing
end

set_cuda_backend()

include(joinpath(@__DIR__, "..", "src", "precompile.jl"))

function __init__()
    return set_cuda_backend()
end

end
