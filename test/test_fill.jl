using Test
using KernelAbstractions
using MicroMagnetic

# Interface unit tests for the O(1) uniform array.
@testset "Fill interface" begin
    f = MicroMagnetic.Fill(0.1, 7)

    @test size(f) == (7,)
    @test axes(f) == (Base.OneTo(7),)
    @test length(f) == 7
    @test eltype(f) == Float64
    @test f[1] == 0.1 && f[7] == 0.1
    @test f[3, 1] == 0.1  # linear indexing
    @test Base.IndexStyle(typeof(f)) == Base.IndexLinear()
    @test collect(f) == fill(0.1, 7)
    @test Array(f) == fill(0.1, 7)
    @test Array{Float32}(f) == fill(0.1f0, 7)
    @test sum(f) ≈ 0.7
    @test sum(f) == 0.1 * 7          # O(1), exact product
    @test maximum(f) == 0.1
    @test [x for x in f] == fill(0.1, 7)
    @test f == MicroMagnetic.Fill(0.1, 7)
    @test f != MicroMagnetic.Fill(0.2, 7)
    @test f != MicroMagnetic.Fill(0.1, 8)

    # broadcasting materialises a dense array with the right values
    bc = f .+ 1
    @test bc isa Array && bc == fill(1.1, 7)

    # immutability: every in-place write must throw, never silently corrupt
    @test_throws Exception copyto!(f, zeros(7))
    @test_throws Exception (f .= 0.3)
    @test_throws Exception fill!(f, 0.0)
    @test_throws Exception (f[1] = 0.5)
    @test sum(f) ≈ 0.7  # unchanged (0.1*7 rounds to 0.7000000000000001)

    # isbits: the property that lets a Fill travel into a GPU kernel by value
    @test isbits(f)
    @test !isbits(MicroMagnetic.Fill(:x, 3))

    # copyto! into a dense dest (the IO / Array() path) works and is exact
    dest = zeros(7)
    copyto!(dest, f)
    @test dest == fill(0.1, 7)

    # transpose used by average_m: b .* ms' stays lazy and correct
    b = reshape(collect(1.0:14.0), 2, 7)
    ms = MicroMagnetic.Fill(2.0, 7)
    @test sum(b .* ms'; dims=2)[:, 1] == [2 * sum(b[1, :]), 2 * sum(b[2, :])]

    # copyto!(Array, Fill) hook serves Array(sim.mu0_Ms)-style calls on N-dim fills
    f2 = MicroMagnetic.Fill{Int,2}(3, (2, 2))
    @test Array(f2) == [3 3; 3 3]
end

# kernel smoke test: a Fill as a kernel argument must produce the same result as
# a dense array with the same values.
@kernel function fill_read_kernel!(out, @Const(x), p, n::Int)
    I = @index(Global)
    out[I] = x[I] * p[I]
end

@testset "Fill in KA kernel" begin
    n = 64
    x = collect(1.0:n)
    out_dense = zeros(n)
    fill_read_kernel!(CPU(), 32)(out_dense, x, x, n; ndrange=n)
    @test out_dense == x .^ 2

    out_fill = zeros(n)
    fill_read_kernel!(CPU(), 32)(out_fill, MicroMagnetic.Fill(2.0, n), x, n; ndrange=n)
    @test out_fill == 2 .* x

    # bit-identical to the dense path for the same values
    xf = fill(2.0, n)
    out_dense2 = zeros(n)
    fill_read_kernel!(CPU(), 32)(out_dense2, xf, x, n; ndrange=n)
    @test out_dense2 == out_fill
end
