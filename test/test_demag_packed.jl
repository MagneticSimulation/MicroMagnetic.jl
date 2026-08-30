#Tests for the packed-FFT demag (tensors stored on the (y,z) parity
#fundamental domain, see pack_demag_tensor in src/micro/demag.jl).
#
#The reference is DirectDemag (add_demag(...; fft=false)), which builds the
#dense N x N demag matrix from the *same* mesh kernels but with a completely
#independent sign/index logic (no FFT, no padding, no parity tricks).  It
#supports the same macro PBC, so every case below compares
#   field  :  max |H_fft - H_direct|  <=  tol * max |H_direct|
#   energy :  sum of the per-site demag energies
#for random magnetization textures and spatially varying Ms, evaluated twice
#with different m (the second evaluation guards the reused padded buffers
#against stale/garbage padding).

using MicroMagnetic
using Test
using Random

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

function compare_demag_fft_direct(; nx, ny, nz, dx=2e-9, dy=2e-9, dz=2e-9,
                                  Nx=0, Ny=0, Nz=0, seed=1234)
    ntotal = nx * ny * nz
    mesh = FDMesh(; dx=dx, dy=dy, dz=dz, nx=nx, ny=ny, nz=nz)

    sim_fft = Sim(mesh)
    sim_dir = Sim(mesh)
    for sim in (sim_fft, sim_dir)
        Random.seed!(seed)
        m = 2.0 .* rand(3 * ntotal) .- 1.0
        init_m0(sim, m; norm=false)
        Random.seed!(seed + 7919)
        set_Ms(sim, 6e5 .+ 3e5 .* rand(ntotal))
    end

    add_demag(sim_fft; Nx=Nx, Ny=Ny, Nz=Nz)
    add_demag(sim_dir; fft=false, Nx=Nx, Ny=Ny, Nz=Nz)

    for trial in 1:2
        if trial == 2 #different texture for the second call
            for sim in (sim_fft, sim_dir)
                Random.seed!(seed + 104729)
                m2 = 2.0 .* rand(3 * ntotal) .- 1.0
                init_m0(sim, m2; norm=false)
            end
        end

        MicroMagnetic.effective_field(sim_fft, sim_fft.spin, 0.0)
        MicroMagnetic.effective_field(sim_dir, sim_dir.spin, 0.0)

        T = eltype(sim_fft.spin)
        h_fft = Array(sim_fft.field)
        h_dir = Array(sim_dir.field)
        scale = maximum(abs.(h_dir))
        tol = T == Float32 ? 2e-3 : 1e-8
        @test maximum(abs.(h_fft .- h_dir)) <= tol * max(scale, 1e-12)

        e_fft = sum(Array(sim_fft.interactions[1].energy))
        e_dir = sum(Array(sim_dir.interactions[1].energy))
        @test isapprox(e_fft, e_dir; rtol=10 * tol, atol=abs(scale) * 1e-12)
    end
end

function test_demag_packed()
    #3D, anisotropic cells (dx != dy != dz), varying Ms
    compare_demag_fft_direct(nx=5, ny=4, nz=3, dx=2e-9, dy=3e-9, dz=4e-9, seed=11)
    #thin film (nz=1: the c-parity axis of the packing degenerates)
    compare_demag_fft_direct(nx=8, ny=6, nz=1, dx=3e-9, dz=5e-9, seed=22)
    #small dims -> odd padded sizes (nx_fft=3, ny_fft=5, nz_fft=3)
    compare_demag_fft_direct(nx=2, ny=3, nz=2, seed=33)
    #even padded sizes everywhere -> Nyquist lines are exercised
    compare_demag_fft_direct(nx=4, ny=4, nz=4, seed=44)
    #macro PBC on x and y (DirectDemag replicates the same kernels)
    compare_demag_fft_direct(nx=4, ny=3, nz=2, Nx=2, Ny=3, Nz=1, seed=55)
    #monolayer: yz (and xz) tensor components are physically ~0, which stresses
    #the tolerance logic of the packed path against roundoff-level garbage
    compare_demag_fft_direct(nx=6, ny=5, nz=1, seed=66)
end

@using_gpu()
test_functions("DemagPacked", test_demag_packed)
