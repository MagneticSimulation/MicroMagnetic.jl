#Tests for the true 3D-periodic demag (the mesh pbc="xyz" dispatches to
#pbc3d via the bare add_demag, see
#src/micro/demag_pbc3d.jl).  The reference is analytic:
#  * the uniform (k=0) mode produces no field  (tin-foil convention),
#  * a single Fourier mode M = Ms*cos(2*pi*u/L) along axis ax with polarization
#    p gives H = -Ms*sinc^2(pi/n_ax)*cos(2*pi*u/L) along p if p is parallel to
#    the mode direction, and H = 0 if p is transverse (no magnetic charges),
#  * the macro PBC (image summation) converges to pbc3d for zero-mean
#    textures as the number of images grows (the k!=0 modes agree; only the
#    k=0 convention differs, which is projected out here).

using MicroMagnetic
using Test
using Random

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

#spin array of a single Fourier mode along axis ax (1=x, 2=y, 3=z),
#polarized along axis pol
function mode_spin(nx, ny, nz, dx, dy, dz, ax::Int, pol::Int)
    d = (dx, dy, dz)
    n = (nx, ny, nz)
    m = zeros(3 * nx * ny * nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        I0 = nx * ny * (k - 1) + nx * (j - 1) + (i - 1)
        u = ((i, j, k)[ax] - 0.5) * d[ax]
        m[3 * I0 + pol] = cos(2 * pi * u / (n[ax] * d[ax]))
    end
    return m
end

function test_pbc3d_analytic(T = nothing)
    nx, ny, nz = 8, 6, 4
    dx, dy, dz = 2e-9, 3e-9, 5e-9
    n = (nx, ny, nz)
    d = (dx, dy, dz)
    Ms = 8e5
    mesh = FDMesh(; dx=dx, dy=dy, dz=dz, nx=nx, ny=ny, nz=nz, pbc="xyz")
    sim = Sim(mesh)
    set_Ms(sim, Ms)
    add_demag(sim)

    for ax in 1:3, pol in 1:3
        init_m0(sim, mode_spin(nx, ny, nz, dx, dy, dz, ax, pol); norm=false)
        MicroMagnetic.effective_field(sim, sim.spin, 0.0)
        h = Array(sim.field)
        S = (sin(pi / n[ax]) / (pi / n[ax]))^2
        href = zeros(3 * nx * ny * nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            I0 = nx * ny * (k - 1) + nx * (j - 1) + (i - 1)
            u = ((i, j, k)[ax] - 0.5) * d[ax]
            amp = ax == pol ? -Ms * S * cos(2 * pi * u / (n[ax] * d[ax])) : 0.0
            href[3 * I0 + pol] = amp
        end
        Tf = eltype(sim.spin)
        tol = Tf == Float32 ? 1e-4 : 1e-9
        @test maximum(abs.(h .- href)) <= tol * Ms
    end
end

function test_pbc3d_uniform()
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=6, ny=5, nz=3, pbc="xyz")
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, (0.3, -0.5, 0.81); norm=false)
    add_demag(sim)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    Tf = eltype(sim.spin)
    tol = Tf == Float32 ? 1e-3 : 1e-9
    @test maximum(abs.(Array(sim.field))) <= tol * 8e5  #tin-foil: <H> = 0
end

#The band-limited spectral kernel and the top-hat (Newell/macro-PBC) kernel are
#two discretizations of the same continuum problem: for a fixed PHYSICAL mode
#(wavelength = Lx/2) their per-mode amplitudes differ by 1 - sinc^2(pi*fx/nx),
#which shrinks ~1/nx^2 under grid refinement.  Both describe the STACKED
#geometry (the grid repeated with no vacuum), so the macro reference uses
#images in all three directions (Nx=Ny=Nz=6, converged to <2% here).
function test_pbc3d_refinement()
    Ms = 8e5
    errs = Float64[]
    for nx in (12, 24)
        dx = 2e-9
        mesh = FDMesh(; dx=dx, dy=dx, dz=dx, nx=nx, ny=2, nz=1, pbc="xyz")
        sim = Sim(mesh)
        set_Ms(sim, Ms)
        m = zeros(3 * nx * 2)
        for j in 1:2, i in 1:nx
            I0 = nx * (j - 1) + (i - 1)
            m[3 * I0 + 1] = cos(2pi * (i - 0.5) * 2 / nx)
        end
        init_m0(sim, m; norm=false)

        add_demag(sim)
        MicroMagnetic.effective_field(sim, sim.spin, 0.0)
        h_pbc = Array(sim.field)
        S = (sin(2pi / nx) / (2pi / nx))^2
        href = zeros(3 * nx * 2)
        for j in 1:2, i in 1:nx
            I0 = nx * (j - 1) + (i - 1)
            href[3 * I0 + 1] = -Ms * S * cos(2pi * (i - 0.5) * 2 / nx)
        end
        Tf = eltype(sim.spin)
        @test maximum(abs.(h_pbc .- href)) <= (Tf == Float32 ? 1e-4 : 1e-9) * Ms

        #macro reference (converged image sum, z open like the texture)
        sim2 = Sim(mesh)
        set_Ms(sim2, Ms)
        init_m0(sim2, m; norm=false)
        add_demag(sim2; macroPBC=true, Nx=6, Ny=6, Nz=6)   #same stacked geometry as pbc3d
        MicroMagnetic.effective_field(sim2, sim2.spin, 0.0)
        push!(errs, maximum(abs.(Array(sim2.field) .- h_pbc)))
    end
    @test errs[2] < errs[1] / 2.5        #the discretization difference shrinks ~1/nx^2
    @test errs[2] < 3e-2 * Ms            #... and is small at nx=24
end

function test_demag_pbc3d()
    test_pbc3d_uniform()
    test_pbc3d_analytic()
    test_pbc3d_refinement()
end

@using_gpu()
test_functions("DemagPBC3D", test_demag_pbc3d)
