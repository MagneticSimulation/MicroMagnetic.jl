#Tests for the true 2D-periodic demag (add_demag(...; pbc2d=true), see
#src/micro/demag_pbc2d.jl).  References:
#  * uniform in-plane M  ->  H = 0        (no charges, in-plane uniform)
#  * uniform M_z         ->  H_z = -M_z  (the film demag; the k=0 column)
#  * brute force: explicit in-plane image sum of the Newell tensor (z open),
#    assembled independently of the spectral pipeline
#  * macro PBC (Nx=Ny=N images, z open) converges to pbc2d as N grows

using MicroMagnetic
using Random
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

#explicit in-plane image sum |p|,|q| <= P of the Newell tensor at signed
#displacements; the field convention matches DirectDemag (h = -sum N * mu0_Ms/mu_0 * m)
function pbc2d_bruteforce_field(sim; P=5)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    maxs = max(mesh.dx, mesh.dy, mesh.dz)
    dxn = Float64(mesh.dx / maxs); dyn = Float64(mesh.dy / maxs); dzn = Float64(mesh.dz / maxs)
    m = Array(sim.spin)
    mu0_Ms = Array(sim.mu0_Ms)
    h = zeros(3 * nx * ny * nz)
    nxy = nx * ny
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = nxy * (k - 1) + nx * (j - 1) + (i - 1)
        for q in -P:P, p in -P:P, kp in 1:nz, jp in 1:ny, ip in 1:nx
            idp = nxy * (kp - 1) + nx * (jp - 1) + (ip - 1)
            xa = (ip - i + p * nx) * dxn; ya = (jp - j + q * ny) * dyn; za = (kp - k) * dzn
            f = -mu0_Ms[idp+1] / MicroMagnetic.mu_0
            mx = m[3*idp+1]; my = m[3*idp+2]; mz = m[3*idp+3]
            xx = MicroMagnetic.demag_tensor_xx(xa, ya, za, dxn, dyn, dzn)
            xy = MicroMagnetic.demag_tensor_xy(xa, ya, za, dxn, dyn, dzn)
            xz = MicroMagnetic.demag_tensor_xz(xa, ya, za, dxn, dyn, dzn)
            yy = MicroMagnetic.demag_tensor_yy(xa, ya, za, dxn, dyn, dzn)
            yz = MicroMagnetic.demag_tensor_yz(xa, ya, za, dxn, dyn, dzn)
            zz = MicroMagnetic.demag_tensor_zz(xa, ya, za, dxn, dyn, dzn)
            h[3*id+1] += f * (xx*mx + xy*my + xz*mz)
            h[3*id+2] += f * (xy*mx + yy*my + yz*mz)
            h[3*id+3] += f * (xz*mx + yz*my + zz*mz)
        end
    end
    return h
end

function test_pbc2d_uniform()
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=6, ny=5, nz=3)
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    init_m0(sim, (0.6, -0.8, 0.0); norm=false)
    add_demag(sim; pbc2d=true)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    Tf = eltype(sim.spin)
    tol = Tf == Float32 ? 1e-5 : 1e-9
    @test maximum(abs.(Array(sim.field))) <= tol * 8e5   #in-plane uniform: H = 0

    init_m0(sim, (0.0, 0.0, 1.0); norm=false)
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    h = Array(sim.field)
    @test maximum(abs.(h[1:3:end])) <= tol * 8e5          #H_x = 0
    @test maximum(abs.(h[2:3:end])) <= tol * 8e5          #H_y = 0
    @test maximum(abs.(h[3:3:end] .+ 8e5)) <= tol * 8e5   #H_z = -M_z
end

function test_pbc2d_bruteforce()
    nx, ny, nz = 8, 6, 2
    mesh = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=nx, ny=ny, nz=nz)
    Random.seed!(11)
    m = randn(3 * nx * ny * nz)
    for c in 1:3
        idx = c:3:length(m)
        m[idx] .-= sum(m[idx]) / length(idx)
    end
    #random texture, two tail ranges: the default (Ic=2, ~1e-3 field error)
    #and a tight one (Ic=4, ~2e-4)
    for (Ic, tol64) in ((2, 1e-2), (4, 1e-3))
        sim = Sim(mesh)
        set_Ms(sim, 8e5)
        init_m0(sim, m; norm=false)
        add_demag(sim; pbc2d=true, Ic=Ic, Jc=Ic)
        MicroMagnetic.effective_field(sim, sim.spin, 0.0)
        h = Array(sim.field)
        h_bf = pbc2d_bruteforce_field(sim; P=Ic + 1)
        scale = maximum(abs, h_bf)
        Tf = eltype(sim.spin)
        tol = Tf == Float32 ? 5 * tol64 : tol64
        @test maximum(abs, h .- h_bf) <= tol * scale
    end
end

function test_pbc2d_macro_convergence()
    nx, ny, nz = 8, 6, 2
    dx, dy, dz = 2e-9, 3e-9, 4e-9
    mesh = FDMesh(; dx=dx, dy=dy, dz=dz, nx=nx, ny=ny, nz=nz)
    Random.seed!(11)
    m = randn(3 * nx * ny * nz)
    for c in 1:3
        idx = c:3:length(m)
        m[idx] .-= sum(m[idx]) / length(idx)
    end
    function build(; pbc2d=false, N=0, Ic=2)
        sim = Sim(mesh)
        set_Ms(sim, 8e5)
        init_m0(sim, m; norm=false)
        add_demag(sim; pbc2d=pbc2d, Nx=N, Ny=N, Nz=0, Ic=Ic, Jc=Ic)
        MicroMagnetic.effective_field(sim, sim.spin, 0.0)
        return Array(sim.field)
    end
    #tight pbc2d reference (Ic=6, tail error ~6e-5) so the macro truncation
    #dominates the comparison
    h_pbc = build(pbc2d=true, Ic=6)
    scale = maximum(abs.(h_pbc))
    e = [maximum(abs.(build(N=N) .- h_pbc)) for N in (2, 4, 8)]
    @test e[3] < e[1]                    #macro converges to pbc2d
    @test e[3] < 5e-3 * scale            #... and is close at N=8
    @test e[2] < 2e-2 * scale
end

function test_demag_pbc2d()
    test_pbc2d_uniform()
    test_pbc2d_bruteforce()
    test_pbc2d_macro_convergence()
end

@using_gpu()
test_functions("DemagPBC2D", test_demag_pbc2d)
