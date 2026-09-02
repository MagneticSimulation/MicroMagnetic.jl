#Tests for the true 2D-periodic demag for the direction pairs xy/xz/yz: the
#mesh declares the periodic pair (FDMesh(...; pbc=...)) and the bare add_demag
#dispatches to pbc2d; see src/micro/demag_pbc2d.jl.  References:
#  * uniform M in the periodic plane -> H = 0   (no bound charges)
#  * uniform M along the open axis  -> H = -M   (the film demag; the DC column)
#  * brute force: explicit image sum over the 2D image lattice
#  * macro (macroPBC=true, Nx=N, Ny=N) converges to pbc2d

using MicroMagnetic
using Random
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

const P2_NX, P2_NY, P2_NZ = 6, 5, 3

#pair -> (the two periodic axes, the open axis) as 1=x, 2=y, 3=z
pair_axes(pair) = pair == 1 ? (1, 2, 3) : pair == 2 ? (1, 3, 2) : (2, 3, 1)
pbc2d_mesh(pair::Int) = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=P2_NX, ny=P2_NY,
                               nz=P2_NZ, pbc = pair == 1 ? "xy" : pair == 2 ? "xz" : "yz")

function pbc2d_build(pair::Int, m; Ic::Int=2, Jc::Int=2, macroN::Int=0)
    sim = Sim(pbc2d_mesh(pair))
    set_Ms(sim, 8e5)
    init_m0(sim, m; norm=false)
    if macroN > 0
        N = pair == 1 ? (macroN, macroN, 0) : pair == 2 ? (macroN, 0, macroN) :
            (0, macroN, macroN)
        add_demag(sim; macroPBC=true, Nx=N[1], Ny=N[2], Nz=N[3])
    else
        add_demag(sim; Ic=Ic, Jc=Jc)
    end
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    return Array(sim.field)
end

#explicit image sum over the 2D image lattice (|p| <= P, |q| <= Q) of the
#Newell tensor at signed displacements; the field convention matches DirectDemag
function pbc2d_bruteforce_field(m; pair=1, P=3, Q=3)
    nx, ny, nz = P2_NX, P2_NY, P2_NZ
    maxs = 4e-9
    d = (Float64(2e-9 / maxs), Float64(3e-9 / maxs), Float64(4e-9 / maxs))
    a1, a2, ao = pair_axes(pair)
    n1 = a1 == 1 ? nx : a1 == 2 ? ny : nz
    n2 = a2 == 1 ? nx : a2 == 2 ? ny : nz
    mu0_Ms = fill(MicroMagnetic.mu_0 * 8e5, nx * ny * nz)
    h = zeros(3 * nx * ny * nz)
    nxy = nx * ny
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = nxy * (k - 1) + nx * (j - 1) + (i - 1)
        t = (i, j, k)
        for q in -Q:Q, p in -P:P, kp in 1:nz, jp in 1:ny, ip in 1:nx
            idp = nxy * (kp - 1) + nx * (jp - 1) + (ip - 1)
            s = (ip, jp, kp)
            o1 = (s[a1] - t[a1] + p * n1) * d[a1]
            o2 = (s[a2] - t[a2] + q * n2) * d[a2]
            o3 = (s[ao] - t[ao]) * d[ao]
            xa = a1 == 1 ? o1 : a2 == 1 ? o2 : o3
            ya = a1 == 2 ? o1 : a2 == 2 ? o2 : o3
            za = a1 == 3 ? o1 : a2 == 3 ? o2 : o3
            f = -mu0_Ms[idp+1] / MicroMagnetic.mu_0
            mx = m[3*idp+1]; my = m[3*idp+2]; mz = m[3*idp+3]
            xx = MicroMagnetic.demag_tensor_xx(xa, ya, za, d[1], d[2], d[3])
            xy = MicroMagnetic.demag_tensor_xy(xa, ya, za, d[1], d[2], d[3])
            xz = MicroMagnetic.demag_tensor_xz(xa, ya, za, d[1], d[2], d[3])
            yy = MicroMagnetic.demag_tensor_yy(xa, ya, za, d[1], d[2], d[3])
            yz = MicroMagnetic.demag_tensor_yz(xa, ya, za, d[1], d[2], d[3])
            zz = MicroMagnetic.demag_tensor_zz(xa, ya, za, d[1], d[2], d[3])
            h[3*id+1] += f * (xx*mx + xy*my + xz*mz)
            h[3*id+2] += f * (xy*mx + yy*my + yz*mz)
            h[3*id+3] += f * (xz*mx + yz*my + zz*mz)
        end
    end
    return h
end

function test_pbc2d_uniform()
    for pair in (1, 2, 3)
        a1, a2, ao = pair_axes(pair)
        #in-plane uniform: H = 0
        m = zeros(3 * P2_NX * P2_NY * P2_NZ)
        m[a1:3:end] .= 0.6
        m[a2:3:end] .= -0.8
        h = pbc2d_build(pair, m)
        tol = eltype(h) == Float32 ? 1e-2 : 1e-3
        @test maximum(abs, h) <= tol * 8e5
        #open-axis uniform: H = -M along the open axis, ~0 elsewhere
        m2 = zeros(3 * P2_NX * P2_NY * P2_NZ)
        m2[ao:3:end] .= 1.0
        h2 = pbc2d_build(pair, m2)
        target = zeros(3 * P2_NX * P2_NY * P2_NZ)
        target[ao:3:end] .= -8e5
        @test maximum(abs, h2 - target) <= tol * 8e5
    end
end

function test_pbc2d_bruteforce()
    for pair in (1, 2, 3)
        Random.seed!(41 + pair)
        m = randn(3 * P2_NX * P2_NY * P2_NZ)
        for c in 1:3
            idx = c:3:length(m)
            m[idx] .-= sum(m[idx]) / length(idx)
        end
        h = pbc2d_build(pair, m; Ic=2, Jc=2)
        h_bf = pbc2d_bruteforce_field(m; pair=pair, P=3, Q=3)
        scale = maximum(abs, h_bf)
        tol = eltype(h) == Float32 ? 1e-2 : 2e-3
        @test maximum(abs, h .- h_bf) <= tol * scale
    end
end

function test_pbc2d_macro_convergence()
    Random.seed!(11)
    m = randn(3 * P2_NX * P2_NY * P2_NZ)
    for c in 1:3
        idx = c:3:length(m)
        m[idx] .-= sum(m[idx]) / length(idx)
    end
    #tight pbc2d reference (Ic=Jc=6) so the macro truncation dominates
    h_pbc = pbc2d_build(1, m; Ic=6, Jc=6)
    scale = maximum(abs, h_pbc)
    e = [maximum(abs, pbc2d_build(1, m; macroN=N) - h_pbc) for N in (2, 4, 8)]
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
