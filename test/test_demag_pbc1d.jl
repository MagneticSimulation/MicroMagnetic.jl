#Tests for the true 1D-periodic demag: the mesh declares the periodic axis
#(FDMesh(...; pbc=...)) and the bare add_demag dispatches to pbc1d; see
#src/micro/demag_pbc1d.jl.  References:
#  * brute force: explicit image sum |p| <= P along the periodic axis
#  * macro (macroPBC=true, Nx=N) converges to pbc1d
#  * uniform magnetization along the periodic axis: the true response is zero
#    (no bound charges); the computed field is the tail truncation error

using MicroMagnetic
using Random
using Test

if !isdefined(Main, :test_functions)
    include("test_utils.jl")
end

const P1_NX, P1_NY, P1_NZ = 8, 5, 3

pbc1d_mesh(axis::Int) = FDMesh(; dx=2e-9, dy=3e-9, dz=4e-9, nx=P1_NX, ny=P1_NY,
                               nz=P1_NZ, pbc = axis == 1 ? "x" : axis == 2 ? "y" : "z")

function pbc1d_build(axis::Int, m; Ic::Int=2, macroN::Int=0)
    sim = Sim(pbc1d_mesh(axis))
    set_Ms(sim, 8e5)
    init_m0(sim, m; norm=false)
    if macroN > 0
        N = axis == 1 ? (macroN, 0, 0) : axis == 2 ? (0, macroN, 0) : (0, 0, macroN)
        add_demag(sim; macroPBC=true, Nx=N[1], Ny=N[2], Nz=N[3])
    else
        add_demag(sim; Ic=Ic)
    end
    MicroMagnetic.effective_field(sim, sim.spin, 0.0)
    return Array(sim.field)
end

#explicit image sum |p| <= P along the periodic axis of the Newell tensor at
#signed displacements; the field convention matches DirectDemag
function pbc1d_bruteforce_field(m; axis=1, P=5)
    nx, ny, nz = P1_NX, P1_NY, P1_NZ
    maxs = 4e-9
    dxn = Float64(2e-9 / maxs); dyn = Float64(3e-9 / maxs); dzn = Float64(4e-9 / maxs)
    na = axis == 1 ? nx : axis == 2 ? ny : nz
    mu0_Ms = fill(MicroMagnetic.mu_0 * 8e5, nx * ny * nz)
    h = zeros(3 * nx * ny * nz)
    nxy = nx * ny
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = nxy * (k - 1) + nx * (j - 1) + (i - 1)
        for p in -P:P, kp in 1:nz, jp in 1:ny, ip in 1:nx
            idp = nxy * (kp - 1) + nx * (jp - 1) + (ip - 1)
            off = axis == 1 ? (ip - i + p * na) * dxn :
                  axis == 2 ? (jp - j + p * na) * dyn : (kp - k + p * na) * dzn
            xa = axis == 1 ? off : (ip - i) * dxn
            ya = axis == 2 ? off : (jp - j) * dyn
            za = axis == 3 ? off : (kp - k) * dzn
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

function test_pbc1d_uniform()
    for axis in (1, 2, 3)
        m = zeros(3 * P1_NX * P1_NY * P1_NZ)
        m[axis:3:end] .= 1.0
        #the uniform-along-axis state carries no bound charges, so the true
        #periodic response is zero; the computed field is the kernel tail
        #truncation error and must converge with Ic
        h2 = pbc1d_build(axis, m; Ic=2)
        h6 = pbc1d_build(axis, m; Ic=6)
        tol = eltype(h6) == Float32 ? 1e-4 : 1e-5
        @test maximum(abs, h6) <= tol * 8e5          #converged: H ~ 0
        @test maximum(abs, h2) <= 1e-3 * 8e5         #Ic=2 already small
    end
end

function test_pbc1d_bruteforce()
    for axis in (1, 2, 3)
        Random.seed!(21 + axis)
        m = randn(3 * P1_NX * P1_NY * P1_NZ)
        for c in 1:3
            idx = c:3:length(m)
            m[idx] .-= sum(m[idx]) / length(idx)
        end
        h = pbc1d_build(axis, m; Ic=2)
        h_bf = pbc1d_bruteforce_field(m; axis=axis, P=4)
        scale = maximum(abs, h_bf)
        tol = eltype(h) == Float32 ? 1e-2 : 2e-3
        @test maximum(abs, h .- h_bf) <= tol * scale
    end
end

function test_pbc1d_macro_convergence()
    Random.seed!(31)
    m = randn(3 * P1_NX * P1_NY * P1_NZ)
    for c in 1:3
        idx = c:3:length(m)
        m[idx] .-= sum(m[idx]) / length(idx)
    end
    #tight pbc1d reference (Ic=6) so the macro truncation dominates
    h_pbc = pbc1d_build(1, m; Ic=6)
    scale = maximum(abs, h_pbc)
    e = [maximum(abs, pbc1d_build(1, m; macroN=N) - h_pbc) for N in (2, 4, 8)]
    @test e[3] < e[1]                    #macro converges to pbc1d
    @test e[3] < 5e-3 * scale            #... and is close at N=8
    @test e[2] < 2e-2 * scale
end

function test_demag_pbc1d()
    test_pbc1d_uniform()
    test_pbc1d_bruteforce()
    test_pbc1d_macro_convergence()
end

@using_gpu()
test_functions("DemagPBC1D", test_demag_pbc1d)
