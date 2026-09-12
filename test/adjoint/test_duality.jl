# Duality tests for the adjoint layer.
#
#   ⟨Df·x, λ⟩ = ⟨x, Dfᵀ·λ⟩
#
# 1. Per-operator 3N duality: production forward field kernels vs the VJP
#    layer (`apply_KT!`).  Every supported interaction is linear in m, so
#    K·u is just the field computed at u — no tangent-space subtleties, and
#    open-boundary one-sided bonds / periodic wraps are exercised exactly.
# 2. Full-RHS duality: JVP from the existing LLGJacOperator (symbolic, 2N
#    tangent frame) vs `mul!(v, adjoint(op), λ)`; plus an independent check
#    against the transpose of the dense symbolic Jacobian.
# 3. Theta-space duality: central finite differences of the production field
#    in per-cell design parameters vs the grad_theta kernels.
# 4. loss_grad_kernel!: squared-error gradient with tangential projection.

using Test
using MicroMagnetic
using LinearAlgebra
using Random
using Printf
import KernelAbstractions

MicroMagnetic.set_backend("cpu")

function rand_unit_m0(N; seed = 1)
    Random.seed!(seed)
    m0 = randn(3N)
    for i in 1:N
        k = 3 * (i - 1) + 1
        nrm = sqrt(m0[k]^2 + m0[k+1]^2 + m0[k+2]^2)
        m0[k] /= nrm
        m0[k+1] /= nrm
        m0[k+2] /= nrm
    end
    return m0
end

function build_sim(interactions!::Function; nx = 5, ny = 4, nz = 3, pbc = "open",
                   mode = Float64)
    MicroMagnetic.set_precision(mode)
    mesh = FDMesh(nx = nx, ny = ny, nz = nz, dx = 5e-9, dy = 5e-9, dz = 5e-9, pbc = pbc)
    sim = Sim(mesh)
    set_Ms(sim, 8e5)
    interactions!(sim)
    return sim
end

@testset "adjoint: per-operator Kᵀ duality (3N)" begin
    cases = [
        ("exchange uniform", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_exch(sim, 1.3e-11)),
        ("exchange map", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_exch(sim, (x, y, z) -> 1.3e-11 + 0.3e-11 * sin(1e9 * x))),
        ("bulk DMI uniform (open)", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_dmi(sim, 3e-3; type = "bulk")),
        ("bulk DMI map (open)", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_dmi(sim, (x, y, z) -> 3e-3 + 0.6e-3 * cos(1e9 * y); type = "bulk")),
        ("bulk DMI (pbc x)", Dict(:nx => 5, :ny => 4, :nz => 3, :pbc => "x"),
         sim -> add_dmi(sim, 3e-3; type = "bulk")),
        ("bulk DMI (pbc xyz)", Dict(:nx => 5, :ny => 4, :nz => 3, :pbc => "xyz"),
         sim -> add_dmi(sim, 3e-3; type = "bulk")),
        ("interfacial DMI (open)", Dict(:nx => 5, :ny => 4, :nz => 2),
         sim -> add_dmi(sim, 3e-3; type = "interfacial")),
        ("interfacial DMI (pbc xy)", Dict(:nx => 5, :ny => 4, :nz => 2, :pbc => "xy"),
         sim -> add_dmi(sim, (x, y, z) -> 3e-3 + 0.5e-3 * sin(1e9 * x);
                        type = "interfacial")),
        ("uniaxial anisotropy", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_anis(sim, 5e4; axis = (0.2, 0.3, 1.0))),
        ("demag FFT (open)", Dict(:nx => 5, :ny => 4, :nz => 3),
         sim -> add_demag(sim)),
        ("demag pbc1d (pbc x)", Dict(:nx => 5, :ny => 4, :nz => 3, :pbc => "x"),
         sim -> add_demag(sim)),
        ("demag pbc3d (pbc xyz)", Dict(:nx => 5, :ny => 4, :nz => 3, :pbc => "xyz"),
         sim -> add_demag(sim)),
        ("interlayer exchange", Dict(:nx => 5, :ny => 4, :nz => 2),
         sim -> add_exch_int(sim, 1e-3)),
        ("interlayer DMI", Dict(:nx => 5, :ny => 4, :nz => 2),
         sim -> add_dmi_int(sim, (0.0, 0.0, 1e-3))),
    ]

    for (name, kwargs, inters!) in cases
        sim = build_sim(inters!; kwargs...)
        N = sim.n_total
        Random.seed!(11)
        u = randn(3N)
        w = randn(3N)

        MicroMagnetic.effective_field(sim, u, 0.0)           # K·u via the production path
        Ku = Float64.(sim.field)
        ws = MicroMagnetic.AdjointWorkspace(sim)
        KtW = apply_KT!(zeros(3N), sim, w, ws) # Kᵀ·w via the adjoint layer

        d1 = dot(Ku, w)
        d2 = dot(u, KtW)
        scale = max(norm(Ku) * norm(w), norm(u) * norm(KtW))
        @test isapprox(d1, d2; rtol = 1e-10, atol = 1e-10 * scale)
        println(rpad("  ✓ " * name, 30),
                "  <Ku,w> = ", @sprintf("%.6e", d1),
                "   rel gap = ", @sprintf("%.2e", abs(d1 - d2) / max(abs(d1), 1e-30)))
    end
end

@testset "adjoint: full-RHS duality vs LLGJacOperator" begin
    function rhs_duality(name, inters!; alpha, nx = 6, ny = 5, nz = 2, seed = 3)
        sim = build_sim(inters!; nx = nx, ny = ny, nz = nz, mode = AbstractFloat)
        N = sim.n_total
        sim.spin .= rand_unit_m0(N; seed = seed)
        op = build_matrix(sim; matrixfree = true, alpha = alpha, gamma = 2.21e5)

        Random.seed!(seed + 100)
        x = randn(2N)
        y = randn(2N)
        jvp = zeros(2N)
        mul!(jvp, op, x)                       # Df · x   (existing operator)
        vjp = zeros(2N)
        mul!(vjp, adjoint(op), y)              # Dfᵀ · λ (new adjoint mul!)

        d1 = dot(jvp, y)
        d2 = dot(x, vjp)
        scale = max(norm(jvp) * norm(y), norm(x) * norm(vjp))
        @test isapprox(d1, d2; rtol = 1e-9, atol = 1e-9 * scale)
        println(rpad("  ✓ " * name, 40),
                "  alpha = ", @sprintf("%.2f", alpha),
                "  rel gap = ", @sprintf("%.2e", abs(d1 - d2) / max(abs(d1), 1e-30)))
        return op
    end

    op = rhs_duality("exchange", sim -> add_exch(sim, 1.3e-11); alpha = 0.0)
    rhs_duality("exchange", sim -> add_exch(sim, 1.3e-11); alpha = 0.01)
    rhs_duality("exch + anis + zeeman",
                sim -> begin
                    add_exch(sim, 1.3e-11)
                    add_anis(sim, 5e4; axis = (0.2, 0.3, 1.0))
                    add_zeeman(sim, (0, 0, 800))
                end; alpha = 0.01)
    rhs_duality("exch + bulk DMI",
                sim -> begin
                    add_exch(sim, 1.3e-11)
                    add_dmi(sim, 3e-3; type = "bulk")
                end; alpha = 0.0)
    rhs_duality("exch + interfacial DMI + zeeman",
                sim -> begin
                    add_exch(sim, 1.3e-11)
                    add_dmi(sim, 3e-3; type = "interfacial")
                    add_zeeman(sim, (0, 0, 800))
                end; alpha = 0.01)
    rhs_duality("exch + demag (direct)",
                sim -> begin
                    add_exch(sim, 1.3e-11)
                    add_demag(sim)
                end; alpha = 0.01)
    rhs_duality("exch + interlayer exch/DMI",
                sim -> begin
                    add_exch(sim, 1.3e-11)
                    add_exch_int(sim, 1e-3)
                    add_dmi_int(sim, (0.0, 0.0, 1e-3))
                end; alpha = 0.01)

    # Independent check: the adjoint mul! equals the transpose of the dense
    # symbolic Jacobian (a completely separate assembly path).
    Random.seed!(7)
    y = randn(2 * op.N)
    B = Matrix(op)
    ref = adjoint(B) * y
    vjp = zeros(2 * op.N)
    mul!(vjp, adjoint(op), y)
    @test isapprox(vjp, ref; rtol = 1e-9, atol = 1e-9 * norm(B) * norm(y))
    println("  ✓ adjoint mul! == dense Bᵀ·y  (exchange, α=0)")
end

@testset "adjoint: grad_theta duality (per-cell FD)" begin
    # Central FD of ⟨K(θ±ε·δθ)·u, w⟩ against ⟨δθ, G_θ(w)⟩ with the production
    # forward kernels as ground truth.  Scalar per-cell parameters feed all
    # three direction components (Ax = Ay = Az share one array), matching the
    # grad kernels' single-scalar-per-cell convention.
    function theta_check(name, inters!, get_inter!, perturb!, gradG!;
                         nx = 5, ny = 4, nz = 3, seed = 5, eps_abs = 1e-5)
        sim = build_sim(inters!; nx = nx, ny = ny, nz = nz, mode = Float64)
        inter = get_inter!(sim)
        N = sim.n_total
        Random.seed!(seed)
        u = randn(3N)
        w = randn(3N)
        dtheta = 0.5 .+ 0.5 * randn(N)        # per-cell perturbation pattern

        perturb!(inter, eps_abs, dtheta)       # θ + ε·δθ
        MicroMagnetic.effective_field(sim, u, 0.0)
        Kp = Float64.(sim.field)
        perturb!(inter, -2 * eps_abs, dtheta)  # θ − ε·δθ  (the perturbs are additive)
        MicroMagnetic.effective_field(sim, u, 0.0)
        Km = Float64.(sim.field)
        perturb!(inter, eps_abs, dtheta)       # restore θ

        d_fd = dot(Kp - Km, w) / (2 * eps_abs)

        G = zeros(N)
        gradG!(G, sim, inter, u, w)
        d_ad = dot(dtheta, G)

        scale = max(abs(d_fd), abs(d_ad))
        @test isapprox(d_fd, d_ad; rtol = 1e-6, atol = 1e-8 * scale)
        println(rpad("  ✓ " * name, 30),
                "  d<H,w>/de = ", @sprintf("%.6e", d_fd),
                "   rel gap = ", @sprintf("%.2e", abs(d_fd - d_ad) / max(abs(d_fd), 1e-30)))
    end

    theta_check("grad Ku", sim -> add_anis(sim, (x, y, z) -> 5e4 + 1e4 * sin(1e9 * x);
                                               axis = (0.2, 0.3, 1.0)),
                sim -> sim.interactions[1],
                (anis, eps, dtheta) -> (anis.Ku .+= eps * dtheta),
                (G, sim, anis, u, w) -> begin
                    MicroMagnetic.grad_ku_kernel!(KernelAbstractions.CPU(), 512)(
                        G, u, w, anis.axis_x, anis.axis_y, anis.axis_z, sim.mu0_Ms;
                        ndrange = sim.n_total)
                end; eps_abs = 5e-1)

    theta_check("grad Aex", sim -> add_exch(sim, (x, y, z) -> 1.3e-11 + 0.3e-11 * sin(1e9 * x)),
                sim -> sim.interactions[1],
                (exch, eps, dtheta) -> (exch.Ax .+= eps * dtheta),
                (G, sim, exch, u, w) -> begin
                    mesh = sim.mesh
                    MicroMagnetic.grad_aex_kernel!(KernelAbstractions.CPU(), (128, 4, 1))(
                        G, u, w, sim.mu0_Ms, sim.inv_ms, exch.Ax, exch.Ay, exch.Az,
                        Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                        mesh.nx, mesh.ny, mesh.nz,
                        mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                        ndrange = (mesh.nx, mesh.ny, mesh.nz))
                end; eps_abs = 1.3e-16)

    theta_check("grad D (bulk)", sim -> add_dmi(sim, (x, y, z) -> 3e-3 + 0.6e-3 * cos(1e9 * y);
                                                    type = "bulk"),
                sim -> sim.interactions[1],
                (dmi, eps, dtheta) -> (dmi.Dx .+= eps * dtheta),
                (G, sim, dmi, u, w) -> begin
                    mesh = sim.mesh
                    MicroMagnetic.grad_dmi_kernel!(KernelAbstractions.CPU(), (128, 4, 1))(
                        G, u, w, sim.mu0_Ms, sim.inv_ms, dmi.Dx, dmi.Dy, dmi.Dz,
                        Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                        mesh.nx, mesh.ny, mesh.nz,
                        mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                        ndrange = (mesh.nx, mesh.ny, mesh.nz))
                end; eps_abs = 3e-8)

    theta_check("grad D (interfacial)", sim -> add_dmi(sim, (x, y, z) -> 3e-3 + 0.5e-3 * sin(1e9 * x);
                                                           type = "interfacial"),
                sim -> sim.interactions[1],
                (dmi, eps, dtheta) -> (dmi.Dx .+= eps * dtheta),
                (G, sim, dmi, u, w) -> begin
                    mesh = sim.mesh
                    MicroMagnetic.grad_dmi_interfacial_kernel!(KernelAbstractions.CPU(), (128, 4, 1))(
                        G, u, w, sim.mu0_Ms, sim.inv_ms, dmi.Dx,
                        Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                        mesh.nx, mesh.ny, mesh.nz,
                        mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                        ndrange = (mesh.nx, mesh.ny, mesh.nz))
                end; eps_abs = 3e-8)
end

@testset "adjoint: loss gradient projection" begin
    N = 60
    Random.seed!(9)
    m = rand_unit_m0(N; seed = 9)
    m_target = rand_unit_m0(N; seed = 10)

    out = zeros(3N)
    MicroMagnetic.loss_grad_kernel!(KernelAbstractions.CPU(), 512)(out, m, m_target;
                                                                   ndrange = N)
    expected = similar(out)
    for i in 1:N
        j = 3 * (i - 1)
        gx = m[j+1] - m_target[j+1]
        gy = m[j+2] - m_target[j+2]
        gz = m[j+3] - m_target[j+3]
        d = m[j+1] * gx + m[j+2] * gy + m[j+3] * gz
        expected[j+1] = gx - d * m[j+1]
        expected[j+2] = gy - d * m[j+2]
        expected[j+3] = gz - d * m[j+3]
    end
    @test maximum(abs.(out - expected)) ≤ 1e-14 * max(norm(out), 1.0)
    # tangentiality: out ⊥ m per cell
    perp = 0.0
    for i in 1:N
        j = 3 * (i - 1)
        perp = max(perp, abs(out[j+1] * m[j+1] + out[j+2] * m[j+2] + out[j+3] * m[j+3]))
    end
    @test perp ≤ 1e-12
    println("  ✓ loss_grad_kernel!: P_t(m−m*) matches and is tangential (max |·| = ",
            @sprintf("%.2e", perp), ")")
end

MicroMagnetic.set_precision(Float64)
println("test_duality.jl done.")
