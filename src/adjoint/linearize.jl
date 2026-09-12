# Linearization / VJP assembly.
#
# Everything here is Float64-only and CPU-first; kernels live in kernels.jl.
# Zero edits to src/micro, src/eigen, src/atomistic: the field-operator VJP
# reuses the existing `effective_field` launchers, and the adjoint of the
# matrix-free operator is added as new mul! methods for the lazy
# `Adjoint`/`Transpose` wrappers of `LLGJacOperator`.
#
# Field-operator transpose rule (see kernels.jl header): every linear stencil
# operator here has K = diag(1/Ms)·Ã with Ã symmetric, so
#     Kᵀ w = Ms ⊙ Op(w ⊙ 1/Ms),
# implemented for exchange/anisotropy/demag by the "sandwich" around their
# existing launchers, and encoded directly in the dedicated DMI/interlayer VJP
# kernels (per-bond source-cell weights).

import LinearAlgebra

export AdjointWorkspace, apply_KT!, llg_rhs_vjp!, llg_adjoint_mul!

# ----------------------------------------------------------------------------------
# Workspace
# ----------------------------------------------------------------------------------

"""
    AdjointWorkspace(sim)

Pre-allocated Float64 scratch buffers for the adjoint layer, sized from `sim`.
Holds the combined Kᵀ stencil input, the sandwich input/output buffers, the
m×H scratch and the per-cell Ms/inv_ms maps tiled to the 3N interleaved layout.
"""
mutable struct AdjointWorkspace
    u::Vector{Float64}        # combined Kᵀ stencil input  −γ′(λ×m + α·m×(m×λ))
    ums::Vector{Float64}      # u ⊙ 1/Ms — input for reuse-launcher operators
    buf::Vector{Float64}      # operator output scratch
    mxH::Vector{Float64}      # m × H scratch
    ms3::Vector{Float64}      # Ms tiled ×3 (interleaved layout)
    inv_ms3::Vector{Float64}  # 1/Ms tiled ×3 (0 at vacuum)
end

function AdjointWorkspace(sim::AbstractSim)
    n3 = 3 * sim.n_total
    ms = Float64.(Array(sim.mu0_Ms))
    inv_ms = Float64.(Array(sim.inv_ms))
    AdjointWorkspace(zeros(n3), zeros(n3), zeros(n3), zeros(n3),
                     repeat(ms; inner = 3), repeat(inv_ms; inner = 3))
end

# ----------------------------------------------------------------------------------
# Per-interaction Kᵀ application:  out3 += Kᵀ_int · w
# ----------------------------------------------------------------------------------

# Reuse-launcher operators: Kᵀ w = Ms ⊙ Op(w ⊙ 1/Ms) around the existing
# effective_field launcher (zero new stencil kernels).
function _apply_KT_one!(out3, inter::Union{Exchange,Anisotropy}, sim::AbstractSim,
                        w::Vector{Float64}, ws::AdjointWorkspace)
    ws.ums .= w .* ws.inv_ms3
    effective_field(inter, sim, ws.ums, 0.0)
    out3 .+= ws.ms3 .* inter.field
    return nothing
end

function _apply_KT_one!(out3, demag::Union{Demag,DirectDemag,DemagPBC3D}, sim::AbstractSim,
                        w::Vector{Float64}, ws::AdjointWorkspace)
    ws.ums .= w .* ws.inv_ms3
    effective_field(demag, sim, ws.ums, 0.0; output = ws.buf)
    out3 .+= ws.ms3 .* ws.buf
    return nothing
end

# Dedicated VJP kernels: source-cell weights are built in, w is fed directly.
function _apply_KT_one!(out3, dmi::DMI, sim::AbstractSim, w::Vector{Float64},
                        ws::AdjointWorkspace)
    mesh = sim.mesh
    back = get_backend(sim.spin)
    on_cpu = back isa KernelAbstractions.CPU
    tfac = Float64(dmi.ft(0.0))
    fill!(ws.buf, 0.0)
    if dmi.type === :bulk
        if on_cpu
            dmi_bulk_vjp_ngbs_kernel!(back, 512)(w, ws.buf, sim.mu0_Ms, sim.inv_ms,
                dmi.Dx, dmi.Dy, dmi.Dz, Float64(mesh.dx), Float64(mesh.dy),
                Float64(mesh.dz), mesh.ngbs, tfac; ndrange = sim.n_total)
        else
            dmi_bulk_vjp_kernel!(back, _stencil_wg(mesh))(w, ws.buf, sim.mu0_Ms, sim.inv_ms,
                dmi.Dx, dmi.Dy, dmi.Dz, Float64(mesh.dx), Float64(mesh.dy),
                Float64(mesh.dz), tfac, mesh.nx, mesh.ny, mesh.nz,
                mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                ndrange = (mesh.nx, mesh.ny, mesh.nz))
        end
    else
        if on_cpu
            dmi_interfacial_vjp_ngbs_kernel!(back, 512)(w, ws.buf, sim.mu0_Ms, sim.inv_ms,
                dmi.Dx, Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                mesh.ngbs, tfac; ndrange = sim.n_total)
        else
            dmi_interfacial_vjp_kernel!(back, _stencil_wg(mesh))(w, ws.buf, sim.mu0_Ms,
                sim.inv_ms, dmi.Dx, Float64(mesh.dx), Float64(mesh.dy), Float64(mesh.dz),
                tfac, mesh.nx, mesh.ny, mesh.nz,
                mesh.xperiodic, mesh.yperiodic, mesh.zperiodic;
                ndrange = (mesh.nx, mesh.ny, mesh.nz))
        end
    end
    out3 .+= ws.buf
    return nothing
end

function _apply_KT_one!(out3, exch::InterlayerExchange, sim::AbstractSim,
                        w::Vector{Float64}, ws::AdjointWorkspace)
    mesh = sim.mesh
    fill!(ws.buf, 0.0)
    interlayer_exch_vjp_kernel!(get_backend(sim.spin), groupsize[])(w, ws.buf, sim.mu0_Ms,
        sim.inv_ms, exch.Js, exch.k1, exch.k2, Int32(mesh.nx), Int32(mesh.ny),
        Float64(mesh.dz); ndrange = (mesh.nx, mesh.ny))
    out3 .+= ws.buf
    return nothing
end

function _apply_KT_one!(out3, dmi::InterlayerDMI, sim::AbstractSim, w::Vector{Float64},
                        ws::AdjointWorkspace)
    mesh = sim.mesh
    fill!(ws.buf, 0.0)
    interlayer_dmi_vjp_kernel!(get_backend(sim.spin), groupsize[])(w, ws.buf, sim.mu0_Ms,
        sim.inv_ms, Float64(dmi.Dx), Float64(dmi.Dy), Float64(dmi.Dz), dmi.k1, dmi.k2,
        Int32(mesh.nx), Int32(mesh.ny), Float64(mesh.dz); ndrange = (mesh.nx, mesh.ny))
    out3 .+= ws.buf
    return nothing
end

# Zeeman: ∂H/∂m = 0 — no contribution to Kᵀ.
_apply_KT_one!(out3, ::Zeeman, sim::AbstractSim, w::Vector{Float64}, ws::AdjointWorkspace) = nothing

function _apply_KT_one!(out3, inter, sim::AbstractSim, w::Vector{Float64},
                        ws::AdjointWorkspace)
    throw(ArgumentError("adjoint layer: no Kᵀ for interaction $(nameof(typeof(inter))); " *
                        "supported: Exchange, Anisotropy (uniaxial), DMI, InterlayerExchange, " *
                        "InterlayerDMI, Demag, DirectDemag, Zeeman."))
end

"""
    apply_KT!(out3, sim, w[, ws]) -> out3

Accumulate `out3 += Kᵀ·w`, where K = ∂H_eff/∂m is the total field Jacobian of
`sim` (sum over its interactions; Zeeman contributes zero).  `w` and `out3` are
3N interleaved Float64 vectors.  This is the raw adjoint of the field
operators — the duality ⟨K·u, v⟩ = ⟨u, Kᵀ·v⟩ holds for arbitrary u, v because
every supported interaction is linear in m.
"""
function apply_KT!(out3::Vector{Float64}, sim::AbstractSim, w::Vector{Float64},
                   ws::AdjointWorkspace = AdjointWorkspace(sim))
    for inter in sim.interactions
        _apply_KT_one!(out3, inter, sim, w, ws)
    end
    return out3
end

# ----------------------------------------------------------------------------------
# Full LLG right-hand-side VJP in Cartesian 3N space
# ----------------------------------------------------------------------------------

"""
    llg_rhs_vjp!(out3, sim, m0, H0, lam[, ws; alpha=0.0, gamma=2.21e5]) -> out3

Transpose of the LLG right-hand side linearisation at the state `m0` with total
baseline field `H0`, evaluated in Cartesian 3N space:

    Dfᵀ[λ] = −γ′[ H×λ + Kᵀ(λ×m) + α( (m×H)×λ + (m×λ)×H + Kᵀ(m×(m×λ)) ) ]

with γ′ = γ/(1+α²).  `H0` must be the total effective field at `m0` (all
interactions).  Sign conventions match `LLGJacOperator`'s Cartesian
linearisation exactly; the duality tests are the final judge.
"""
function llg_rhs_vjp!(out3::Vector{Float64}, sim::AbstractSim, m0::Vector{Float64},
                      H0::Vector{Float64}, lam::Vector{Float64},
                      ws::AdjointWorkspace = AdjointWorkspace(sim);
                      alpha::Real = 0.0, gamma::Real = 2.21e5)
    n_total = sim.n_total
    length(m0) == 3 * n_total && length(H0) == 3 * n_total &&
        length(lam) == 3 * n_total && length(out3) == 3 * n_total ||
        throw(DimensionMismatch("llg_rhs_vjp!: expected 3N vectors with N = $n_total"))

    back = get_backend(sim.spin)
    cross_field_kernel!(back, groupsize[])(ws.mxH, m0, H0; ndrange = n_total)
    gamma_eff = -Float64(gamma) / (1 + Float64(alpha)^2)
    llg_rhs_vjp_kernel!(back, groupsize[])(out3, ws.u, m0, H0, ws.mxH, lam,
                                           Float64(alpha), gamma_eff; ndrange = n_total)
    apply_KT!(out3, sim, ws.u, ws)
    return out3
end

# ----------------------------------------------------------------------------------
# Adjoint of the matrix-free LLGJacOperator
#
# The operator implements B = W ∘ L ∘ U with U the tangent unwrap (2N → 3N) and
# W the tangent-frame wrap (3N → 2N); W = Uᵀ (both are "apply the tangent
# frame"), so Bᵀ = W ∘ Lᵀ ∘ U.  Lᵀ is the Cartesian VJP above, built from the
# operator's own cached baselines (m0, H0_local, per-demag H0) — no field
# recomputation beyond the Kᵀ stencil passes.
# ----------------------------------------------------------------------------------

"""
    llg_adjoint_mul!(out, op, x) -> out

Apply the adjoint of the matrix-free LLG Jacobian: `out = Bᵀ·x` for
`op::LLGJacOperator` (2N tangent vectors in/out).  Also reachable as
`mul!(out, adjoint(op), x)` / `mul!(out, transpose(op), x)`, which makes
`adjoint(op)` usable with Arpack/KrylovKit (e.g. left eigenvectors, P2a).

Allocates O(N) scratch per call (a future revision will attach a cached workspace).
"""
function llg_adjoint_mul!(out::Vector{Float64}, op::LLGJacOperator, x::Vector)
    N = op.N
    length(x) == 2N && length(out) == 2N ||
        throw(DimensionMismatch("length(out)=$(length(out)), length(x)=$(length(x)), expected $(2N)"))

    T = Float64
    n3 = 3N
    ws = AdjointWorkspace(op.sim)
    lam3 = zeros(T, n3)
    out3 = zeros(T, n3)
    H0 = zeros(T, n3)

    # 1. unwrap λ: tangent frame → Cartesian (same rotation as mul! step 1)
    @inbounds for i in 1:N
        xi = 3 * (i - 1) + 1
        R = op.Rs[i]
        x1 = x[2i - 1]
        x2 = x[2i]
        lam3[xi] = R[1, 1] * x1 + R[1, 2] * x2
        lam3[xi + 1] = R[2, 1] * x1 + R[2, 2] * x2
        lam3[xi + 2] = R[3, 1] * x1 + R[3, 2] * x2
    end

    # 2. total baseline field H0 = H0_local (local + Zeeman) + Σ per-demag H0
    copyto!(H0, op.H0_local)
    for dd in op.demag_data
        H0 .+= dd.H0
    end

    # 3. Cartesian Dfᵀ: pointwise kernel + one Kᵀ pass over the combined input
    llg_rhs_vjp!(out3, op.sim, op.m0_F64, H0, lam3, ws;
                 alpha = op.alpha, gamma = op.gamma)

    # 4. wrap: project onto the tangent frame (same R_inv dots as mul! step 4;
    #    the −γ′ scaling is already inside the Cartesian VJP)
    @inbounds for i in 1:N
        xi = 3 * (i - 1) + 1
        R_inv = op.R_invs[i]
        w1, w2, w3 = out3[xi], out3[xi + 1], out3[xi + 2]
        out[2i - 1] = R_inv[1, 1] * w1 + R_inv[1, 2] * w2 + R_inv[1, 3] * w3
        out[2i] = R_inv[2, 1] * w1 + R_inv[2, 2] * w2 + R_inv[2, 3] * w3
    end
    return out
end

# Lazy-wrapper dispatch: mul!(y, A', x) and mul!(y, transpose(A), x).
function LinearAlgebra.mul!(out::Vector{Float64},
                            op::LinearAlgebra.Adjoint{<:Any,<:LLGJacOperator},
                            x::Vector)
    return llg_adjoint_mul!(out, op.parent, x)
end

function LinearAlgebra.mul!(out::Vector{Float64},
                            op::LinearAlgebra.Transpose{<:Any,<:LLGJacOperator},
                            x::Vector)
    return llg_adjoint_mul!(out, op.parent, x)
end
