export build_matrix, build_demag_matrix, is_long_range

import SparseArrays

# Interaction classification for Route Y.
#   Local  → symbolic Epsilon kernel.
#   Demag  → Float64 impulse-columns via effective_field (backend-agnostic).
is_long_range(interaction) = false
is_long_range(interaction::Demag)       = true
is_long_range(interaction::DirectDemag) = true
is_long_range(interaction::DemagFE)     = true


"""
    build_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01)

Build the linearised 2N×2N LLG Jacobian matrix: dm̄/dt = B · m̄.

# Required: precision setup

Local interactions use the symbolic `Epsilon` / `FlatTerm` dual types defined
in `eigen/util.jl`.  The pre-allocated `sim.field` / `interaction.field`
buffers must therefore be capable of holding these values, which means the
caller must switch to the `AbstractFloat` precision **before** constructing
the sim (the struct's type parameter is baked at construction time):

```julia
MicroMagnetic.set_precision(AbstractFloat)    # BEFORE Sim(...) / Sim(...)
sim = Sim(mesh);  set_Ms(sim, 8e5);  init_m0(sim, m0);  add_exchange(…); …
B   = build_matrix(sim)
```

# Route Y

- **Local** block (Exchange, Anisotropy, DMI, Zeeman, Interlayer, torques):
  single O(N) symbolic pass → all 2N columns extracted in parallel via
  `collect_terms`.
- **Demag** block (any backend): 2N Float64 `effective_field` calls on
  impulse vectors, accumulated in-place.

When `sparse=true` the dense result is returned as `SparseMatrixCSC`.
`alpha` controls the full damped-LLG linearisation (Gilbert form).
"""
function build_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01)
    @info("Building LLG Jacobian (Route Y: Epsilon local + Float64 demag, dense 2N×2N) ...")
    if eltype(sim.spin) !== AbstractFloat
        @warn("""build_matrix requires MicroMagnetic.set_precision(AbstractFloat)
                 BEFORE constructing the sim; current eltype(sim.spin) = $(eltype(sim.spin)).
                 Otherwise symbolic writes fail with MethodError: Float64(::AddExpr).""")
    end
    return _build_matrix_impl(sim; gamma=gamma, sparse=sparse, alpha=alpha)
end


function _build_matrix_impl(sim; gamma, sparse, alpha)
    N = sim.n_total
    local SymT = Union{Float64, Epsilon, AddExpr}   # AddExpr = FlatTerm alias

    m0_F64::Vector{Float64} = convert(Vector{Float64}, Array(sim.spin))
    Ms_raw = isa(sim, MicroSim) ? Array(sim.mu0_Ms)::Vector{<:AbstractFloat} :
                                  Array(sim.mu_s)::Vector{<:AbstractFloat}
    Ms_F64::Vector{Float64} = convert(Vector{Float64}, Ms_raw)

    # Pre-compute forward / inverse rotation frames once (Householder stable).
    Rs     = Vector{Matrix{Float64}}(undef, N)
    R_invs = similar(Rs)
    @inbounds for i in 1:N
        x = 3 * (i - 1) + 1
        R = rotation_matrix(m0_F64[x], m0_F64[x+1], m0_F64[x+2])
        Rs[i]     = R
        R_invs[i] = [R[1,1] R[2,1] R[3,1];
                     R[1,2] R[2,2] R[3,2];
                     R[1,3] R[2,3] R[3,3]]
    end

    # Step 1: spin = m0 + Σ_j (ε_{2j-1} t_{j,1} + ε_{2j} t_{j,2})
    spin  = zeros(AbstractFloat, 3N)     # AbstractFloat eltype: kernels do T(0) / convert(T, _)
    dm_dt = Vector{SymT}(undef, 3N)
    fill!(dm_dt, 0.0)

    for i in 1:N
        if @inbounds Ms_F64[i] == 0
            @inbounds(spin[3i-2] = 0.0); @inbounds(spin[3i-1] = 0.0); @inbounds(spin[3i] = 0.0)
            continue
        end
        x = 3 * (i - 1)
        R = Rs[i]
        v1 = Epsilon(2i - 1, 1.0)
        v2 = Epsilon(2i,     1.0)
        @inbounds begin
            spin[x+1] = R[1,1]*v1 + R[1,2]*v2 + R[1,3]
            spin[x+2] = R[2,1]*v1 + R[2,2]*v2 + R[2,3]
            spin[x+3] = R[3,1]*v1 + R[3,2]*v2 + R[3,3]
        end
    end
    for k in eachindex(spin); spin[k] = simplify(spin[k]); end

    # Step 2a: local interactions → symbolic H_eff
    @info("  Local block: symbolic ε-pass (Exchange/Anisotropy/DMI/Zeeman …) ...")
    effective_field_local(sim, spin)
    H_eff = Vector{SymT}(undef, 3N)
    @inbounds for k in eachindex(H_eff)
        H_eff[k] = simplify(convert(SymT, sim.field[k]))
    end

    # Step 3: dmdt = -γ (m × H + α m × (m × H))  — scaled by pre_damp at the end
    pre_damp = 1.0 / (1.0 + alpha*alpha)
    use_damping = abs(alpha) > eps(Float64)

    @info("  Assembling m × H_eff (symbolic) ...")
    for i in 1:N
        x = 3 * (i - 1)
        @inbounds sx, sy, sz = spin[x+1], spin[x+2], spin[x+3]
        @inbounds hx, hy, hz = H_eff[x+1], H_eff[x+2], H_eff[x+3]
        px = cross_x(sx, sy, sz, hx, hy, hz)
        py = cross_y(sx, sy, sz, hx, hy, hz)
        pz = cross_z(sx, sy, sz, hx, hy, hz)
        if use_damping
            qx = cross_x(sx, sy, sz, px, py, pz)
            qy = cross_y(sx, sy, sz, px, py, pz)
            qz = cross_z(sx, sy, sz, px, py, pz)
            dm_dt[x+1] = -gamma * (px + alpha * qx)
            dm_dt[x+2] = -gamma * (py + alpha * qy)
            dm_dt[x+3] = -gamma * (pz + alpha * qz)
        else
            dm_dt[x+1] = -gamma * px
            dm_dt[x+2] = -gamma * py
            dm_dt[x+3] = -gamma * pz
        end
    end
    for k in eachindex(dm_dt); dm_dt[k] = simplify(dm_dt[k]); end

    # Step 4: project dmdt onto 2-D tangent frame
    for i in 1:N
        x = 3 * (i - 1)
        R_inv = R_invs[i]
        @inbounds v1, v2, v3 = dm_dt[x+1], dm_dt[x+2], dm_dt[x+3]
        @inbounds begin
            spin[x+1] = R_inv[1,1]*v1 + R_inv[1,2]*v2 + R_inv[1,3]*v3
            spin[x+2] = R_inv[2,1]*v1 + R_inv[2,2]*v2 + R_inv[2,3]*v3
            spin[x+3] = R_inv[3,1]*v1 + R_inv[3,2]*v2 + R_inv[3,3]*v3
        end
    end
    for k in eachindex(spin); spin[k] = simplify(spin[k]); end

    # Step 5: collect ε-coefficients → dense B_loc
    @info("  Extracting local ε-coefficients → 2N×2N Jacobian ...")
    B = zeros(Float64, 2N, 2N)
    for i in 1:N
        x = 3 * (i - 1)
        terms1 = collect_terms(spin[x+1])
        @inbounds for j in 1:(2N)
            v = get(terms1, j, 0.0)::Float64
            abs(v) > eps(Float64) && (B[2i-1, j] = v)
        end
        terms2 = collect_terms(spin[x+2])
        @inbounds for j in 1:(2N)
            v = get(terms2, j, 0.0)::Float64
            abs(v) > eps(Float64) && (B[2i, j] = v)
        end
    end
    use_damping && (B .*= pre_damp)

    # Step 2b: demag blocks — Float64, in-place
    lr_list = [inter for inter in sim.interactions if is_long_range(inter)]
    if !isempty(lr_list)
        @info("  Demag block: Float64 impulse-columns (2N × effective_field, backend-agnostic) ...")
        for demag in lr_list
            add_demag_block!(B, demag, sim, m0_F64, Ms_F64, Rs, R_invs;
                             gamma=gamma, alpha=alpha, pre_damp=pre_damp,
                             use_damping=use_damping)
        end
    end

    sparse && @info("  Converting dense result to SparseMatrixCSC ...")
    return sparse ? SparseArrays.sparse(B) : B
end


"""
Accumulate effective fields of LOCAL interactions (skip `is_long_range`) into
`sim.field`.  Requires sim created with `set_precision(AbstractFloat)`.
"""
function effective_field_local(sim, spin)
    fill!(sim.field, 0.0)
    for interaction in sim.interactions
        is_long_range(interaction) && continue
        effective_field(interaction, sim, spin, 0.0)
        vector_add(sim.field, interaction.field)
    end
    nothing
end


"""
In-place add the demag contribution of a single `Demag` backend to `matrix`.

Uses linearised damped-LLG:

    δṁ = -γ/(1+α²) · ( m0 × (δH + α m0×δH)  +  δm × (H0 + α m0×H0) )

where δH = T·δm comes directly from `effective_field(demag, sim, dm; output=δH)`.
`Rs` and `R_invs` are precomputed rotation frames (one per cell).
"""
function add_demag_block!(matrix::Matrix{Float64}, demag, sim,
                          m0_F64::Vector{Float64}, Ms_F64::Vector{Float64},
                          Rs::Vector{Matrix{Float64}}, R_invs::Vector{Matrix{Float64}};
                          gamma=2.21e5, alpha=0.0, pre_damp=1.0, use_damping=false)
    N = sim.n_total

    H0 = zeros(Float64, 3N)
    effective_field(demag, sim, m0_F64, 0.0; output=H0)

    m0xH0 = if use_damping
        out = zeros(Float64, 3N)
        @inbounds for i in 1:N
            xi = 3*(i-1) + 1
            a, b, c = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
            d, e, f = H0[xi], H0[xi+1], H0[xi+2]
            out[xi]   = cross_x(a, b, c, d, e, f)
            out[xi+1] = cross_y(a, b, c, d, e, f)
            out[xi+2] = cross_z(a, b, c, d, e, f)
        end
        out
    else
        nothing
    end

    dm = zeros(Float64, 3N)
    dH = zeros(Float64, 3N)
    gamma_eff = -gamma * pre_damp
    alpha_eff = alpha

    @inbounds for j in 1:N
        Ms_F64[j] == 0 && continue
        xj = 3 * (j - 1) + 1
        Rj = Rs[j]

        for alpha_axis in 1:2
            c = 2 * (j - 1) + alpha_axis
            t1u = Rj[1, alpha_axis]
            t2u = Rj[2, alpha_axis]
            t3u = Rj[3, alpha_axis]
            fill!(dm, 0.0)
            dm[xj]   = t1u
            dm[xj+1] = t2u
            dm[xj+2] = t3u

            fill!(dH, 0.0)
            effective_field(demag, sim, dm, 0.0; output=dH)

            for i in 1:N
                xi = 3 * (i - 1) + 1
                m0x, m0y, m0z = m0_F64[xi], m0_F64[xi+1], m0_F64[xi+2]
                dhx, dhy, dhz = dH[xi], dH[xi+1], dH[xi+2]

                px = cross_x(m0x, m0y, m0z, dhx, dhy, dhz)
                py = cross_y(m0x, m0y, m0z, dhx, dhy, dhz)
                pz = cross_z(m0x, m0y, m0z, dhx, dhy, dhz)
                if use_damping
                    qx = cross_x(m0x, m0y, m0z, px, py, pz)
                    qy = cross_y(m0x, m0y, m0z, px, py, pz)
                    qz = cross_z(m0x, m0y, m0z, px, py, pz)
                    v1 = px + alpha_eff * qx
                    v2 = py + alpha_eff * qy
                    v3 = pz + alpha_eff * qz
                else
                    v1, v2, v3 = px, py, pz
                end

                if i == j
                    h0x, h0y, h0z = H0[xi], H0[xi+1], H0[xi+2]
                    v1 += cross_x(t1u, t2u, t3u, h0x, h0y, h0z)
                    v2 += cross_y(t1u, t2u, t3u, h0x, h0y, h0z)
                    v3 += cross_z(t1u, t2u, t3u, h0x, h0y, h0z)
                    if use_damping
                        @inbounds ax, ay, az = m0xH0[xi], m0xH0[xi+1], m0xH0[xi+2]
                        v1 += alpha_eff * cross_x(t1u, t2u, t3u, ax, ay, az)
                        v2 += alpha_eff * cross_y(t1u, t2u, t3u, ax, ay, az)
                        v3 += alpha_eff * cross_z(t1u, t2u, t3u, ax, ay, az)
                    end
                end

                R_inv = R_invs[i]
                fx = gamma_eff * (R_inv[1,1]*v1 + R_inv[1,2]*v2 + R_inv[1,3]*v3)
                fy = gamma_eff * (R_inv[2,1]*v1 + R_inv[2,2]*v2 + R_inv[2,3]*v3)
                matrix[2*i - 1, c] += fx
                matrix[2*i,     c] += fy
            end
        end
    end
    nothing
end


"""
    build_demag_matrix(demag, sim; gamma=2.21e5, sparse=false, alpha=0.01)

Demag-only 2N×2N LLG Jacobian.  Pure Float64, zero symbolic overhead;
`sparse=true` returns a `SparseMatrixCSC`.
"""
function build_demag_matrix(demag, sim; gamma=2.21e5, sparse=false, alpha=0.01)
    @info("Building demag-only LLG Jacobian (impulse-columns, dense 2N×2N) ...")
    N = sim.n_total
    m0_F64::Vector{Float64} = convert(Vector{Float64}, Array(sim.spin))
    Ms_raw = isa(sim, MicroSim) ? Array(sim.mu0_Ms)::Vector{<:AbstractFloat} :
                                  Array(sim.mu_s)::Vector{<:AbstractFloat}
    Ms_F64::Vector{Float64} = convert(Vector{Float64}, Ms_raw)

    Rs     = [rotation_matrix(m0_F64[3i-2], m0_F64[3i-1], m0_F64[3i]) for i in 1:N]
    R_invs = map(R -> [R[1,1] R[2,1] R[3,1];
                       R[1,2] R[2,2] R[3,2];
                       R[1,3] R[2,3] R[3,3]], Rs)

    pre_damp    = 1.0 / (1.0 + alpha*alpha)
    use_damping = abs(alpha) > eps(Float64)

    B = zeros(Float64, 2N, 2N)
    add_demag_block!(B, demag, sim, m0_F64, Ms_F64, Rs, R_invs;
                     gamma=gamma, alpha=alpha, pre_damp=pre_damp,
                     use_damping=use_damping)

    sparse && @info("  Converting dense result to SparseMatrixCSC ...")
    return sparse ? SparseArrays.sparse(B) : B
end
