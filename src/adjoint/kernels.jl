# Adjoint-layer kernels: transpose algebra, field-operator VJPs, and grad_theta kernels.
#
# All kernels are Float64-only, KernelAbstractions-based (GPU-ready: no scalar
# indexing, @Const on read-only arrays), and mirror the production kernels in
# micro/kernels.jl line by line so the duality tests in test/adjoint/ compare
# against the exact forward operators.
#
# Transpose algebra (verified term by term with scalar triple products;
# the duality tests are the final judge):
#
#   Df[δm]  = −γ′ [ δm×H + m×(Kδm) + α( δm×(m×H) + m×(δm×H) + m×(m×Kδm) ) ]
#   Dfᵀ[λ]  = −γ′ [ H×λ + Kᵀ(λ×m) + α( (m×H)×λ + (m×λ)×H + Kᵀ(m×(m×λ)) ) ]
#
# The two Kᵀ terms share the same stencil input, so the pointwise kernel emits
# the combined input  u = −γ′( λ×m + α·m×(m×λ)  )  and the driver applies each
# interaction's Kᵀ to u exactly once (Kᵀ is linear).
#
# Field-operator transpose rule: every linear stencil operator
# in this package has K = diag(1/Ms)·Ã with Ã symmetric as a full 3N×3N matrix
# — for the DMI family the antisymmetry of the [e_d]× cross-product blocks
# cancels the ±bond sign flip of the central-difference stencil, so Kᵀ = K
# modulo the Ms weights.  Hence
#
#   Kᵀ w = Ms ⊙ Op(w ⊙ 1/Ms)        (uniform Ms reduces to Op(w)),
#
# which the reuse-launcher operators (exchange/anisotropy/demag) get via the
# sandwich in linearize.jl, and the dedicated DMI/interlayer VJP kernels below
# encode directly: they are the forward kernels with (i) the per-bond 1/Ms
# weight taken from the bond SOURCE cell instead of the target cell and (ii)
# no final 1/Ms scaling of the output.  A degenerate n=1 periodic direction
# produces paired ± self-bonds whose antisymmetric blocks cancel identically
# in both the forward operator and these mirrors.

# ----------------------------------------------------------------------------------
# Pointwise helpers
# ----------------------------------------------------------------------------------

# out = a × b (pointwise per cell); used for m×H baselines.
@kernel function cross_field_kernel!(out, @Const(a), @Const(b))
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds ax, ay, az = a[j + 1], a[j + 2], a[j + 3]
    @inbounds bx, by, bz = b[j + 1], b[j + 2], b[j + 3]
    @inbounds out[j + 1] = cross_x(ax, ay, az, bx, by, bz)
    @inbounds out[j + 2] = cross_y(ax, ay, az, bx, by, bz)
    @inbounds out[j + 3] = cross_z(ax, ay, az, bx, by, bz)
end

"""
    llg_rhs_vjp_kernel!(out, u, m, H, mxH, lam, alpha, gamma_eff)

Pointwise part of the LLG right-hand-side VJP:

    out = γ_eff · [ H×λ + α( (m×H)×λ + (m×λ)×H ) ]          (γ_eff = −γ/(1+α²))
    u   = γ_eff · ( λ×m + α·m×(m×λ) )                       (combined Kᵀ stencil input)

`m`, `H`, `mxH = m×H`, `lam` are 3N interleaved; `out` and `u` are 3N outputs.
Vacuum cells (m = 0) reduce to out = γ_eff·H×λ, which is exactly the transpose
of the operator's δm×H rows there.
"""
@kernel function llg_rhs_vjp_kernel!(out, u, @Const(m), @Const(H), @Const(mxH),
                                     @Const(lam), alpha::T, gamma_eff::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds mx, my, mz = m[j + 1], m[j + 2], m[j + 3]
    @inbounds lx, ly, lz = lam[j + 1], lam[j + 2], lam[j + 3]
    @inbounds hx, hy, hz = H[j + 1], H[j + 2], H[j + 3]

    # w1 = λ×m ; w2 = m×(λ×m) = −(m×(m×λ))  (the Kᵀ input needs m×(m×λ),
    # so the α term below enters with a minus sign)
    w1x = cross_x(lx, ly, lz, mx, my, mz)
    w1y = cross_y(lx, ly, lz, mx, my, mz)
    w1z = cross_z(lx, ly, lz, mx, my, mz)
    w2x = cross_x(mx, my, mz, w1x, w1y, w1z)
    w2y = cross_y(mx, my, mz, w1x, w1y, w1z)
    w2z = cross_z(mx, my, mz, w1x, w1y, w1z)

    @inbounds px, py, pz = mxH[j + 1], mxH[j + 2], mxH[j + 3]
    tx = cross_x(px, py, pz, lx, ly, lz)
    ty = cross_y(px, py, pz, lx, ly, lz)
    tz = cross_z(px, py, pz, lx, ly, lz)
    mlx, mly, mlz = -w1x, -w1y, -w1z          # m×λ = −(λ×m)
    sx = cross_x(mlx, mly, mlz, hx, hy, hz)
    sy = cross_y(mlx, mly, mlz, hx, hy, hz)
    sz = cross_z(mlx, mly, mlz, hx, hy, hz)

    @inbounds out[j + 1] = gamma_eff * (cross_x(hx, hy, hz, lx, ly, lz) + alpha * (tx + sx))
    @inbounds out[j + 2] = gamma_eff * (cross_y(hx, hy, hz, lx, ly, lz) + alpha * (ty + sy))
    @inbounds out[j + 3] = gamma_eff * (cross_z(hx, hy, hz, lx, ly, lz) + alpha * (tz + sz))

    @inbounds u[j + 1] = gamma_eff * (w1x - alpha * w2x)
    @inbounds u[j + 2] = gamma_eff * (w1y - alpha * w2y)
    @inbounds u[j + 3] = gamma_eff * (w1z - alpha * w2z)
end

"""
    llg_param_weight_kernel!(r, m, lam, alpha, gamma_prime)

Reduced adjoint weight for the ∂f/∂θ assembly, where **θ denotes the
design variables (Ms, Aex, Ku, D) — the inverse-design parameter vector — NOT a
spherical polar angle**: with v = λ×m,

    ∂J/∂θ = γ′ ⟨ v − α·m×v , ∂H/∂θ ⟩        (θ enters f only through H)

so `r = γ′(v − α·m×v)` is the vector the grad_theta kernels dot with ∂H/∂θ.
`r` is tangential by construction (v = λ×m ⊥ m and m×v ⊥ m) — Cartesian cross
products only, no angles, no chart.
`gamma_prime` is the positive γ/(1+α²).
"""
@kernel function llg_param_weight_kernel!(r, @Const(m), @Const(lam), alpha::T,
                                          gamma_prime::T) where {T<:AbstractFloat}
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds mx, my, mz = m[j + 1], m[j + 2], m[j + 3]
    @inbounds lx, ly, lz = lam[j + 1], lam[j + 2], lam[j + 3]
    vx = cross_x(lx, ly, lz, mx, my, mz)
    vy = cross_y(lx, ly, lz, mx, my, mz)
    vz = cross_z(lx, ly, lz, mx, my, mz)
    mvx = cross_x(mx, my, mz, vx, vy, vz)
    mvy = cross_y(mx, my, mz, vx, vy, vz)
    mvz = cross_z(mx, my, mz, vx, vy, vz)
    @inbounds r[j + 1] = gamma_prime * (vx - alpha * mvx)
    @inbounds r[j + 2] = gamma_prime * (vy - alpha * mvy)
    @inbounds r[j + 3] = gamma_prime * (vz - alpha * mvz)
end

"""
    loss_grad_kernel!(out, m, m_target)

Gradient of the squared-error loss g = ½‖m − m*‖² with the tangential projection
P_t = I − m mᵀ applied per cell (singularity handling):

    out = P_t (m − m*) = (m − m*) − m (m·(m − m*))
"""
@kernel function loss_grad_kernel!(out, @Const(m), @Const(m_target))
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds mx, my, mz = m[j + 1], m[j + 2], m[j + 3]
    @inbounds gx = m[j + 1] - m_target[j + 1]
    @inbounds gy = m[j + 2] - m_target[j + 2]
    @inbounds gz = m[j + 3] - m_target[j + 3]
    d = mx * gx + my * gy + mz * gz
    @inbounds out[j + 1] = gx - d * mx
    @inbounds out[j + 2] = gy - d * my
    @inbounds out[j + 3] = gz - d * mz
end

# ----------------------------------------------------------------------------------
# DMI VJP kernels — dedicated transpose stencils
#
# Each bond of the forward kernel contributes block  (1/Ms_target)·coef·[a]×
# acting on the neighbour spin; the transpose gathers (1/Ms_source)·coef·[a]×
# — identical cross-product shape and bond signs, with the 1/Ms weight moved
# from the target cell to the bond source cell and no final output scaling.
# ----------------------------------------------------------------------------------

# CPU (ngbs-table) variant of the bulk-DMI VJP; mirrors bulkdmi_ngbs_kernel!.
@kernel function dmi_bulk_vjp_ngbs_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                           Dxs, Dys, Dzs, dx::T, dy::T, dz::T,
                                           @Const(ngbs), tfac::T) where {T<:AbstractFloat}
    I = @index(Global)
    i = 3 * I - 2
    @inbounds Ms_local = mu0_Ms[I]

    axes = (T(1 / dx), T(-1 / dx), T(1 / dy), T(-1 / dy), T(1 / dz), T(-1 / dz))

    if Ms_local == T(0)
        @inbounds h[i] = T(0)
        @inbounds h[i+1] = T(0)
        @inbounds h[i+2] = T(0)
    else
        @inbounds Dx_I = tfac * Dxs[I]
        @inbounds Dy_I = tfac * Dys[I]
        @inbounds Dz_I = tfac * Dzs[I]

        fx = T(0)
        fy = T(0)
        fz = T(0)

        # ---- (±x): block = (1/Ms_nb)·D_eff·ax·[e_x]× · w_nb ----
        for j in 1:2
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                Dx_n = tfac * Dxs[id]
                if Dx_I != 0 && Dx_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dx_I * Dx_n / (Dx_I + Dx_n)
                    ax = axes[j]
                    wgt = inv_ms[id]           # source-cell weight (transpose)
                    fy += wgt * D_eff * (-ax * w[k+2])
                    fz += wgt * D_eff * ( ax * w[k+1])
                end
            end
        end

        # ---- (±y) ----
        for j in 3:4
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                Dy_n = tfac * Dys[id]
                if Dy_I != 0 && Dy_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dy_I * Dy_n / (Dy_I + Dy_n)
                    ay = axes[j]
                    wgt = inv_ms[id]
                    fx += wgt * D_eff * ( ay * w[k+2])
                    fz += wgt * D_eff * (-ay * w[k]  )
                end
            end
        end

        # ---- (±z) ----
        for j in 5:6
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                Dz_n = tfac * Dzs[id]
                if Dz_I != 0 && Dz_n != 0
                    k = 3 * id - 2
                    D_eff = 2 * Dz_I * Dz_n / (Dz_I + Dz_n)
                    az = axes[j]
                    wgt = inv_ms[id]
                    fx += wgt * D_eff * (-az * w[k+1])
                    fy += wgt * D_eff * ( az * w[k]  )
                end
            end
        end

        @inbounds h[i] = fx
        @inbounds h[i+1] = fy
        @inbounds h[i+2] = fz
    end
end

# GPU-ready (arithmetic addressing) variant of the bulk-DMI VJP; mirrors
# bulkdmi_kernel! with mu0_Ms guards in place of mat_class (mat_class == 0
# ⟺ Ms == 0) and the source-cell 1/Ms weight per bond.
@kernel function dmi_bulk_vjp_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                      Dxs, Dys, Dzs, dx::T, dy::T, dz::T, tfac::T,
                                      nx, ny, nz, px::Bool, py::Bool, pz::Bool) where {T<:AbstractFloat}
    ci, cj, ck = @index(Global, NTuple)
    lid = (ck - 1) * nx * ny + (cj - 1) * nx + ci
    i = 3 * lid - 2
    @inbounds Ms_local = mu0_Ms[lid]

    axes = (T(1 / dx), T(-1 / dx), T(1 / dy), T(-1 / dy), T(1 / dz), T(-1 / dz))

    if Ms_local == T(0)
        @inbounds h[i] = T(0)
        @inbounds h[i+1] = T(0)
        @inbounds h[i+2] = T(0)
    else
        @inbounds Dx_I = tfac * Dxs[lid]
        @inbounds Dy_I = tfac * Dys[lid]
        @inbounds Dz_I = tfac * Dzs[lid]

        fx = T(0)
        fy = T(0)
        fz = T(0)

        # ---- (±x) ----
        @inbounds begin
            id = _ngb_neg(ci, lid, nx, px, 1)
            Dx_n = tfac * Dxs[id]
            if _valid_neg(ci, px) && mu0_Ms[id] > 0 && Dx_I != 0 && Dx_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dx_I * Dx_n / (Dx_I + Dx_n)
                ax = axes[1]
                wgt = inv_ms[id]
                fy += wgt * D_eff * (-ax * w[k+2])
                fz += wgt * D_eff * ( ax * w[k+1])
            end
            id = _ngb_pos(ci, lid, nx, px, 1)
            Dx_n = tfac * Dxs[id]
            if _valid_pos(ci, nx, px) && mu0_Ms[id] > 0 && Dx_I != 0 && Dx_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dx_I * Dx_n / (Dx_I + Dx_n)
                ax = axes[2]
                wgt = inv_ms[id]
                fy += wgt * D_eff * (-ax * w[k+2])
                fz += wgt * D_eff * ( ax * w[k+1])
            end
        end

        # ---- (±y) ----
        @inbounds begin
            id = _ngb_neg(cj, lid, ny, py, nx)
            Dy_n = tfac * Dys[id]
            if _valid_neg(cj, py) && mu0_Ms[id] > 0 && Dy_I != 0 && Dy_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dy_I * Dy_n / (Dy_I + Dy_n)
                ay = axes[3]
                wgt = inv_ms[id]
                fx += wgt * D_eff * ( ay * w[k+2])
                fz += wgt * D_eff * (-ay * w[k]  )
            end
            id = _ngb_pos(cj, lid, ny, py, nx)
            Dy_n = tfac * Dys[id]
            if _valid_pos(cj, ny, py) && mu0_Ms[id] > 0 && Dy_I != 0 && Dy_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dy_I * Dy_n / (Dy_I + Dy_n)
                ay = axes[4]
                wgt = inv_ms[id]
                fx += wgt * D_eff * ( ay * w[k+2])
                fz += wgt * D_eff * (-ay * w[k]  )
            end
        end

        # ---- (±z) ----
        @inbounds begin
            id = _ngb_neg(ck, lid, nz, pz, nx * ny)
            Dz_n = tfac * Dzs[id]
            if _valid_neg(ck, pz) && mu0_Ms[id] > 0 && Dz_I != 0 && Dz_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dz_I * Dz_n / (Dz_I + Dz_n)
                az = axes[5]
                wgt = inv_ms[id]
                fx += wgt * D_eff * (-az * w[k+1])
                fy += wgt * D_eff * ( az * w[k]  )
            end
            id = _ngb_pos(ck, lid, nz, pz, nx * ny)
            Dz_n = tfac * Dzs[id]
            if _valid_pos(ck, nz, pz) && mu0_Ms[id] > 0 && Dz_I != 0 && Dz_n != 0
                k = 3 * id - 2
                D_eff = 2 * Dz_I * Dz_n / (Dz_I + Dz_n)
                az = axes[6]
                wgt = inv_ms[id]
                fx += wgt * D_eff * (-az * w[k+1])
                fy += wgt * D_eff * ( az * w[k]  )
            end
        end

        @inbounds h[i] = fx
        @inbounds h[i+1] = fy
        @inbounds h[i+2] = fz
    end
end

# CPU (ngbs-table) variant of the interfacial-DMI VJP; mirrors
# interfacial_dmi_ngbs_kernel! (4 bonds: −x, +x, −y, +y).
@kernel function dmi_interfacial_vjp_ngbs_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                                  Ds, dx::T, dy::T, dz::T, @Const(ngbs),
                                                  tfac::T) where {T<:AbstractFloat}
    I = @index(Global)
    @inbounds Ms_local = mu0_Ms[I]
    @inbounds D_I = tfac * Ds[I]

    Dd = (T(1 / dx), T(1 / dx), T(1 / dy), T(1 / dy))
    ax = (T(0), T(0), T(-1), T(1))
    ay = (T(1), T(-1), T(0), T(0))
    az = (T(0), T(0), T(0), T(0))

    i = 3 * I - 2
    if Ms_local == T(0) || D_I == T(0)
        @inbounds h[i] = T(0)
        @inbounds h[i+1] = T(0)
        @inbounds h[i+2] = T(0)
    else
        fx = T(0)
        fy = T(0)
        fz = T(0)
        for j in 1:4
            @inbounds id = ngbs[j, I]
            @inbounds if id > 0 && mu0_Ms[id] > 0
                D_nb = tfac * Ds[id]
                if D_nb != T(0)
                    k = 3 * id - 2
                    D_eff = 2 * D_I * D_nb / (D_I + D_nb)
                    wgt = inv_ms[id]           # source-cell weight (transpose)
                    coeff = wgt * D_eff * Dd[j]
                    fx += coeff * cross_x(ax[j], ay[j], az[j], w[k], w[k+1], w[k+2])
                    fy += coeff * cross_y(ax[j], ay[j], az[j], w[k], w[k+1], w[k+2])
                    fz += coeff * cross_z(ax[j], ay[j], az[j], w[k], w[k+1], w[k+2])
                end
            end
        end
        @inbounds h[i] = fx
        @inbounds h[i+1] = fy
        @inbounds h[i+2] = fz
    end
end

# GPU-ready (arithmetic addressing) variant of the interfacial-DMI VJP; mirrors
# interfacial_dmi_kernel! with mu0_Ms guards and source-cell 1/Ms weights.
@kernel function dmi_interfacial_vjp_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                             Ds, dx::T, dy::T, dz::T, tfac::T,
                                             nx, ny, nz, px::Bool, py::Bool, pz::Bool) where {T<:AbstractFloat}
    ci, cj, ck = @index(Global, NTuple)
    lid = (ck - 1) * nx * ny + (cj - 1) * nx + ci
    i = 3 * lid - 2
    @inbounds Ms_local = mu0_Ms[lid]
    @inbounds D_I = tfac * Ds[lid]

    Dd = (T(1 / dx), T(1 / dx), T(1 / dy), T(1 / dy))
    axv = (T(0), T(0), T(-1), T(1))
    ayv = (T(1), T(-1), T(0), T(0))
    azv = (T(0), T(0), T(0), T(0))

    if Ms_local == T(0) || D_I == T(0)
        @inbounds h[i] = T(0)
        @inbounds h[i+1] = T(0)
        @inbounds h[i+2] = T(0)
    else
        fx = T(0)
        fy = T(0)
        fz = T(0)
        @inbounds begin
            id = _ngb_neg(ci, lid, nx, px, 1)
            if _valid_neg(ci, px) && mu0_Ms[id] > 0
                D_nb = tfac * Ds[id]
                if D_nb != T(0)
                    k = 3 * id - 2
                    D_eff = 2 * D_I * D_nb / (D_I + D_nb)
                    coeff = inv_ms[id] * D_eff * Dd[1]
                    fx += coeff * cross_x(axv[1], ayv[1], azv[1], w[k], w[k+1], w[k+2])
                    fy += coeff * cross_y(axv[1], ayv[1], azv[1], w[k], w[k+1], w[k+2])
                    fz += coeff * cross_z(axv[1], ayv[1], azv[1], w[k], w[k+1], w[k+2])
                end
            end
            id = _ngb_pos(ci, lid, nx, px, 1)
            if _valid_pos(ci, nx, px) && mu0_Ms[id] > 0
                D_nb = tfac * Ds[id]
                if D_nb != T(0)
                    k = 3 * id - 2
                    D_eff = 2 * D_I * D_nb / (D_I + D_nb)
                    coeff = inv_ms[id] * D_eff * Dd[2]
                    fx += coeff * cross_x(axv[2], ayv[2], azv[2], w[k], w[k+1], w[k+2])
                    fy += coeff * cross_y(axv[2], ayv[2], azv[2], w[k], w[k+1], w[k+2])
                    fz += coeff * cross_z(axv[2], ayv[2], azv[2], w[k], w[k+1], w[k+2])
                end
            end
            id = _ngb_neg(cj, lid, ny, py, nx)
            if _valid_neg(cj, py) && mu0_Ms[id] > 0
                D_nb = tfac * Ds[id]
                if D_nb != T(0)
                    k = 3 * id - 2
                    D_eff = 2 * D_I * D_nb / (D_I + D_nb)
                    coeff = inv_ms[id] * D_eff * Dd[3]
                    fx += coeff * cross_x(axv[3], ayv[3], azv[3], w[k], w[k+1], w[k+2])
                    fy += coeff * cross_y(axv[3], ayv[3], azv[3], w[k], w[k+1], w[k+2])
                    fz += coeff * cross_z(axv[3], ayv[3], azv[3], w[k], w[k+1], w[k+2])
                end
            end
            id = _ngb_pos(cj, lid, ny, py, nx)
            if _valid_pos(cj, ny, py) && mu0_Ms[id] > 0
                D_nb = tfac * Ds[id]
                if D_nb != T(0)
                    k = 3 * id - 2
                    D_eff = 2 * D_I * D_nb / (D_I + D_nb)
                    coeff = inv_ms[id] * D_eff * Dd[4]
                    fx += coeff * cross_x(axv[4], ayv[4], azv[4], w[k], w[k+1], w[k+2])
                    fy += coeff * cross_y(axv[4], ayv[4], azv[4], w[k], w[k+1], w[k+2])
                    fz += coeff * cross_z(axv[4], ayv[4], azv[4], w[k], w[k+1], w[k+2])
                end
            end
        end
        @inbounds h[i] = fx
        @inbounds h[i+1] = fy
        @inbounds h[i+2] = fz
    end
end

# ----------------------------------------------------------------------------------
# Interlayer VJP kernels — layer-index / Ms-weight swap
#
# Forward (exchange): h_k1 = J/(Ms₁dz)·m_k2, h_k2 = J/(Ms₂dz)·m_k1, so
# Kᵀw: h_k1 = J/(dz·Ms₂)·w_k2, h_k2 = J/(dz·Ms₁)·w_k1 — same sign pattern with
# the two Ms weights swapped (source-side).  For the DMI form the [D]× block
# is antisymmetric, which cancels the sign flip of the transpose exactly like
# the stencil DMI above: same signs, swapped weights.
# ----------------------------------------------------------------------------------

@kernel function interlayer_exch_vjp_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                             @Const(Js), K1::Int32, K2::Int32,
                                             nx::Int32, ny::Int32, dz::T) where {T<:AbstractFloat}
    i, j = @index(Global, NTuple)

    id1 = (K1 - 1) * nx * ny + (j - 1) * nx + i
    id2 = (K2 - 1) * nx * ny + (j - 1) * nx + i
    id = (j - 1) * nx + i

    k1 = 3 * id1 - 2
    k2 = 3 * id2 - 2

    @inbounds Ms1 = mu0_Ms[id1]
    @inbounds Ms2 = mu0_Ms[id2]
    @inbounds J = Js[id]

    if Ms1 > 0 && Ms2 > 0
        # (Kᵀw)_k1 = J/(dz·Ms₂)·w_k2 ; (Kᵀw)_k2 = J/(dz·Ms₁)·w_k1
        @inbounds h[k1] = (J / dz) * inv_ms[id2] * w[k2]
        @inbounds h[k1 + 1] = (J / dz) * inv_ms[id2] * w[k2 + 1]
        @inbounds h[k1 + 2] = (J / dz) * inv_ms[id2] * w[k2 + 2]
        @inbounds h[k2] = (J / dz) * inv_ms[id1] * w[k1]
        @inbounds h[k2 + 1] = (J / dz) * inv_ms[id1] * w[k1 + 1]
        @inbounds h[k2 + 2] = (J / dz) * inv_ms[id1] * w[k1 + 2]
    end
end

@kernel function interlayer_dmi_vjp_kernel!(@Const(w), h, @Const(mu0_Ms), @Const(inv_ms),
                                            Dx::T, Dy::T, Dz::T, K1::Int32, K2::Int32,
                                            nx::Int32, ny::Int32, dz::T) where {T<:AbstractFloat}
    i, j = @index(Global, NTuple)

    id1 = (K1 - 1) * nx * ny + (j - 1) * nx + i
    id2 = (K2 - 1) * nx * ny + (j - 1) * nx + i

    k1 = 3 * id1 - 2
    k2 = 3 * id2 - 2

    @inbounds Ms1 = mu0_Ms[id1]
    @inbounds Ms2 = mu0_Ms[id2]

    if Ms1 > 0 && Ms2 > 0
        # (Kᵀw)_k1 = +(1/(dz·Ms₂))·D×w_k2 ; (Kᵀw)_k2 = −(1/(dz·Ms₁))·D×w_k1
        @inbounds h[k1] = inv_ms[id2] * cross_x(Dx, Dy, Dz, w[k2], w[k2+1], w[k2+2]) / dz
        @inbounds h[k1 + 1] = inv_ms[id2] * cross_y(Dx, Dy, Dz, w[k2], w[k2+1], w[k2+2]) / dz
        @inbounds h[k1 + 2] = inv_ms[id2] * cross_z(Dx, Dy, Dz, w[k2], w[k2+1], w[k2+2]) / dz
        @inbounds h[k2] = -inv_ms[id1] * cross_x(Dx, Dy, Dz, w[k1], w[k1+1], w[k1+2]) / dz
        @inbounds h[k2 + 1] = -inv_ms[id1] * cross_y(Dx, Dy, Dz, w[k1], w[k1+1], w[k1+2]) / dz
        @inbounds h[k2 + 2] = -inv_ms[id1] * cross_z(Dx, Dy, Dz, w[k1], w[k1+1], w[k1+2]) / dz
    end
end

# ----------------------------------------------------------------------------------
# grad_theta kernels: (∂H/∂θ)ᵀ applied to an adjoint vector w.
#
# The per-cell derivatives include the full harmonic-mean chain
# ∂(2ab/(a+b))/∂a = 2b²/(a+b)²  — the honest per-site gradient for spatial
# parameter maps.  Summed over all cells the kernels reproduce the global
# (uniform-parameter) derivative exactly; per cell they are the gradient for
# map-valued design variables (not yet implemented).  The pair-table fast path's θ-chain
# is a future add-on.
# ----------------------------------------------------------------------------------

"""
    grad_ku_kernel!(G, m, w, axis_x, axis_y, axis_z, mu0_Ms)

∂H/∂Ku per cell: H_i = (2Ku_i/Ms_i)(m_i·a_i)a_i is linear in Ku, so
G[i] += (2/Ms_i)(m_i·a_i)(a_i·w_i).  Accumulates into G (zero-init by caller).
"""
@kernel function grad_ku_kernel!(G, @Const(m), @Const(w), @Const(axis_x), @Const(axis_y),
                                 @Const(axis_z), @Const(mu0_Ms))
    I = @index(Global)
    j = 3 * (I - 1)
    @inbounds Ms_local = mu0_Ms[I]
    if Ms_local != 0.0
        @inbounds ax = axis_x[I]
        @inbounds ay = axis_y[I]
        @inbounds az = axis_z[I]
        sa = m[j + 1] * ax + m[j + 2] * ay + m[j + 3] * az
        wa = w[j + 1] * ax + w[j + 2] * ay + w[j + 3] * az
        @inbounds G[I] += (2.0 / Ms_local) * sa * wa
    end
end

"""
    grad_aex_kernel!(G, m, w, mu0_Ms, inv_ms, Axs, Ays, Azs, dx, dy, dz, nx, ny, nz, px, py, pz)

Per-cell exchange-stiffness gradient with the harmonic-mean chain.  A_p enters
the field through every bond incident to p: as the cell's own stiffness (row p)
and as the neighbour stiffness (rows of p's neighbours):

    G[p] = Σ_{j∈N(p)} inv_ms[p]·c_d·(2A_j²/(A_p+A_j)²)·w_p·(m_j−m_p)
         + Σ_{j∈N(p)} inv_ms[j]·c_d·(2A_j²/(A_p+A_j)²)·w_j·(m_p−m_j)

with c_d = 2/d_d² per bond direction.  Accumulates into G.
"""
@kernel function grad_aex_kernel!(G, @Const(m), @Const(w), @Const(mu0_Ms), @Const(inv_ms),
                                  Axs, Ays, Azs, dx::T, dy::T, dz::T,
                                  nx, ny, nz, px::Bool, py::Bool, pz::Bool) where {T<:AbstractFloat}
    ci, cj, ck = @index(Global, NTuple)
    lid = (ck - 1) * nx * ny + (cj - 1) * nx + ci
    i = 3 * lid - 2
    @inbounds Ms_I = mu0_Ms[lid]

    nxs = T(2) / (dx * dx)
    nys = T(2) / (dy * dy)
    nzs = T(2) / (dz * dz)

    if Ms_I != T(0)
        @inbounds im_I = inv_ms[lid]
        @inbounds Ax_I = Axs[lid]
        @inbounds Ay_I = Ays[lid]
        @inbounds Az_I = Azs[lid]
        g = T(0)
        # ---- (±x) ----
        @inbounds begin
            id = _ngb_neg(ci, lid, nx, px, 1)
            Ax_n = Axs[id]
            if _valid_neg(ci, px) && mu0_Ms[id] > 0 && Ax_I != 0 && Ax_n != 0
                k = 3 * id - 2
                d_own = 2 * Ax_n^2 / (Ax_I + Ax_n)^2
                d_nb = 2 * Ax_n^2 / (Ax_I + Ax_n)^2
                g += im_I * nxs * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nxs * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
            id = _ngb_pos(ci, lid, nx, px, 1)
            Ax_n = Axs[id]
            if _valid_pos(ci, nx, px) && mu0_Ms[id] > 0 && Ax_I != 0 && Ax_n != 0
                k = 3 * id - 2
                d_own = 2 * Ax_n^2 / (Ax_I + Ax_n)^2
                d_nb = 2 * Ax_n^2 / (Ax_I + Ax_n)^2
                g += im_I * nxs * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nxs * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
        end
        # ---- (±y) ----
        @inbounds begin
            id = _ngb_neg(cj, lid, ny, py, nx)
            Ay_n = Ays[id]
            if _valid_neg(cj, py) && mu0_Ms[id] > 0 && Ay_I != 0 && Ay_n != 0
                k = 3 * id - 2
                d_own = 2 * Ay_n^2 / (Ay_I + Ay_n)^2
                d_nb = 2 * Ay_n^2 / (Ay_I + Ay_n)^2
                g += im_I * nys * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nys * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
            id = _ngb_pos(cj, lid, ny, py, nx)
            Ay_n = Ays[id]
            if _valid_pos(cj, ny, py) && mu0_Ms[id] > 0 && Ay_I != 0 && Ay_n != 0
                k = 3 * id - 2
                d_own = 2 * Ay_n^2 / (Ay_I + Ay_n)^2
                d_nb = 2 * Ay_n^2 / (Ay_I + Ay_n)^2
                g += im_I * nys * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nys * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
        end
        # ---- (±z) ----
        @inbounds begin
            id = _ngb_neg(ck, lid, nz, pz, nx * ny)
            Az_n = Azs[id]
            if _valid_neg(ck, pz) && mu0_Ms[id] > 0 && Az_I != 0 && Az_n != 0
                k = 3 * id - 2
                d_own = 2 * Az_n^2 / (Az_I + Az_n)^2
                d_nb = 2 * Az_n^2 / (Az_I + Az_n)^2
                g += im_I * nzs * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nzs * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
            id = _ngb_pos(ck, lid, nz, pz, nx * ny)
            Az_n = Azs[id]
            if _valid_pos(ck, nz, pz) && mu0_Ms[id] > 0 && Az_I != 0 && Az_n != 0
                k = 3 * id - 2
                d_own = 2 * Az_n^2 / (Az_I + Az_n)^2
                d_nb = 2 * Az_n^2 / (Az_I + Az_n)^2
                g += im_I * nzs * d_own *
                     (w[i] * (m[k] - m[i]) + w[i+1] * (m[k+1] - m[i+1]) + w[i+2] * (m[k+2] - m[i+2]))
                g += inv_ms[id] * nzs * d_nb *
                     (w[k] * (m[i] - m[k]) + w[k+1] * (m[i+1] - m[k+1]) + w[k+2] * (m[i+2] - m[k+2]))
            end
        end
        @inbounds G[lid] += g
    end
end

"""
    grad_dmi_kernel!(G, m, w, mu0_Ms, inv_ms, Dxs, Dys, Dzs, dx, dy, dz, nx, ny, nz, px, py, pz)

Per-cell bulk-DMI gradient (scalar D per cell feeding Dx = Dy = Dz) with the
harmonic-mean chain.  Each bond contributes block (1/Ms)·ax·D_eff·[e_d]× with
ax = ±1/d the signed bond factor, so

    G[p] = Σ_{j∈N(p)} inv_ms[p]·ax_j·(2D_j²/(D_p+D_j)²)·w_p·(e_d×m_j)
         − Σ_{j∈N(p)} inv_ms[j]·ax_j·(2D_j²/(D_p+D_j)²)·w_j·(e_d×m_p)

(the neighbour-row bond carries the opposite signed factor −ax_j).  A degenerate
n=1 periodic direction yields paired ± self-bonds whose contributions cancel.
Accumulates into G.
"""
@kernel function grad_dmi_kernel!(G, @Const(m), @Const(w), @Const(mu0_Ms), @Const(inv_ms),
                                  Dxs, Dys, Dzs, dx::T, dy::T, dz::T,
                                  nx, ny, nz, px::Bool, py::Bool, pz::Bool) where {T<:AbstractFloat}
    ci, cj, ck = @index(Global, NTuple)
    lid = (ck - 1) * nx * ny + (cj - 1) * nx + ci
    i = 3 * lid - 2
    @inbounds Ms_I = mu0_Ms[lid]

    axes = (T(1 / dx), T(-1 / dx), T(1 / dy), T(-1 / dy), T(1 / dz), T(-1 / dz))

    if Ms_I != T(0)
        @inbounds im_I = inv_ms[lid]
        @inbounds Dx_I = Dxs[lid]
        @inbounds Dy_I = Dys[lid]
        @inbounds Dz_I = Dzs[lid]
        g = T(0)
        # ---- (±x): cross pattern e_x×v = (0, −v_z, +v_y) ----
        @inbounds begin
            id = _ngb_neg(ci, lid, nx, px, 1)
            Dx_n = Dxs[id]
            if _valid_neg(ci, px) && mu0_Ms[id] > 0 && Dx_I != 0 && Dx_n != 0
                k = 3 * id - 2
                ax = axes[1]
                g += im_I * ax * (2 * Dx_n^2 / (Dx_I + Dx_n)^2) *
                     (-w[i+1] * m[k+2] + w[i+2] * m[k+1])
                g -= inv_ms[id] * ax * (2 * Dx_n^2 / (Dx_I + Dx_n)^2) *
                     (-w[k+1] * m[i+2] + w[k+2] * m[i+1])
            end
            id = _ngb_pos(ci, lid, nx, px, 1)
            Dx_n = Dxs[id]
            if _valid_pos(ci, nx, px) && mu0_Ms[id] > 0 && Dx_I != 0 && Dx_n != 0
                k = 3 * id - 2
                ax = axes[2]
                g += im_I * ax * (2 * Dx_n^2 / (Dx_I + Dx_n)^2) *
                     (-w[i+1] * m[k+2] + w[i+2] * m[k+1])
                g -= inv_ms[id] * ax * (2 * Dx_n^2 / (Dx_I + Dx_n)^2) *
                     (-w[k+1] * m[i+2] + w[k+2] * m[i+1])
            end
        end
        # ---- (±y): e_y×v = (+v_z, 0, −v_x) ----
        @inbounds begin
            id = _ngb_neg(cj, lid, ny, py, nx)
            Dy_n = Dys[id]
            if _valid_neg(cj, py) && mu0_Ms[id] > 0 && Dy_I != 0 && Dy_n != 0
                k = 3 * id - 2
                ay = axes[3]
                g += im_I * ay * (2 * Dy_n^2 / (Dy_I + Dy_n)^2) *
                     (w[i] * m[k+2] - w[i+2] * m[k])
                g -= inv_ms[id] * ay * (2 * Dy_n^2 / (Dy_I + Dy_n)^2) *
                     (w[k] * m[i+2] - w[k+2] * m[i])
            end
            id = _ngb_pos(cj, lid, ny, py, nx)
            Dy_n = Dys[id]
            if _valid_pos(cj, ny, py) && mu0_Ms[id] > 0 && Dy_I != 0 && Dy_n != 0
                k = 3 * id - 2
                ay = axes[4]
                g += im_I * ay * (2 * Dy_n^2 / (Dy_I + Dy_n)^2) *
                     (w[i] * m[k+2] - w[i+2] * m[k])
                g -= inv_ms[id] * ay * (2 * Dy_n^2 / (Dy_I + Dy_n)^2) *
                     (w[k] * m[i+2] - w[k+2] * m[i])
            end
        end
        # ---- (±z): e_z×v = (−v_y, +v_x, 0) ----
        @inbounds begin
            id = _ngb_neg(ck, lid, nz, pz, nx * ny)
            Dz_n = Dzs[id]
            if _valid_neg(ck, pz) && mu0_Ms[id] > 0 && Dz_I != 0 && Dz_n != 0
                k = 3 * id - 2
                az = axes[5]
                g += im_I * az * (2 * Dz_n^2 / (Dz_I + Dz_n)^2) *
                     (-w[i] * m[k+1] + w[i+1] * m[k])
                g -= inv_ms[id] * az * (2 * Dz_n^2 / (Dz_I + Dz_n)^2) *
                     (-w[k] * m[i+1] + w[k+1] * m[i])
            end
            id = _ngb_pos(ck, lid, nz, pz, nx * ny)
            Dz_n = Dzs[id]
            if _valid_pos(ck, nz, pz) && mu0_Ms[id] > 0 && Dz_I != 0 && Dz_n != 0
                k = 3 * id - 2
                az = axes[6]
                g += im_I * az * (2 * Dz_n^2 / (Dz_I + Dz_n)^2) *
                     (-w[i] * m[k+1] + w[i+1] * m[k])
                g -= inv_ms[id] * az * (2 * Dz_n^2 / (Dz_I + Dz_n)^2) *
                     (-w[k] * m[i+1] + w[k+1] * m[i])
            end
        end
        @inbounds G[lid] += g
    end
end

"""
    grad_dmi_interfacial_kernel!(G, m, w, mu0_Ms, inv_ms, Ds, dx, dy, dz, nx, ny, nz, px, py, pz)

Per-cell interfacial-DMI gradient; same structure as `grad_dmi_kernel!` with
the interfacial bond set (−x, +x, −y, +y), fixed cross axes a_j and factors
Dd_j = 1/d.  Accumulates into G.
"""
@kernel function grad_dmi_interfacial_kernel!(G, @Const(m), @Const(w), @Const(mu0_Ms),
                                              @Const(inv_ms), Ds, dx::T, dy::T, dz::T,
                                              nx, ny, nz, px::Bool, py::Bool, pz::Bool) where {T<:AbstractFloat}
    ci, cj, ck = @index(Global, NTuple)
    lid = (ck - 1) * nx * ny + (cj - 1) * nx + ci
    i = 3 * lid - 2
    @inbounds Ms_I = mu0_Ms[lid]
    @inbounds D_I = Ds[lid]

    Dd = (T(1 / dx), T(1 / dx), T(1 / dy), T(1 / dy))
    axv = (T(0), T(0), T(-1), T(1))
    ayv = (T(1), T(-1), T(0), T(0))
    azv = (T(0), T(0), T(0), T(0))

    if Ms_I != T(0) && D_I != T(0)
        @inbounds im_I = inv_ms[lid]
        g = T(0)
        # x-bonds carry axis a = (0, ay, 0): a×m = ay·(m_z, 0, −m_x)
        # y-bonds carry axis a = (ax, 0, 0): a×m = ax·(0, −m_z, m_y)
        @inbounds begin
            id = _ngb_neg(ci, lid, nx, px, 1)
            if _valid_neg(ci, px) && mu0_Ms[id] > 0 && Ds[id] != 0
                k = 3 * id - 2
                g += im_I * Dd[1] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (ayv[1] * (w[i] * m[k+2] - w[i+2] * m[k]))
                g -= inv_ms[id] * Dd[1] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (ayv[1] * (w[k] * m[i+2] - w[k+2] * m[i]))
            end
            id = _ngb_pos(ci, lid, nx, px, 1)
            if _valid_pos(ci, nx, px) && mu0_Ms[id] > 0 && Ds[id] != 0
                k = 3 * id - 2
                g += im_I * Dd[2] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (ayv[2] * (w[i] * m[k+2] - w[i+2] * m[k]))
                g -= inv_ms[id] * Dd[2] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (ayv[2] * (w[k] * m[i+2] - w[k+2] * m[i]))
            end
            id = _ngb_neg(cj, lid, ny, py, nx)
            if _valid_neg(cj, py) && mu0_Ms[id] > 0 && Ds[id] != 0
                k = 3 * id - 2
                g += im_I * Dd[3] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (axv[3] * (-w[i+1] * m[k+2] + w[i+2] * m[k+1]))
                g -= inv_ms[id] * Dd[3] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (axv[3] * (-w[k+1] * m[i+2] + w[k+2] * m[i+1]))
            end
            id = _ngb_pos(cj, lid, ny, py, nx)
            if _valid_pos(cj, ny, py) && mu0_Ms[id] > 0 && Ds[id] != 0
                k = 3 * id - 2
                g += im_I * Dd[4] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (axv[4] * (-w[i+1] * m[k+2] + w[i+2] * m[k+1]))
                g -= inv_ms[id] * Dd[4] * (2 * Ds[id]^2 / (D_I + Ds[id])^2) *
                     (axv[4] * (-w[k+1] * m[i+2] + w[k+2] * m[i+1]))
            end
        end
        @inbounds G[lid] += g
    end
end
