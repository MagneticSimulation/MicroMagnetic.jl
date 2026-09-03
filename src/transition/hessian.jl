export HessianModes, compute_hessian_modes

# The tangent-space Hessian construction follows the systematic saddle-point
# search formulation for spin systems introduced in
# G. P. Müller et al., Phys. Rev. Lett. 121, 197202 (2018),
# https://doi.org/10.1103/PhysRevLett.121.197202.

"""
    HessianModes

Lowest eigenpairs of the Riemannian energy Hessian on the product of spin
spheres. Eigenvectors are returned on the CPU so they can be inspected and
saved independently of the simulation backend. `eigenvalues` and `residuals`
use MicroMagnetic's native energy unit, joules. `converged` applies
`residual_tolerance + relative_residual_tolerance * abs(eigenvalue)`. This
flag measures the quality of an initial Ritz vector; transition acceptance
does not depend on it alone because MMF updates the mode and the final GNEB
image is independently checked for force convergence and Hessian index one.
"""
struct HessianModes
    eigenvalues::Vector{Float64}
    eigenvectors::Matrix{Float64}
    residuals::Vector{Float64}
    converged::BitVector
end

mutable struct _TransitionBackend{F,S,W,A}
    create_sim::F
    evaluator::S
    weights::W
    active::A
    active_cpu::BitVector
end

@kernel function _transition_gradient_kernel!(gradient, @Const(field), @Const(spin),
                                              @Const(weights), @Const(active))
    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        gx = -weights[site] * field[j]
        gy = -weights[site] * field[j + 1]
        gz = -weights[site] * field[j + 2]
        radial = gx * spin[j] + gy * spin[j + 1] + gz * spin[j + 2]
        gradient[j] = gx - radial * spin[j]
        gradient[j + 1] = gy - radial * spin[j + 1]
        gradient[j + 2] = gz - radial * spin[j + 2]
    else
        gradient[j] = zero(eltype(gradient))
        gradient[j + 1] = zero(eltype(gradient))
        gradient[j + 2] = zero(eltype(gradient))
    end
end

@kernel function _transition_euclidean_gradient_kernel!(
    gradient, @Const(field), @Const(weights), @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        gradient[j] = -weights[site] * field[j]
        gradient[j + 1] = -weights[site] * field[j + 1]
        gradient[j + 2] = -weights[site] * field[j + 2]
    else
        gradient[j] = zero(eltype(gradient))
        gradient[j + 1] = zero(eltype(gradient))
        gradient[j + 2] = zero(eltype(gradient))
    end
end

@kernel function _transition_affine_hessian_kernel!(
    output, @Const(field_direction), @Const(constant_field),
    @Const(base_gradient), @Const(spin), @Const(direction),
    @Const(weights), @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        dgx = -weights[site] *
              (field_direction[j] - constant_field[j])
        dgy = -weights[site] *
              (field_direction[j + 1] - constant_field[j + 1])
        dgz = -weights[site] *
              (field_direction[j + 2] - constant_field[j + 2])
        derivative_radial =
            dgx * spin[j] + dgy * spin[j + 1] + dgz * spin[j + 2]
        base_radial =
            base_gradient[j] * spin[j] +
            base_gradient[j + 1] * spin[j + 1] +
            base_gradient[j + 2] * spin[j + 2]
        output[j] = dgx - derivative_radial * spin[j] -
                    base_radial * direction[j]
        output[j + 1] = dgy - derivative_radial * spin[j + 1] -
                        base_radial * direction[j + 1]
        output[j + 2] = dgz - derivative_radial * spin[j + 2] -
                        base_radial * direction[j + 2]
    else
        output[j] = zero(eltype(output))
        output[j + 1] = zero(eltype(output))
        output[j + 2] = zero(eltype(output))
    end
end

@kernel function _project_tangent_kernel!(vector, @Const(spin), @Const(active))
    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        radial = vector[j] * spin[j] + vector[j + 1] * spin[j + 1] +
                 vector[j + 2] * spin[j + 2]
        vector[j] -= radial * spin[j]
        vector[j + 1] -= radial * spin[j + 1]
        vector[j + 2] -= radial * spin[j + 2]
    else
        vector[j] = zero(eltype(vector))
        vector[j + 1] = zero(eltype(vector))
        vector[j + 2] = zero(eltype(vector))
    end
end

@kernel function _site_norm_kernel!(site_norms, @Const(vector), @Const(active))
    site = @index(Global)
    j = 3 * site - 2
    site_norms[site] = active[site] ?
                       sqrt(vector[j]^2 + vector[j + 1]^2 + vector[j + 2]^2) :
                       zero(eltype(site_norms))
end

@kernel function _geodesic_step_kernel!(output, @Const(spin), @Const(direction),
                                        step_size, @Const(active))
    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        dx, dy, dz = direction[j], direction[j + 1], direction[j + 2]
        tangent_norm = sqrt(dx^2 + dy^2 + dz^2)
        if tangent_norm > eps(eltype(output))
            angle = step_size * tangent_norm
            c = cos(angle)
            s = sin(angle) / tangent_norm
            x = c * spin[j] + s * dx
            y = c * spin[j + 1] + s * dy
            z = c * spin[j + 2] + s * dz
            local_norm = sqrt(x^2 + y^2 + z^2)
            output[j] = x / local_norm
            output[j + 1] = y / local_norm
            output[j + 2] = z / local_norm
        else
            output[j] = spin[j]
            output[j + 1] = spin[j + 1]
            output[j + 2] = spin[j + 2]
        end
    else
        output[j] = spin[j]
        output[j + 1] = spin[j + 1]
        output[j + 2] = spin[j + 2]
    end
end

@kernel function _transport_tangent_kernel!(output, @Const(vector), @Const(spin_from),
                                            @Const(spin_to), @Const(active))
    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        denominator = one(eltype(output)) + spin_from[j] * spin_to[j] +
                      spin_from[j + 1] * spin_to[j + 1] +
                      spin_from[j + 2] * spin_to[j + 2]
        if denominator > eltype(output)(1e-12)
            vb = vector[j] * spin_to[j] + vector[j + 1] * spin_to[j + 1] +
                 vector[j + 2] * spin_to[j + 2]
            factor = vb / denominator
            output[j] = vector[j] - factor * (spin_from[j] + spin_to[j])
            output[j + 1] = vector[j + 1] -
                            factor * (spin_from[j + 1] + spin_to[j + 1])
            output[j + 2] = vector[j + 2] -
                            factor * (spin_from[j + 2] + spin_to[j + 2])
        else
            output[j] = vector[j]
            output[j + 1] = vector[j + 1]
            output[j + 2] = vector[j + 2]
        end
    else
        output[j] = zero(eltype(output))
        output[j + 1] = zero(eltype(output))
        output[j + 2] = zero(eltype(output))
    end
end

function _transition_weights(sim)
    if hasproperty(sim, :mu_s)
        return Array(sim.mu_s)
    elseif hasproperty(sim, :mu0_Ms)
        mesh = sim.mesh
        volume = hasproperty(mesh, :volume) ? mesh.volume :
                 mesh.dx * mesh.dy * mesh.dz
        return Array(sim.mu0_Ms) .* volume
    end
    error("The simulation does not provide mu_s or mu0_Ms energy weights.")
end

function _TransitionBackend(create_sim::F) where {F}
    sim = create_sim()
    sim isa AbstractSim || error("create_sim must return an AbstractSim.")
    T = eltype(sim.spin)
    weights_cpu = T.(_transition_weights(sim))
    active_cpu = BitVector(weights_cpu .> 0)
    if hasproperty(sim, :pins)
        pins = BitVector(Array(sim.pins))
        length(pins) == length(active_cpu) ||
            error("The simulation pins and magnetic weights disagree.")
        active_cpu .&= .!pins
    end
    any(active_cpu) || error("The simulation contains no active magnetic sites.")
    weights = kernel_array(weights_cpu)
    active = kernel_array(Vector{Bool}(active_cpu))
    return _TransitionBackend{F,typeof(sim),typeof(weights),typeof(active)}(
        create_sim, sim, weights, active, active_cpu)
end

const _AffineAtomisticInteraction = Union{
    Zeeman,
    HeisenbergExchange,
    SpatialHeisenberg,
    HeisenbergDMI,
    SpatialHeisenbergDMI,
    HeisenbergCantedDMI,
    HeisenbergTubeBulkDMI,
}

function _supports_affine_transition_hessian(backend::_TransitionBackend)
    sim = backend.evaluator
    return sim isa AtomisticSim &&
           all(weight -> weight > 0, Array(sim.mu_s)) &&
           all(interaction -> interaction isa _AffineAtomisticInteraction,
               sim.interactions)
end

function _prepare_transition_hessian(backend::_TransitionBackend, spin)
    _supports_affine_transition_hessian(backend) || return nothing
    sim = backend.evaluator
    zero_spin = create_zeros(length(spin))
    _set_transition_spin!(sim, zero_spin)
    effective_field_energy(sim, sim.spin, 0.0)
    constant_field = copy(sim.field)

    _set_transition_spin!(sim, spin)
    effective_field_energy(sim, sim.spin, 0.0)
    base_gradient = create_zeros(length(spin))
    kernel! = _transition_euclidean_gradient_kernel!(
        get_backend(spin), groupsize[])
    kernel!(base_gradient, sim.field, backend.weights, backend.active;
            ndrange=length(backend.active_cpu))
    return (constant_field=constant_field, base_gradient=base_gradient)
end

function _copy_transition!(destination, source)
    prepared = get_backend(source) == CPU() ? Array(source) : source
    copyto!(destination, prepared)
    return destination
end

function _device_vector(values)
    output = create_zeros(length(values))
    _copy_transition!(output, values)
    return output
end

function _check_transition_spin(backend::_TransitionBackend, spin)
    expected = 3 * length(backend.active_cpu)
    length(spin) == expected ||
        throw(DimensionMismatch("Expected $expected spin components, got $(length(spin))."))
    all(isfinite, spin) || error("The spin state contains NaN or Inf.")
    return nothing
end

function _set_transition_spin!(sim, spin)
    _copy_transition!(sim.spin, spin)
    _copy_transition!(sim.prespin, spin)
    return sim
end

function _transition_energy(backend::_TransitionBackend, spin)
    _check_transition_spin(backend, spin)
    sim = backend.evaluator
    _set_transition_spin!(sim, spin)
    effective_field_energy(sim, sim.spin, 0.0)
    return Float64(sum(sim.energy))
end

function _transition_gradient!(gradient, backend::_TransitionBackend, spin)
    _check_transition_spin(backend, spin)
    length(gradient) == length(spin) ||
        throw(DimensionMismatch("Gradient and spin arrays must have equal length."))
    sim = backend.evaluator
    _set_transition_spin!(sim, spin)
    effective_field_energy(sim, sim.spin, 0.0)
    kernel! = _transition_gradient_kernel!(get_backend(spin), groupsize[])
    kernel!(gradient, sim.field, spin, backend.weights, backend.active;
            ndrange=length(backend.active_cpu))
    return gradient
end

function _transition_energy_gradient!(gradient, backend::_TransitionBackend, spin)
    _check_transition_spin(backend, spin)
    length(gradient) == length(spin) ||
        throw(DimensionMismatch("Gradient and spin arrays must have equal length."))
    sim = backend.evaluator
    _set_transition_spin!(sim, spin)
    effective_field_energy(sim, sim.spin, 0.0)
    kernel! = _transition_gradient_kernel!(get_backend(spin), groupsize[])
    kernel!(gradient, sim.field, spin, backend.weights, backend.active;
            ndrange=length(backend.active_cpu))
    return Float64(sum(sim.energy))
end

function _project_tangent!(vector, spin, backend::_TransitionBackend)
    length(vector) == length(spin) ||
        throw(DimensionMismatch("Vector and spin arrays must have equal length."))
    kernel! = _project_tangent_kernel!(get_backend(spin), groupsize[])
    kernel!(vector, spin, backend.active; ndrange=length(backend.active_cpu))
    return vector
end

function _maximum_site_norm(vector, backend::_TransitionBackend)
    site_norms = create_zeros(length(backend.active_cpu))
    kernel! = _site_norm_kernel!(get_backend(vector), groupsize[])
    kernel!(site_norms, vector, backend.active; ndrange=length(backend.active_cpu))
    return Float64(maximum(site_norms))
end

function _normalize_tangent!(vector, spin, backend::_TransitionBackend)
    _project_tangent!(vector, spin, backend)
    nrm = norm(vector)
    nrm > 0 || error("Cannot normalize a zero tangent vector.")
    vector ./= nrm
    return vector
end

function _scale_by_maximum_site!(vector, backend::_TransitionBackend)
    scale = _maximum_site_norm(vector, backend)
    scale > 0 || return 0.0
    vector ./= scale
    return scale
end

function _geodesic_step!(output, spin, direction, step_size,
                         backend::_TransitionBackend)
    kernel! = _geodesic_step_kernel!(get_backend(spin), groupsize[])
    kernel!(output, spin, direction, eltype(spin)(step_size), backend.active;
            ndrange=length(backend.active_cpu))
    return output
end

function _transport_tangent!(output, vector, spin_from, spin_to,
                             backend::_TransitionBackend)
    kernel! = _transport_tangent_kernel!(get_backend(spin_from), groupsize[])
    kernel!(output, vector, spin_from, spin_to, backend.active;
            ndrange=length(backend.active_cpu))
    _project_tangent!(output, spin_to, backend)
    return output
end

# ---------------------------------------------------------------------------
# 2N tangent-space basis (Spirit-aligned precision mode)
# ---------------------------------------------------------------------------
#
# For each active spin site ``m_i = (x, y, z)`` we build an orthonormal pair
# ``(e_1, e_2)`` spanning the plane orthogonal to ``m_i``. The full system
# tangent-space basis is the 3N×2N block-diagonal matrix ``B`` whose i-th 3×2
# block is ``[e_1, e_2]``. With this basis the Riemannian Hessian eigenvalue
# problem can be solved in the 2N-dimensional tangent space:
#
#     H_2N = B^T (projected Hessian) B,   v_3N = B v_2N.
#
# This is exactly what Spirit does with its ``basis_3Nx2N`` and matches the
# PRL 121, 197202 (2018) reference implementation. The default
# ``:projective`` mode avoids constructing B explicitly and works directly in
# 3N with repeated radial projections; ``:full_2N`` constructs B and runs
# Lanczos in 2N, giving identical eigenvalues but better numerical stability
# at the cost of storing the basis.

@kernel function _tangent_basis_kernel!(
    basis, @Const(spin), @Const(active))
    # basis has shape (3, 2, n_sites) laid out as a flat (6, n_sites) array
    # so that basis[:, :, i] = [e1 e2] for site i.
    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        x, y, z = spin[j], spin[j + 1], spin[j + 2]
        # Pick an arbitrary vector not parallel to m, then Gram-Schmidt.
        ax, ay, az = if abs(x) <= abs(y) && abs(x) <= abs(z)
            (one(x), zero(x), zero(x))
        elseif abs(y) <= abs(z)
            (zero(x), one(x), zero(x))
        else
            (zero(x), zero(x), one(x))
        end
        # e1 = normalize(a - (a·m) m)
        dot_am = ax * x + ay * y + az * z
        e1x = ax - dot_am * x
        e1y = ay - dot_am * y
        e1z = az - dot_am * z
        n1 = sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
        b = 6 * site - 5
        if n1 > eps(eltype(basis))
            e1x, e1y, e1z = e1x / n1, e1y / n1, e1z / n1
            # e2 = m × e1 (already unit and orthogonal to both m and e1).
            e2x = y * e1z - z * e1y
            e2y = z * e1x - x * e1z
            e2z = x * e1y - y * e1x
            # Store as (3, 2) block: columns are e1, e2.
            basis[b] = e1x
            basis[b + 1] = e1y
            basis[b + 2] = e1z
            basis[b + 3] = e2x
            basis[b + 4] = e2y
            basis[b + 5] = e2z
        else
            # Degenerate (spin is zero or NaN): write NaN sentinel so the
            # CPU-side wrapper can detect and report it. error() is
            # unsupported inside GPU kernels because it requires string
            # formatting machinery the GPU compiler cannot lower.
            for k in 0:5
                basis[b + k] = eltype(basis)(NaN)
            end
        end
    else
        b = 6 * site - 5
        for k in 0:5
            basis[b + k] = zero(eltype(basis))
        end
    end
end

# Project a 3N tangent vector down to 2N coordinates: v_2N = B^T v_3N.
@kernel function _project_to_2N_kernel!(
    output, @Const(vector), @Const(basis), @Const(active), n_sites)
    # output has length 2*n_sites. Each thread handles one (site, component).
    index = @index(Global)
    site = div(index - 1, 2) + 1
    component = mod(index - 1, 2) + 1
    if site <= n_sites && active[site]
        j = 3 * site - 2
        b = 6 * site - 5 + 3 * (component - 1)
        output[index] = vector[j] * basis[b] +
                        vector[j + 1] * basis[b + 1] +
                        vector[j + 2] * basis[b + 2]
    else
        output[index] = zero(eltype(output))
    end
end

# Lift a 2N vector back to 3N: v_3N = B v_2N. The result is automatically
# tangent to the spin sphere at every active site.
@kernel function _lift_from_2N_kernel!(
    output, @Const(vector_2N), @Const(basis), @Const(active), n_sites)
    site = @index(Global)
    j = 3 * site - 2
    if site <= n_sites && active[site]
        b = 6 * site - 5
        c1 = vector_2N[2 * site - 1]
        c2 = vector_2N[2 * site]
        output[j] = c1 * basis[b] + c2 * basis[b + 3]
        output[j + 1] = c1 * basis[b + 1] + c2 * basis[b + 4]
        output[j + 2] = c1 * basis[b + 2] + c2 * basis[b + 5]
    else
        output[j] = zero(eltype(output))
        output[j + 1] = zero(eltype(output))
        output[j + 2] = zero(eltype(output))
    end
end

# Build the 2N tangent-space basis for the current spin on the active device.
# Returns a flat device array of length 6*n_sites (layout: 6 per site, first
# three entries are e1, next three are e2).
function _build_tangent_basis!(backend::_TransitionBackend, spin)
    n_sites = length(backend.active_cpu)
    basis = create_zeros(6 * n_sites)
    kernel! = _tangent_basis_kernel!(get_backend(spin), groupsize[])
    kernel!(basis, spin, backend.active; ndrange=n_sites)
    # GPU kernels cannot call error() (string formatting is unsupported on
    # device). The kernel writes NaN sentinels for degenerate sites; check
    # here on the host side. all(isfinite, ...) triggers an implicit
    # synchronize on GPU backends.
    all(isfinite, basis) ||
        error("Tangent basis degenerate: spin is zero or NaN at an active site.")
    return basis
end

# Project a 3N device vector to 2N coordinates.
function _project_vector_to_2N!(output_2N, vector_3N, basis, backend::_TransitionBackend)
    n_sites = length(backend.active_cpu)
    kernel! = _project_to_2N_kernel!(get_backend(vector_3N), groupsize[])
    kernel!(output_2N, vector_3N, basis, backend.active, n_sites;
            ndrange=2 * n_sites)
    return output_2N
end

# Lift a 2N device vector back to 3N.
function _lift_vector_from_2N!(output_3N, vector_2N, basis, backend::_TransitionBackend)
    n_sites = length(backend.active_cpu)
    kernel! = _lift_from_2N_kernel!(get_backend(vector_2N), groupsize[])
    kernel!(output_3N, vector_2N, basis, backend.active, n_sites;
            ndrange=n_sites)
    return output_3N
end

# Apply the Hessian to a 2N tangent vector: given v_2N, lift to 3N, apply the
# 3N Hessian action, project the result back to 2N. This is the
# matrix-vector product used by the 2N Lanczos iteration.
function _apply_transition_hessian_2N!(
    output_2N, backend::_TransitionBackend, spin, tangent_2N, basis,
    finite_difference_step, hessian_context, symmetry)

    lifted = create_zeros(length(spin))
    _lift_vector_from_2N!(lifted, tangent_2N, basis, backend)
    if !isnothing(symmetry)
        source = copy(lifted)
        _project_transition_vector!(lifted, source, spin, symmetry, backend)
    end
    applied_3N = create_zeros(length(spin))
    _apply_constrained_transition_hessian!(
        applied_3N, backend, spin, lifted, finite_difference_step,
        hessian_context, symmetry)
    _project_vector_to_2N!(output_2N, applied_3N, basis, backend)
    return output_2N
end

# Apply the (constrained) Riemannian Hessian to a 3N tangent vector. Uses the
# exact affine Hessian-vector action when the backend supports it; otherwise
# falls back to centered geodesic finite differences of the projected gradient.
function _apply_transition_hessian!(
    output, backend::_TransitionBackend, spin, tangent,
    finite_difference_step, hessian_context)

    finite_difference_step > 0 ||
        error("finite_difference_step must be positive.")
    direction = similar(spin)
    copyto!(direction, tangent)
    _project_tangent!(direction, spin, backend)
    direction_scale = _maximum_site_norm(direction, backend)
    direction_scale > 0 || error("Cannot apply the Hessian to a zero tangent vector.")
    direction ./= direction_scale

    if !isnothing(hessian_context)
        sim = backend.evaluator
        _set_transition_spin!(sim, direction)
        effective_field_energy(sim, sim.spin, 0.0)
        kernel! = _transition_affine_hessian_kernel!(
            get_backend(sim.spin), groupsize[])
        kernel!(
            output, sim.field, hessian_context.constant_field,
            hessian_context.base_gradient, spin, direction,
            backend.weights, backend.active;
            ndrange=length(backend.active_cpu))
        output .*= direction_scale
        _project_tangent!(output, spin, backend)
        return output
    end

    spin_plus, spin_minus = similar(spin), similar(spin)
    _geodesic_step!(spin_plus, spin, direction, finite_difference_step, backend)
    _geodesic_step!(spin_minus, spin, direction, -finite_difference_step, backend)
    gradient_plus, gradient_minus = similar(spin), similar(spin)
    _transition_gradient!(gradient_plus, backend, spin_plus)
    _transition_gradient!(gradient_minus, backend, spin_minus)
    transported_plus, transported_minus = similar(spin), similar(spin)
    _transport_tangent!(transported_plus, gradient_plus, spin_plus, spin, backend)
    _transport_tangent!(transported_minus, gradient_minus, spin_minus, spin, backend)
    @. output = direction_scale * (transported_plus - transported_minus) /
                (2 * finite_difference_step)
    _project_tangent!(output, spin, backend)
    return output
end

function _apply_constrained_transition_hessian!(
    output, backend::_TransitionBackend, spin, tangent,
    finite_difference_step, hessian_context, symmetry)

    _apply_transition_hessian!(
        output, backend, spin, tangent, finite_difference_step,
        hessian_context)
    if !isnothing(symmetry)
        source = copy(output)
        _project_transition_vector!(
            output, source, spin, symmetry, backend)
    end
    return output
end

function _orthogonalize!(vector, basis)
    for _ in 1:2, q in basis
        vector .-= dot(q, vector) .* q
    end
    return vector
end

function _orthogonalize_columns!(vector, matrix, columns)
    for _ in 1:2, column in 1:columns
        q = view(matrix, :, column)
        vector .-= dot(q, vector) .* q
    end
    return vector
end

function _hessian_residual_threshold(absolute_tolerance, relative_tolerance,
                                     spectral_scale)
    absolute_tolerance >= 0 ||
        error("residual_tolerance must be nonnegative.")
    relative_tolerance >= 0 ||
        error("relative_residual_tolerance must be nonnegative.")
    return absolute_tolerance + relative_tolerance * spectral_scale
end

function _lanczos_start(start, spin, backend, locked, rng, symmetry)
    q = if isnothing(start)
        _device_vector(randn(rng, eltype(spin), length(spin)))
    else
        _device_vector(start)
    end
    _project_tangent!(q, spin, backend)
    if !isnothing(symmetry)
        source = copy(q)
        _project_transition_vector!(q, source, spin, symmetry, backend)
    end
    _orthogonalize!(q, locked)
    if norm(q) <= 1e-14
        copyto!(q, randn(rng, eltype(spin), length(spin)))
        _project_tangent!(q, spin, backend)
        if !isnothing(symmetry)
            source = copy(q)
            _project_transition_vector!(q, source, spin, symmetry, backend)
        end
        _orthogonalize!(q, locked)
    end
    nrm = norm(q)
    nrm > 0 || error("Unable to construct a nonzero Lanczos start vector.")
    q ./= nrm
    return q
end

function _lowest_hessian_mode(backend, spin, locked, start;
                              finite_difference_step, krylov_dimension,
                              maximum_restarts, residual_tolerance,
                               relative_residual_tolerance,
                              breakdown_tolerance, hessian_context, rng,
                              symmetry)
    q = _lanczos_start(
        start, spin, backend, locked, rng, symmetry)
    best_value, best_residual = Inf, Inf
    best_vector = copy(q)
    n = length(spin)

    for _ in 1:maximum_restarts
        dimension = min(krylov_dimension, n - length(locked))
        basis = create_zeros(n, dimension)
        diagonal = zeros(Float64, dimension)
        off_diagonal = zeros(Float64, max(dimension - 1, 0))
        previous = create_zeros(n)
        previous_beta = 0.0
        actual_dimension = 0

        for column in 1:dimension
            actual_dimension = column
            copyto!(view(basis, :, column), q)
            z = create_zeros(n)
            _apply_constrained_transition_hessian!(
                z, backend, spin, q, finite_difference_step,
                hessian_context, symmetry)
            _orthogonalize!(z, locked)
            column > 1 && (z .-= previous_beta .* previous)
            alpha = Float64(dot(q, z))
            diagonal[column] = alpha
            z .-= alpha .* q
            _orthogonalize!(z, locked)
            _orthogonalize_columns!(z, basis, column)
            beta = Float64(norm(z))
            column < dimension && (off_diagonal[column] = beta)
            if beta <= breakdown_tolerance || column == dimension
                break
            end
            copyto!(previous, q)
            previous_beta = beta
            q .= z ./ beta
        end

        decomposition = eigen(SymTridiagonal(
            diagonal[1:actual_dimension],
            off_diagonal[1:max(actual_dimension - 1, 0)]))
        index = argmin(decomposition.values)
        coefficients = _device_vector(decomposition.vectors[:, index])
        ritz = view(basis, :, 1:actual_dimension) * coefficients
        _project_tangent!(ritz, spin, backend)
        _orthogonalize!(ritz, locked)
        ritz ./= norm(ritz)

        applied = create_zeros(n)
        _apply_constrained_transition_hessian!(
            applied, backend, spin, ritz, finite_difference_step,
            hessian_context, symmetry)
        _orthogonalize!(applied, locked)
        value = Float64(dot(ritz, applied))
        residual_vector = applied .- value .* ritz
        _orthogonalize!(residual_vector, locked)
        residual = Float64(norm(residual_vector))
        if residual < best_residual
            best_value, best_residual = value, residual
            copyto!(best_vector, ritz)
        end
        spectral_scale = abs(value)
        threshold = _hessian_residual_threshold(
            residual_tolerance, relative_residual_tolerance, spectral_scale)
        residual <= threshold && break
        copyto!(q, ritz)
        residual_norm = norm(residual_vector)
        if residual_norm > breakdown_tolerance
            q .+= 0.1 .* residual_vector ./ residual_norm
            _project_tangent!(q, spin, backend)
            _orthogonalize!(q, locked)
            q ./= norm(q)
        end
    end
    return best_value, best_vector, best_residual
end

# Build a 2N Lanczos start vector that lies in the (symmetry-constrained)
# tangent subspace. A 3N seed (or random vector) is projected to the tangent
# space, symmetry-projected when required, and finally compressed to 2N
# coordinates through the orthonormal basis B. Because B^T B = I, this is
# mathematically equivalent to the 3N start but stays in the reduced space.
function _lanczos_start_2N(start_3N, spin, basis, backend, locked_2N, rng, symmetry)
    function _symmetric_tangent_2N(seed_3N)
        q3 = _device_vector(seed_3N)
        _project_tangent!(q3, spin, backend)
        if !isnothing(symmetry)
            source = copy(q3)
            _project_transition_vector!(q3, source, spin, symmetry, backend)
        end
        q2 = create_zeros(2 * length(backend.active_cpu))
        _project_vector_to_2N!(q2, q3, basis, backend)
        return q2
    end

    q = if isnothing(start_3N)
        _symmetric_tangent_2N(randn(rng, eltype(spin), length(spin)))
    else
        _symmetric_tangent_2N(start_3N)
    end
    _orthogonalize!(q, locked_2N)
    if norm(q) <= 1e-14
        q = _symmetric_tangent_2N(randn(rng, eltype(spin), length(spin)))
        _orthogonalize!(q, locked_2N)
    end
    nrm = norm(q)
    nrm > 0 || error("Unable to construct a nonzero 2N Lanczos start vector.")
    q ./= nrm
    return q
end

# Restarted Lanczos in the 2N tangent space. Mirrors `_lowest_hessian_mode`
# but every vector lives in 2N coordinates and the Hessian action goes through
# `_apply_transition_hessian_2N!`. Because the 3N tangent basis B is
# orthonormal, the tridiagonal eigenvalues are identical to the 3N problem;
# the converged 2N Ritz vector is lifted back to 3N by the caller.
function _lowest_hessian_mode_2N(backend, spin, basis, locked_2N, start_3N;
                                 finite_difference_step, krylov_dimension,
                                 maximum_restarts, residual_tolerance,
                                 relative_residual_tolerance,
                                 breakdown_tolerance, hessian_context, rng,
                                 symmetry)
    n_sites = length(backend.active_cpu)
    dim_2N = 2 * n_sites
    q = _lanczos_start_2N(
        start_3N, spin, basis, backend, locked_2N, rng, symmetry)
    best_value, best_residual = Inf, Inf
    best_vector = copy(q)

    for _ in 1:maximum_restarts
        dimension = min(krylov_dimension, max(dim_2N - length(locked_2N), 1))
        krylov_basis = create_zeros(dim_2N, dimension)
        diagonal = zeros(Float64, dimension)
        off_diagonal = zeros(Float64, max(dimension - 1, 0))
        previous = create_zeros(dim_2N)
        previous_beta = 0.0
        actual_dimension = 0

        for column in 1:dimension
            actual_dimension = column
            copyto!(view(krylov_basis, :, column), q)
            z = create_zeros(dim_2N)
            _apply_transition_hessian_2N!(
                z, backend, spin, q, basis, finite_difference_step,
                hessian_context, symmetry)
            _orthogonalize!(z, locked_2N)
            column > 1 && (z .-= previous_beta .* previous)
            alpha = Float64(dot(q, z))
            diagonal[column] = alpha
            z .-= alpha .* q
            _orthogonalize!(z, locked_2N)
            _orthogonalize_columns!(z, krylov_basis, column)
            beta = Float64(norm(z))
            column < dimension && (off_diagonal[column] = beta)
            if beta <= breakdown_tolerance || column == dimension
                break
            end
            copyto!(previous, q)
            previous_beta = beta
            q .= z ./ beta
        end

        decomposition = eigen(SymTridiagonal(
            diagonal[1:actual_dimension],
            off_diagonal[1:max(actual_dimension - 1, 0)]))
        index = argmin(decomposition.values)
        coefficients = _device_vector(decomposition.vectors[:, index])
        ritz = view(krylov_basis, :, 1:actual_dimension) * coefficients
        _orthogonalize!(ritz, locked_2N)
        ritz ./= norm(ritz)

        applied = create_zeros(dim_2N)
        _apply_transition_hessian_2N!(
            applied, backend, spin, ritz, basis, finite_difference_step,
            hessian_context, symmetry)
        _orthogonalize!(applied, locked_2N)
        value = Float64(dot(ritz, applied))
        residual_vector = applied .- value .* ritz
        _orthogonalize!(residual_vector, locked_2N)
        residual = Float64(norm(residual_vector))
        if residual < best_residual
            best_value, best_residual = value, residual
            copyto!(best_vector, ritz)
        end
        spectral_scale = abs(value)
        threshold = _hessian_residual_threshold(
            residual_tolerance, relative_residual_tolerance, spectral_scale)
        residual <= threshold && break
        copyto!(q, ritz)
        residual_norm = norm(residual_vector)
        if residual_norm > breakdown_tolerance
            q .+= 0.1 .* residual_vector ./ residual_norm
            _orthogonalize!(q, locked_2N)
            q ./= norm(q)
        end
    end
    return best_value, best_vector, best_residual
end

function _compute_hessian_modes(backend::_TransitionBackend, spin;
                                n_modes=6, initial_modes=nothing,
                                finite_difference_step=2e-3,
                                krylov_dimension=28, maximum_restarts=5,
                                residual_tolerance=1e-28,
                                 relative_residual_tolerance=5e-3,
                                breakdown_tolerance=1e-30,
                                random_seed=314159,
                                symmetry=nothing,
                                tangent_basis::Symbol=:full_2N)
    n_modes > 0 || error("n_modes must be positive.")
    krylov_dimension >= n_modes ||
        error("krylov_dimension must be at least n_modes.")
    tangent_basis in (:projective, :full_2N) ||
        error("tangent_basis must be :projective or :full_2N.")
    spin_device = _device_vector(spin)
    hessian_context = _prepare_transition_hessian(backend, spin_device)
    rng = MersenneTwister(random_seed)
    locked = Any[]              # 3N locked eigenvectors (projective path)
    locked_2N = Any[]           # 2N locked eigenvectors (full_2N path)
    basis_2N = tangent_basis === :full_2N ?
               _build_tangent_basis!(backend, spin_device) : nothing
    eigenvalues = zeros(Float64, n_modes)
    eigenvectors = zeros(Float64, length(spin), n_modes)
    residuals = zeros(Float64, n_modes)

    for mode in 1:n_modes
        start = if isnothing(initial_modes)
            nothing
        elseif initial_modes isa AbstractMatrix && mode <= size(initial_modes, 2)
            view(initial_modes, :, mode)
        elseif initial_modes isa AbstractVector && mode == 1
            initial_modes
        else
            nothing
        end
        if tangent_basis === :full_2N
            value, vector_2N, residual = _lowest_hessian_mode_2N(
                backend, spin_device, basis_2N, locked_2N, start;
                finite_difference_step=finite_difference_step,
                krylov_dimension=krylov_dimension,
                maximum_restarts=maximum_restarts,
                residual_tolerance=residual_tolerance,
                relative_residual_tolerance=relative_residual_tolerance,
                breakdown_tolerance=breakdown_tolerance,
                hessian_context=hessian_context,
                rng=rng, symmetry=symmetry)
            vector = create_zeros(length(spin_device))
            _lift_vector_from_2N!(vector, vector_2N, basis_2N, backend)
            _project_tangent!(vector, spin_device, backend)
            vector ./= norm(vector)
            push!(locked_2N, copy(vector_2N))
        else
            value, vector, residual = _lowest_hessian_mode(
                backend, spin_device, locked, start;
                finite_difference_step=finite_difference_step,
                krylov_dimension=krylov_dimension,
                maximum_restarts=maximum_restarts,
                residual_tolerance=residual_tolerance,
                relative_residual_tolerance=relative_residual_tolerance,
                breakdown_tolerance=breakdown_tolerance,
                hessian_context=hessian_context,
                rng=rng, symmetry=symmetry)
            push!(locked, copy(vector))
        end
        eigenvalues[mode] = value
        eigenvectors[:, mode] .= Array(vector)
        residuals[mode] = residual
    end
    order = sortperm(eigenvalues)
    thresholds = residual_tolerance .+
                 relative_residual_tolerance .* abs.(eigenvalues)
    return HessianModes(eigenvalues[order], eigenvectors[:, order],
                        residuals[order],
                        BitVector(residuals[order] .<= thresholds[order]))
end

function _evaluate_hessian_seeds(backend::_TransitionBackend, spin, seeds;
                                 finite_difference_step=2e-3,
                                 residual_tolerance=1e-28,
                                 relative_residual_tolerance=5e-3,
                                 symmetry=nothing)
    matrix = seeds isa AbstractMatrix ? seeds : reshape(seeds, :, 1)
    size(matrix, 1) == length(spin) ||
        throw(DimensionMismatch("Initial mode rows must match the spin length."))
    spin_device = _device_vector(spin)
    hessian_context = _prepare_transition_hessian(backend, spin_device)
    count = size(matrix, 2)
    values = zeros(Float64, count)
    vectors = zeros(Float64, length(spin), count)
    residuals = zeros(Float64, count)
    for index in 1:count
        vector = _device_vector(view(matrix, :, index))
        _normalize_tangent!(vector, spin_device, backend)
        if !isnothing(symmetry)
            source = copy(vector)
            _project_transition_vector!(
                vector, source, spin_device, symmetry, backend)
            norm(vector) > 0 ||
                error("A symmetry-projected Hessian seed vanished.")
            vector ./= norm(vector)
        end
        applied = create_zeros(length(spin))
        _apply_constrained_transition_hessian!(
            applied, backend, spin_device, vector, finite_difference_step,
            hessian_context, symmetry)
        value = Float64(dot(vector, applied))
        residual = Float64(norm(applied .- value .* vector))
        values[index] = value
        vectors[:, index] .= Array(vector)
        residuals[index] = residual
    end
    order = sortperm(values)
    thresholds = residual_tolerance .+
                 relative_residual_tolerance .* abs.(values)
    return HessianModes(values[order], vectors[:, order], residuals[order],
                        BitVector(residuals[order] .<= thresholds[order]))
end

"""
    compute_hessian_modes(create_sim, spin; n_modes=6, ...)

Compute the lowest modes of the Riemannian energy Hessian using restarted
Lanczos iteration. Affine atomistic Hamiltonians use an exact
Hessian-vector action; other models automatically use centered geodesic
finite differences. `create_sim` must return an equivalent fresh simulation.
The selected MicroMagnetic backend is used for field evaluations and vector
operations. `finite_difference_step` is the angular, dimensionless spin
displacement used by the fallback;
`residual_tolerance` is in joules and `relative_residual_tolerance` is
dimensionless.

`tangent_basis` selects the working space for the eigenvalue solve:

  * `:full_2N` (default): build the explicit 3N×2N orthonormal tangent basis `B` and run
    Lanczos in the 2N tangent coordinates, matching Spirit's `basis_3Nx2N`.
    The eigenvalues are identical to `:projective` but the reduced problem is
    better conditioned, which helps small systems near degeneracies.
  * `:projective`: matrix-free Lanczos in 3N with repeated radial projections,
    retained as a compatibility option for large systems.

The tangent-space construction and low-mode search follow the systematic
saddle-point search (SPS) formulation of G. P. Müller et al., *Phys. Rev.
Lett.* **121**, 197202 (2018),
https://doi.org/10.1103/PhysRevLett.121.197202. The `:full_2N` implementation
is also aligned with the open-source Spirit framework described by G. P.
Müller et al., *Phys. Rev. B* **99**, 224414 (2019),
https://doi.org/10.1103/PhysRevB.99.224414.
"""
function compute_hessian_modes(create_sim, spin;
                               n_modes::Int=6,
                               finite_difference_step::Real=2e-3,
                               krylov_dimension::Int=max(28, n_modes + 2),
                               maximum_restarts::Int=5,
                               residual_tolerance::Real=1e-28,
                               relative_residual_tolerance::Real=5e-3,
                               breakdown_tolerance::Real=1e-30,
                               random_seed::Int=314159,
                               initial_modes=nothing,
                               tangent_basis::Symbol=:full_2N)
    backend = _TransitionBackend(create_sim)
    return _compute_hessian_modes(
        backend, spin; n_modes=n_modes, initial_modes=initial_modes,
        finite_difference_step=finite_difference_step,
        krylov_dimension=krylov_dimension,
        maximum_restarts=maximum_restarts,
        residual_tolerance=residual_tolerance,
        relative_residual_tolerance=relative_residual_tolerance,
        breakdown_tolerance=breakdown_tolerance,
        random_seed=random_seed,
        tangent_basis=tangent_basis)
end
