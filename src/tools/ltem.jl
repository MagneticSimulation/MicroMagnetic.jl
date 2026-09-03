using FFTW

# ---------------------------------------------------------------------------
# LTEM (Lorentz Transmission Electron Microscopy) tooling.
#
# Physics
# -------
# For an electron beam propagating along +z (electron charge q = -e) the
# magnetic phase shift is
#
#     φ(r) = (q/ħ) ∫ A_z(r, z) dz = -(e/ħ) ∫ A_z(r, z) dz ,
#
# which, in Fourier space over the image plane (f in cycles per meter),
# depends only on the beam-integrated in-plane magnetization
# T(r) = ∫ M⊥(r, z) dz:
#
#     φ_k(f) = (i μ0 / 2Φ0) (f_y T̂_x − f_x T̂_y) / |f|²,     Φ0 = h/2e.
#
# The equality ∫B⊥dz = μ0∫M⊥dz holds in this particular combination (the
# stray-field/H contributions cancel exactly), so the kernel can be fed with
# the projected magnetization directly.  The deflection ∝ ∇φ then matches
# the Lorentz force q v × B.
#
# Tilted samples are handled by actively rotating the volume (magnetization
# vectors included) with `rotate3d` — sequential bilinear 2-D slice
# rotations, the method used in MagRecon/maglab — while the beam stays along
# z.  The projected in-plane components of the *rotated* field are then fed
# to the kernel in laboratory frame, i.e. without rotating them back into
# the sample frame (maglab applies an extra Euler rotation to the projected
# vector, which mixes in the beam-parallel component and overestimates the
# in-plane induction by 1/cos(tilt) for a uniformly in-plane magnetized
# film; the present convention reproduces the exact Lorentz deflection
# ∫B⊥dz, which is tilt-independent for that geometry).
#
# References
# ----------
# * S. K. Walton et al., "MALTS: A tool to simulate Lorentz Transmission
#   Electron Microscopy from micromagnetic simulations", J. Phys. D (2013).
# * L. Keimpema et al., "Electron Holography Image Simulation of
#   Nanoparticles" (2006) — electric (mean inner potential) phase.
# * https://github.com/MagRecon/maglab — tilt/projection method.
# ---------------------------------------------------------------------------

const magnetic_flux_quantum = pi * h_bar / c_e   # Φ0 = h/2e (Wb)

"""
    compute_magnetic_phase_direct(mx, my, dx, dy, Lz, Ms, i0, j0)

Real-space evaluation of the magnetic phase at pixel `(i0, j0)` for a film
of thickness `Lz` with in-plane unit magnetization `(mx, my)` (2-D arrays):

    φ(r) = (μ0 Ms Lz / 2Φ0) ∫ [(y_r−y') mx(r') − (x_r−x') my(r')] / |r−r'|² dx'dy'

Slow (O(N²) per evaluation point) but exact; used as the reference for
testing the FFT implementation.
"""
function compute_magnetic_phase_direct(mx, my, dx, dy, Lz, Ms, i0, j0)
    Nx, Ny = size(mx)
    phi = 0.0
    for i in 1:Nx, j in 1:Ny
        x = (i0 - i) * dx
        y = (j0 - j) * dy
        r = sqrt(x^2 + y^2)
        if r > 0.001 * min(dx, dy)
            phi += (y * mx[i, j] - x * my[i, j]) / r^2
        end
    end
    return -phi * mu_0 * Ms * Lz / (2 * magnetic_flux_quantum) * dx * dy
end

"""
    magnetic_phase_fft(tbx, tby, dx, dy; Nout=-1)

Magnetic phase (rad, beam along +z, electron charge −e) of a magnetic sheet
with beam-integrated in-plane magnetization `tbx = ∫M_x dz`, `tby = ∫M_y dz`
(A) sampled on a grid with pixel sizes `dx`, `dy` (m):

    φ(r) = -(μ0 / 2Φ0) ∫ [(y_r−y') t_x(r') − (x_r−x') t_y(r')] / |r−r'|² d²r'

The convolution is evaluated by FFT on a linearly padded grid (the sampled
real-space kernel makes the result identical to the direct real-space sum,
up to the residual periodic images) and returned on the same (centered) grid
as the input.
"""
function magnetic_phase_fft(tbx::AbstractMatrix, tby::AbstractMatrix, dx::Real, dy::Real;
                            Nout::Int=-1)
    size(tbx) == size(tby) || throw(ArgumentError("tbx and tby must have the same size"))
    nx, ny = size(tbx)
    N = Nout > 0 ? Nout : max(nx, ny)
    (nx <= N && ny <= N) || throw(ArgumentError("Nout=$N smaller than the input size $(size(tbx))"))
    N2 = 4N                    # generous padding: the 1/r kernel converges slowly
    rng = padding_range(N, N2)

    TX = fft(pad_array(Float64.(tbx), (N2, N2)))
    TY = fft(pad_array(Float64.(tby), (N2, N2)))

    # real-space kernels sampled on the padded grid in FFT (wrapped) layout:
    # zero displacement at index 1, positive displacements towards the
    # center, negative ones wrapped to the end (r = 0 excluded).  The dxdy
    # measure is folded in so that the discrete convolution matches the
    # continuum integral (and compute_magnetic_phase_direct) exactly.
    KY = zeros(N2, N2)
    KX = zeros(N2, N2)
    half = N2 ÷ 2
    @inbounds for j in 1:N2, i in 1:N2
        oi = (i - 1) <= half ? (i - 1) : (i - 1 - N2)
        oj = (j - 1) <= half ? (j - 1) : (j - 1 - N2)
        dxc = oi * dx
        dyc = oj * dy
        r2 = dxc^2 + dyc^2
        if r2 > 0
            KY[i, j] = dx * dy * dyc / r2
            KX[i, j] = dx * dy * dxc / r2
        end
    end

    phi = -(mu_0 / (2 * magnetic_flux_quantum)) .*
          real.(ifft(fft(KY) .* TX .- fft(KX) .* TY))
    return phi[rng, rng]
end

# thickness of the material projected along z (in meters), from a
# (3, nx, ny, nz) array of unit magnetization (vacuum = 0)
function _projected_thickness(m::AbstractArray{T,4}, dz::Real) where {T}
    s = abs.(m[1, :, :, :]) .+ abs.(m[2, :, :, :]) .+ abs.(m[3, :, :, :])
    return dz .* dropdims(sum(s .> 1e-12; dims=3); dims=3)
end

# beam-integrated in-plane magnetization (A) of a possibly rotated volume;
# the mu0 factor belongs to the phase kernel, not to this projection
function _induction_integral(m::AbstractArray{T,4}, Ms::Real, dz::Real) where {T}
    tbx = Ms .* dz .* project3d(m[1, :, :, :])
    tby = Ms .* dz .* project3d(m[2, :, :, :])
    return tbx, tby
end

# tilts implied by the `axis` keyword (which axis of the data points along
# the beam before the explicit tx/ty/tz tilts are applied)
function _axis_tilts(axis::Symbol)
    return axis === :z ? (0.0, 0.0, 0.0) :
           axis === :x ? (0.0, -pi / 2, 0.0) :
           axis === :y ? (pi / 2, 0.0, 0.0) :
           throw(ArgumentError("axis must be :x, :y or :z, got $axis"))
end

# apply the `axis` base rotation first, then the explicit tx/ty/tz tilts
# (two sequential rotate3d calls keep the composition order unambiguous)
function _rotate_with_axis(m, bx, by, bz, tx, ty, tz)
    mr = iszero(bx) && iszero(by) && iszero(bz) ? m : rotate3d(m, bx, by, bz)
    return iszero(tx) && iszero(ty) && iszero(tz) ? mr : rotate3d(mr, tx, ty, tz)
end

# rotation-grid size for base rotation + tilts applied in that order
function _rotation_N(nx::Int, ny::Int, nz::Int, bx, by, bz, tx, ty, tz)
    l = _rotated_dims(Float64(nx), Float64(ny), Float64(nz), bx, by, bz)
    l = _rotated_dims(ceil(l[1]), ceil(l[2]), ceil(l[3]), tx, ty, tz)
    return _nextpow2(ceil(Int, max(l...)))
end

"""
    compute_magnetic_phase(m; Ms=1/mu_0, dx=1, dy=1, dz=1, tx=0, ty=0, tz=0,
                           axis=:z, N=-1)

Magnetic phase shift (rad) for an electron beam along +z through the
magnetization array `m` shaped `(3, nx, ny, nz)` (unit vectors; vacuum = 0).

The sample is tilted by rotating the volume about x/y/z by `tx/ty/tz` (rad)
with sequential bilinear slice rotations (method of MagRecon/maglab); `axis`
declares which data axis points along the beam before those tilts are
applied (`:x` or `:y` rotates the volume accordingly).  The rotated volume is
projected along the beam and the phase is obtained from

    φ_k = (i μ0 / 2Φ0) (f_y T̂_x − f_x T̂_y) / |f|²,   T = Ms ∫ m dz.

`Ms` is the magnetization in A/m (the default `1/mu_0` gives the phase of a
unit induction); `dx, dy, dz` are the voxel sizes in meters (dimensionless
values give the phase in matching units).  `N` sets the size of the
(cubic) rotation/output grid; by default it is chosen to contain the rotated
bounding box.  The returned `N x N` phase is centered on the sample.

# Example
```julia
ovf = read_ovf("m.ovf")
phi = compute_magnetic_phase(ovf; Ms=8e5, ty=deg2rad(20))
```
"""
function compute_magnetic_phase(m::AbstractArray{<:Real,4}; Ms::Real=1 / mu_0,
                                dx::Real=1.0, dy::Real=1.0, dz::Real=1.0,
                                tx::Real=0.0, ty::Real=0.0, tz::Real=0.0,
                                axis::Symbol=:z, N::Int=-1)
    size(m, 1) == 3 || throw(ArgumentError("m must be shaped (3, nx, ny, nz)"))
    bx, by, bz = _axis_tilts(axis)
    nx, ny, nz = size(m)[2:4]
    N = N > 0 ? N : _rotation_N(nx, ny, nz, bx, by, bz, tx, ty, tz)

    if iszero(bx + by + bz + tx + ty + tz)
        tbx, tby = _induction_integral(m, Ms, dz)           # fast path, no rotation
    else
        mr = _rotate_with_axis(vector_padding(m, N), bx, by, bz, tx, ty, tz)
        tbx, tby = _induction_integral(mr, Ms, dz)
    end
    return magnetic_phase_fft(tbx, tby, dx, dy; Nout=N)
end

"""
    compute_magnetic_phase(m, theta, axis; kwargs...)

Tilt the sample by `theta` (rad) about `axis` (`"x"` or `"y"`) and compute
the magnetic phase; a convenience wrapper for the keyword-argument method.
"""
function compute_magnetic_phase(m::AbstractArray{<:Real,4}, theta::Real, axis::String;
                                kwargs...)
    axis in ("x", "X") && return compute_magnetic_phase(m; tx=theta, kwargs...)
    axis in ("y", "Y") && return compute_magnetic_phase(m; ty=theta, kwargs...)
    throw(ArgumentError("axis must be \"x\" or \"y\", got $axis"))
end

"""
    compute_magnetic_phase(ovf::OVF2; Ms=1, tx=0, ty=0, tz=0, axis=:z, N=-1)

Magnetic phase from an `OVF2` file; the voxel sizes are taken from the file.
Pass `Ms` in A/m for a physical phase in rad (the default `Ms=1` gives the
phase of a unit magnetization).
"""
function compute_magnetic_phase(ovf::OVF2; Ms::Real=1, kwargs...)
    m = reshape(ovf.data, (3, ovf.xnodes, ovf.ynodes, ovf.znodes))
    return compute_magnetic_phase(m; Ms=Ms, dx=ovf.xstepsize, dy=ovf.ystepsize,
                                  dz=ovf.zstepsize, kwargs...)
end

# ---------------------------------------------------------------------------
# Electron optics helpers
# ---------------------------------------------------------------------------

"""
    electron_wavelength(V)

Relativistic electron wavelength (m) for an accelerating voltage `V` in volts.
"""
function electron_wavelength(V::Real)
    C = 299792458.0
    E0 = m_e * C^2
    Ek = V * c_e
    P = sqrt(Ek^2 + 2 * Ek * E0) / C
    return 2 * pi * h_bar / P
end

"""
    interaction_constant(V)

Electron interaction constant σ (rad/(V·m)) for an accelerating voltage `V`
in volts; the electric phase is `φ = σ ∫ V0 dz` with `V0` the mean inner
potential.  (≈ 6.53e6 rad/(V·m) at 300 kV.)
"""
function interaction_constant(V::Real)
    lambda = electron_wavelength(V)
    E0 = m_e * 299792458.0^2
    Ek = V * c_e
    return (2 * pi / (lambda * V)) * (Ek + E0) / (Ek + 2 * E0)
end

"""
    compute_electric_phase(V, V0, Lz, beta=0)

Return `(lambda, phi_E)`: the electron wavelength (m) and the electric phase
(rad) of a slab with mean inner potential `V0` (V), thickness `Lz` (m) and
beam tilt `beta` (rad) relative to the slab normal.
"""
function compute_electric_phase(V::Real, V0::Real, Lz::Real, beta::Real=0.0)
    return electron_wavelength(V), interaction_constant(V) * V0 * Lz / cos(beta)
end

# ---------------------------------------------------------------------------
# Full LTEM imaging (phase + Fresnel defocus contrast)
# ---------------------------------------------------------------------------

"""
    LTEM(m; V=300, Ms=1e5, V0=-26, df=1600, alpha=1e-5, tx=0, ty=0, tz=0,
         axis=:z, N=-1, dx=1, dy=1, dz=1)
    LTEM(ovf::OVF2; ...)
    LTEM(fname::String; ...)

Simulate a Lorentz TEM image: tilt the sample (`tx/ty/tz` in rad, see
`compute_magnetic_phase`), compute the magnetic phase, add the electric
(mean inner potential) phase over the projected material footprint, and
propagate with the Fresnel defocus transfer function

    T(f) = exp(-i π λ df |f|²),   E(f) = exp(-(π α df |f|)²)

with `df` the defocus in µm and `alpha` the beam divergence semiangle in rad
(MALTS, Walton et al. 2013).

Returns `(phi_M, intensity)`: the unwrapped magnetic phase (rad) and the
normalized defocus image, both `N x N` and centered on the sample.  `V` is
the accelerating voltage in kV, `V0` the mean inner potential in V, `Ms` the
magnetization in A/m; `dx, dy, dz` are the voxel sizes in meters (defaults
of 1 give results in voxel units).
"""
function LTEM(m::AbstractArray{<:Real,4}; V::Real=300, Ms::Real=1e5, V0::Real=-26,
              df::Real=1600, alpha::Real=1e-5, tx::Real=0.0, ty::Real=0.0, tz::Real=0.0,
              axis::Symbol=:z, N::Int=-1, dx::Real=1.0, dy::Real=1.0, dz::Real=1.0)
    size(m, 1) == 3 || throw(ArgumentError("m must be shaped (3, nx, ny, nz)"))
    bx, by, bz = _axis_tilts(axis)
    nx, ny, nz = size(m)[2:4]
    N = N > 0 ? N : _rotation_N(nx, ny, nz, bx, by, bz, tx, ty, tz)

    if iszero(bx + by + bz + tx + ty + tz)
        tbx, tby = _induction_integral(m, Ms, dz)
        thickness = _projected_thickness(m, dz)
    else
        mr = _rotate_with_axis(vector_padding(m, N), bx, by, bz, tx, ty, tz)
        tbx, tby = _induction_integral(mr, Ms, dz)
        thickness = _projected_thickness(mr, dz)
    end

    phi_M = magnetic_phase_fft(tbx, tby, dx, dy; Nout=N)
    # thickness lives on the (nx, ny) grid for the fast path; center it on
    # the same N x N grid as the phase
    thickness = pad_array(thickness, (N, N))
    phi_E = interaction_constant(1000 * V) * V0 .* thickness
    phi = phi_M .+ phi_E

    # Fresnel defocus imaging on a linearly padded grid
    N2 = 2N
    rng = padding_range(N, N2)
    psi = zeros(ComplexF64, N2, N2)
    psi[rng, rng] .= exp.(1im .* phi)

    dfm = df * 1e-6
    kx = fftfreq(N2; d=dx)
    ky = fftfreq(N2; d=dy)
    lambda = electron_wavelength(1000 * V)
    k2 = [kx[i]^2 + ky[j]^2 for i in 1:N2, j in 1:N2]
    T = exp.(-1im .* pi .* lambda .* dfm .* k2)
    E = exp.(-(pi .* alpha .* dfm)^2 .* k2)

    img = abs2.(ifft(fft(psi) .* E .* T))
    return phi_M, img[rng, rng]
end

function LTEM(ovf::OVF2; kwargs...)
    m = reshape(ovf.data, (3, ovf.xnodes, ovf.ynodes, ovf.znodes))
    return LTEM(m; dx=ovf.xstepsize, dy=ovf.ystepsize, dz=ovf.zstepsize, kwargs...)
end

LTEM(fname::AbstractString; kwargs...) = LTEM(read_ovf(fname); kwargs...)
