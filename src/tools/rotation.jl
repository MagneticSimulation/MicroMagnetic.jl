# ---------------------------------------------------------------------------
# Rigid rotation of 3-D scalar/vector fields on regular grids.
#
# Used by the LTEM tooling to emulate sample tilting: the volume (and the
# magnetization vectors attached to it) is actively rotated while the
# electron beam stays fixed along z.  The rotation is implemented as
# sequential 2-D slice rotations with bilinear interpolation, following the
# method used in MagRecon/maglab (Euler_rotation).  No external dependency
# (no PaddedViews / Images) is required.
#
# Convention
# ----------
# * Rotations are active: a point p moves to R*p, i.e. the field is sampled
#   at R^{-1} * r.
# * `euler_matrix(tx, ty, tz)` returns R = Rz(tz) * Ry(ty) * Rx(tx); the
#   data is rotated by applying Rx first, then Ry, then Rz (slice by slice).
# * Positive angles rotate counter-clockwise when looking down the
#   corresponding positive axis towards the origin.
# ---------------------------------------------------------------------------

"""
    warp2d(img, theta)

Rotate a square 2-D array by `theta` (rad) around its center using bilinear
interpolation.  Positive direction is counter-clockwise in the
(row, column) plane.  Pixels whose source falls outside `img` are zero.
"""
function warp2d(img::AbstractMatrix{T}, theta::Real) where {T}
    n, m = size(img)
    n == m || throw(ArgumentError("warp2d expects a square array, got $(size(img))"))
    out = zeros(T, n, n)
    c, s = cos(theta), sin(theta)
    cen = (n + 1) / 2
    tol = 1e-9    # keep samples that land exactly on the far boundary
    @inbounds for j in 1:n, i in 1:n
        xn = i - cen
        yn = j - cen
        # source location: rotate the target offset by -theta
        xo = c * xn + s * yn + cen
        yo = -s * xn + c * yn + cen
        (1 - tol <= xo <= n + tol && 1 - tol <= yo <= n + tol) || continue
        i1 = clamp(floor(Int, xo), 1, n - 1)
        j1 = clamp(floor(Int, yo), 1, n - 1)
        fx = xo - i1
        fy = yo - j1
        out[i, j] = (1 - fx) * (1 - fy) * img[i1, j1] + fx * (1 - fy) * img[i1+1, j1] +
                    (1 - fx) * fy * img[i1, j1+1] + fx * fy * img[i1+1, j1+1]
    end
    return out
end

"""
    euler_matrix(tx, ty, tz)

Rotation matrix `R = Rz(tz) * Ry(ty) * Rx(tx)` for angles in rad about the
x, y and z axes respectively.  (Named differently from `rotation_matrix` in
eigen/util.jl, which builds a matrix from a direction vector.)
"""
function euler_matrix(tx::Real, ty::Real, tz::Real)
    cx, sx = cos(tx), sin(tx)
    cy, sy = cos(ty), sin(ty)
    cz, sz = cos(tz), sin(tz)
    Rx = [1 0 0; 0 cx -sx; 0 sx cx]
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy]
    Rz = [cz -sz 0; sz cz 0; 0 0 1]
    return Rz * Ry * Rx
end

"""
    rotate3d(vol, tx, ty, tz)

Actively rotate a 3-D scalar field (a cubic array) by the Euler angles
`(tx, ty, tz)` (rad about x, y, z; see `euler_matrix`).  Implemented as
sequential 2-D slice rotations: x-rotation on the (y, z) slices, then
y-rotation on the (x, z) slices, then z-rotation on the (x, y) slices.
"""
function rotate3d(vol::AbstractArray{T,3}, tx::Real, ty::Real, tz::Real) where {T}
    nx, ny, nz = size(vol)
    (nx == ny == nz) || throw(ArgumentError("rotate3d expects a cubic array, got $(size(vol)); pad first (see pad_array)"))

    cur = vol
    if !iszero(tx)
        nxt = similar(cur)
        for i in 1:nx                              # slice plane: (y, z), rows = y
            nxt[i, :, :] = warp2d(cur[i, :, :], tx)
        end
        cur = nxt
    end
    if !iszero(ty)
        nxt = similar(cur)
        for j in 1:ny                              # slice plane: (x, z), rows = x
            nxt[:, j, :] = warp2d(cur[:, j, :], -ty)
        end
        cur = nxt
    end
    if !iszero(tz)
        nxt = similar(cur)
        for k in 1:nz                              # slice plane: (x, y), rows = x
            nxt[:, :, k] = warp2d(cur[:, :, k], tz)
        end
        cur = nxt
    end
    return cur === vol ? copy(cur) : cur
end

"""
    rotate3d(v, tx, ty, tz)

Rotate a 3-D vector field `v` shaped `(3, N, N, N)`: the support of each
component is rotated as in the scalar method, and afterwards the components
are mixed with the rotation matrix (a rotated vector field must satisfy
`m_new(r) = R * m_old(R^{-1} r)`).
"""
function rotate3d(v::AbstractArray{T,4}, tx::Real, ty::Real, tz::Real) where {T}
    size(v, 1) == 3 || throw(ArgumentError("rotate3d expects a (3, N, N, N) array"))
    vx = rotate3d(v[1, :, :, :], tx, ty, tz)
    vy = rotate3d(v[2, :, :, :], tx, ty, tz)
    vz = rotate3d(v[3, :, :, :], tx, ty, tz)
    R = euler_matrix(tx, ty, tz)
    out = similar(v)
    for i in 1:3
        out[i, :, :, :] .= R[i, 1] .* vx .+ R[i, 2] .* vy .+ R[i, 3] .* vz
    end
    return out
end

# ---------------------------------------------------------------------------
# Centered padding / cropping
# ---------------------------------------------------------------------------

"""
    padding_range(n, N)

Index range of a length-`n` block placed symmetrically inside a length-`N`
array (`n <= N`).  Odd/even combinations are centered on half-pixels when
needed so that padding and cropping are exact inverses.
"""
function padding_range(n::Int, N::Int)
    n <= N || throw(ArgumentError("padding target $N is smaller than the array length $n"))
    n1 = div(N - n, 2) + 1
    if iseven(N) && isodd(n)
        n1 += 1
    end
    return n1:(n1 + n - 1)
end

"""
    pad_array(v, dims...)

Center-pad (or crop) `v` to the given size.  Missing entries are zero.
"""
function pad_array(v::AbstractArray{T,N}, dims::NTuple{N,Int}) where {T,N}
    size(v) == dims && return copy(v)
    if all(size(v, i) <= dims[i] for i in 1:N)
        out = zeros(T, dims)
        out[map(padding_range, size(v), dims)...] = v
        return out
    elseif all(size(v, i) >= dims[i] for i in 1:N)
        return v[map(padding_range, dims, size(v))...]
    end
    throw(ArgumentError("pad_array cannot mix padding and cropping: size $(size(v)) vs target $dims"))
end

"""
    vector_padding(m, N)

Center-pad a `(3, nx, ny, nz)` magnetization array to `(3, N, N, N)`.
"""
vector_padding(m::AbstractArray{T,4}, N::Int) where {T} = pad_array(m, (size(m, 1), N, N, N))

# next power of two >= n (nextpow2 is not available in all versions/scopes)
_nextpow2(n::Integer) = n <= 1 ? one(n) : one(n) << (8 * sizeof(n) - leading_zeros(n - one(n)))

"""
    rotation_grid_size(nx, ny, nz, tx, ty, tz)

Smallest power-of-two cubic size that contains the bounding box of a
`(nx, ny, nz)` sample after rotation by `(tx, ty, tz)`.
"""
function rotation_grid_size(nx::Int, ny::Int, nz::Int, tx::Real, ty::Real, tz::Real)
    lx, ly, lz = _rotated_dims(Float64(nx), Float64(ny), Float64(nz), tx, ty, tz)
    return _nextpow2(ceil(Int, max(lx, ly, lz)))
end

# bounding box (floats) of a (lx, ly, lz) block after the sequential
# rotations Rx, then Ry, then Rz
function _rotated_dims(lx::Real, ly::Real, lz::Real, tx::Real, ty::Real, tz::Real)
    if !iszero(tx)
        cx, sx = abs(cos(tx)), abs(sin(tx))
        ly, lz = ly * cx + lz * sx, ly * sx + lz * cx
    end
    if !iszero(ty)
        cy, sy = abs(cos(ty)), abs(sin(ty))
        lx, lz = lx * cy + lz * sy, lx * sy + lz * cy
    end
    if !iszero(tz)
        cz, sz = abs(cos(tz)), abs(sin(tz))
        lx, ly = lx * cz + ly * sz, lx * sz + ly * cz
    end
    return lx, ly, lz
end
