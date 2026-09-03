# ---------------------------------------------------------------------------
# Projection along the electron-beam direction (+z).
# ---------------------------------------------------------------------------

"""
    project3d(s)

Integrate a 3-D scalar field along z (the beam direction): returns
`sum(s, dims=3)` with the singleton dimension removed.
"""
function project3d(s::AbstractArray{T,3}) where {T}
    return dropdims(sum(s; dims=3); dims=3)
end

"""
    project3d(v)

Integrate a `(3, N, N, N)` vector field along z.  Returns a `(3, N, N)`
array; slice `c` is the beam-integrated component `c`.
"""
function project3d(v::AbstractArray{T,4}) where {T}
    size(v, 1) == 3 || throw(ArgumentError("project3d expects a (3, N, N, N) array"))
    out = similar(v, T, (3, size(v, 2), size(v, 3)))
    for c in 1:3
        out[c, :, :] = project3d(v[c, :, :, :])
    end
    return out
end

"""
    vector_field_projection(m, theta, axis; N=-1)

Rotate the `(3, nx, ny, nz)` vector field `m` by `theta` (degrees) about the
given axis (`"x"`, `"y"` or `"z"`) and project along z.  The array is padded
to a cubic `N x N x N` grid first (`N` defaults to a size that contains the
rotated bounding box).  Returns a `(3, N, N)` array of beam-integrated
components; vacuum (`m = 0`) contributes nothing.
"""
function vector_field_projection(m::AbstractArray{T,4}, theta::Real, axis::String;
                                 N::Int=-1) where {T}
    size(m, 1) == 3 || throw(ArgumentError("m must be shaped (3, nx, ny, nz)"))
    t = deg2rad(theta)
    tx, ty, tz = axis in ("x", "X") ? (t, 0, 0) :
                 axis in ("y", "Y") ? (0, t, 0) :
                 axis in ("z", "Z") ? (0, 0, t) :
                 throw(ArgumentError("axis must be \"x\", \"y\" or \"z\", got $axis"))
    nx, ny, nz = size(m)[2:4]
    N = N > 0 ? N : rotation_grid_size(nx, ny, nz, tx, ty, tz)
    mr = rotate3d(vector_padding(m, N), tx, ty, tz)
    return project3d(mr)
end
