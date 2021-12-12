function inbounds(x, y, Nx, Ny)
    x1, x2 = floor(Int, x), floor(Int, x) + 1
    y1, y2 = floor(Int, y), floor(Int, y) + 1

    # out of bounds
    if x1 < 1 || x2 > Nx || y1 < 1 || y2 > Ny
        return false
    end

    return true
end

function inbounds(x, y, z, Nx, Ny, Nz)
    x1, x2 = floor(Int, x), floor(Int, x) + 1
    y1, y2 = floor(Int, y), floor(Int, y) + 1
    z1, z2 = floor(Int, z), floor(Int, z) + 1

    # out of bounds
    if x1 < 1 || x2 > Nx || y1 < 1 || y2 > Ny || z1 < 1 || z2 > Nz
        return false
    end

    return true
end

function bilinear_interpolation(x::Real, y::Real, a::Array{T, 2}) where T<:AbstractFloat
    """
    Return bilinear interpolation value on point (x,y).

    ---------------------
    Parameters
    ---------------------
    x, y : coordinates
    a : Array of values on the four nearest points.
        x1, y1 = floor(x), floor(y)
        x2, y2 = ceil(x), ceil(y)
        a[1,1] = value on (x1, y1)
        a[1,2] = value on (x1, y2)
        a[2,1] = value on (x2, y1)
        a[2,2] = value on (x2, y2)

    ---------------------
    Returns
    ----------------
    v: value on (x,y)
    """

    v = T(0.0)
    xs = zeros(2)
    ys = zeros(2)

    xs[1] = x-floor(x)
    xs[2] = 1-xs[1]
    ys[1] = y-floor(y)
    ys[2] = 1-ys[1]

    v = 0.0

    for i=1:2, j=1:2
        f = xs[i] * ys[j]
        v += f * a[3-i, 3-j]
    end

    return v
end

function trilinear_interpolation(x::Real, y::Real, z::Real, a::Array{T, 3}) where T<:AbstractFloat
    """
    Return bilinear interpolation value on point (x,y,z).

    ---------------------
    Parameters
    ---------------------
    x, y, z : coordinates
    a : Array of values on the eight nearest points.
        x1, y1, z1 = floor(x), floor(y), floor(z)
        a[1,1,1] = value on (x1, y1, z1)
        a[2,1,1] = value on (x1+1, y1, z1)
        a[1,2,1] = value on (x1, y1+1, z1)
        a[1,1,2] = value on (x1, y1, z1+1)
        and so on
    ---------------------
    Returns
    ----------------
    v: value on (x,y)
    """
    v = T(0.0)

    xs = zeros(2)
    ys = zeros(2)
    zs = zeros(2)

    xs[1] = x-floor(x)
    xs[2] = 1-xs[1]
    ys[1] = y-floor(y)
    ys[2] = 1-ys[1]
    zs[1] = z-floor(z)
    zs[2] = 1-zs[1]

    for i =1:2, j=1:2, k=1:2
        f = xs[i] * ys[j] * zs[k]
        v += f * a[3-i, 3-j, 3-k]
    end

    return v
end

function around_points(x::Real, y::Real, z::Real, A::Array{T, 3}) where {T<:AbstractFloat}
    """
    Return 2x2x2 values on the eight nearest points of (x,y,z)
    """
    points = zeros(T,2,2,2)
    x1, y1, z1 = floor(Int, x), floor(Int, y), floor(Int, z)
    for li = 1:2, lj=1:2, lk=1:2
        points[li,lj,lk] = A[x1+li-1, y1+lj-1, z1+lk-1]
    end
    return points
end

function around_points(x::Real, y::Real, A::Array{T, 2}) where {T<:AbstractFloat}
    """
    Return 2x2 values on the four nearest points of (x,y)
    """
    points = zeros(T,2,2)
    x1, y1= floor(Int, x), floor(Int, y)
    for li = 1:2, lj=1:2
        points[li,lj] = A[x1+li-1, y1+lj-1]
    end
    return points
end