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

function bilinear_interpolation(x::Real, y::Real, a::Array{Complex{T}, 2}) where T<:AbstractFloat
    return bilinear_interpolation(x,y,real.(a)) + 1im*bilinear_interpolation(x,y,imag.(a))
end

function trilinear_interpolation(x::Real, y::Real, z::Real, a::Array{Complex{T}, 3}) where T<:AbstractFloat
    return trilinear_interpolation(x,y,z, real.(a)) + 1im*trilinear_interpolation(x,y,z, imag.(a))
end

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

function warp(img::Array{T,2}, rad::Float64) where {T<:AbstractFloat}
    nx, ny = size(img)
    N = max(nx, ny)
    img_new = zeros(T, N, N)
    COS, SIN = cos(-1*rad), sin(-1*rad)

    xc, yc = (N+1)/2, (N+1)/2
    for i = 1:nx, j=1:ny
        # new coordinates
        xn, yn = i-xc, j-yc

        # old coordinates
        xo = COS * xn - SIN * yn + xc
        yo = SIN * xn + COS * yn + yc

        if inbounds(xo, yo, nx, ny)
            values = zeros(2,2)
            for li = 1:2, lj=1:2
                x1, y1 = floor(Int, xo), floor(Int, yo)
                values[li,lj] = img[x1+li-1, y1+lj-1]
            end
            img_new[i, j] += bilinear_interpolation(xo, yo, values)
        end
    end
    return img_new
end