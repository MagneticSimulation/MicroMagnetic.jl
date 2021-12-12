function warp(img::Array{T,2}, rad::Float64) where {T<:AbstractFloat}
    """
    warp image with counter-clockwise
    ij = R * xy
    so xy = R^{-1} ij
    """
    nx, ny = size(img)
    N = max(nx, ny)
    new_img = zeros(T, N, N)
    COS, SIN = cos(rad), sin(rad)

    xc, yc = (N+1)/2, (N+1)/2
    for i = 1:nx, j=1:ny
        # new coordinates
        xn, yn = i-xc, j-yc

        # old coordinates
        xo = COS * xn + SIN * yn + xc
        yo = -1 * SIN * xn + COS * yn + yc

        if inbounds(xo, yo, nx, ny)
            points = around_points(xo, yo, img)
            new_img[i, j] += bilinear_interpolation(xo, yo, points)
        end
    end
    return new_img
end

function RZ(alpha)
    M = zeros(3,3)
    ca, sa = cos(alpha), sin(alpha)
    M[1,1], M[1,2]  = ca, -1*sa
    M[2,1], M[2,2]  = sa, ca
    M[3,3] = 1
    return M
end

function RY(beta)
    M = zeros(3,3)
    cb, sb = cos(beta), sin(beta)
    M[1,1], M[1,3]  = cb, sb
    M[2,2]  = 1
    M[3,1], M[3,3] = -1*sb, cb
    return M
end

function RX(gamma)
    M = zeros(3,3)
    cg, sg = cos(gamma), sin(gamma)
    M[1,1], = 1
    M[2,2], M[2,3] = cg, -1*sg
    M[3,2], M[3,3] = sg, cg
    return M
end

function Euler(alpha::Float64, beta::Float64, gamma::Float64)
    return(RZ(alpha) * RY(beta) * RX(gamma))
end

function tilt(obj::Array{T,3}, alpha::Float64, beta::Float64, gamma::Float64) where {T<:AbstractFloat}
    """
    Return tilted 3-d object.
    --------------------
    Parameters
    --------------------
    obj: 3D array. Scalar field to be tilt.
    alpha, beta, gamma: Euler angles around z,y,x axis respectively.
    """
    nx, ny, nz = size(obj)
    N = max(nx, ny, nz)
    new_obj = zeros(T, N, N, N)

    xc, yc, zc = (N+1)/2, (N+1)/2, (N+1)/2
    M = Euler(alpha, beta, gamma)
    MT = transpose(M)
    for i = 1:nx, j=1:ny, k=1:nz
        # new coordinates
        xn, yn, zn = i-xc, j-yc, k-zc

        # old coordinates
        xo, yo, zo = (MT * [xn, yn, zn])

        io, jo, ko = xo + xc, yo+zc, zo+zc

        if inbounds(io, jo, ko, nx, ny, nz)
            points = around_points(io, jo, ko, obj)
            new_obj[i, j, k] += trilinear_interpolation(io, jo, ko, points)
        end
    end
    return new_obj
end

function tilt(obj::Array{Complex{T},3}, alpha::Float64, beta::Float64, gamma::Float64) where {T<:AbstractFloat}
    return tilt(real.(obj), alpha, beta, gamma) + 1im .* tilt(imag.(obj), alpha, beta, gamma)
end

function warp(img::Array{Complex{T},2}, rad::Float64) where {T<:AbstractFloat}
    return warp(real.(img), rad) + 1im .* warp(imag.(img), rad)
end

function tilt_vecfld(v::Array{T,4}, alpha::Float64, beta::Float64, gamma::Float64) where {T<:Number}
    """
    Tilt a 3-d vector field with Euler angles.
    """
    (dims, nx, ny, nz) = size(v)
    vx, vy, vz = v[1,:,:,:], v[2,:,:,:], v[3,:,:,:]
    vxn = tilt(vx, alpha, beta, gamma)
    vyn = tilt(vy, alpha, beta, gamma)
    vzn = tilt(vz, alpha, beta, gamma)
    
    M = Euler(alpha, beta, gamma)
    vnew_x = M[1,1] .* vxn + M[1,2] .* vyn + M[1,3] .* vzn
    vnew_y = M[2,1] .* vxn + M[2,2] .* vyn + M[2,3] .* vzn
    vnew_z = M[3,1] .* vxn + M[3,2] .* vyn + M[3,3] .* vzn

    vnew = zeros(T, size(v))
    copyto!(view(vnew, 1, :, :, :), vnew_x)
    copyto!(view(vnew, 2, :, :, :), vnew_y)
    copyto!(view(vnew, 3, :, :, :), vnew_z)
    return vnew
end