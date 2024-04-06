"""
    warp(img::Array{T,2}, theta::Real) where {T<:Number}

Warp a 2D scalar field.
positive direction: counter-clockwise

Parameters
------------------------------
img: 2D array sized (N, N). 
theta: Polar angle in rad.

Outputs
----------------------------
new_img: 2D array sized (N, N)

"""
function warp(img::Array{T,2}, theta::Real) where {T<:Number}
    """
    warp image with counter-clockwise
    ij = R * xy
    so xy = R^{-1} ij
    """
    nx, ny = size(img)
    N = max(nx, ny)
    new_img = zeros(T, N, N)
    COS, SIN = cos(theta), sin(theta)

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




"""
    tilt_x(obj::Array{T,3}, gamma::Real) where {T<:Number}

Tilt a 3D scalar field around x-axis, following Euler angle

Parameters
------------------------------
obj: 3D array sized (N, N, N). 
gamma: Euler angles in rad.

Outputs
----------------------------
new_obj: 3D array sized (N,N,N). Scalar field after rotation

"""
function tilt_x(obj::Array{T,3}, gamma::Real) where {T<:Number}
    nx, ny, nz = size(obj)
    N = max(nx, ny, nz)
    new_obj = zeros(T, N, N, N)
    for i = 1:N
        new_obj[i, :, :] = warp(obj[i, :, :], gamma)
    end

    return new_obj
end





"""
    tilt_y(obj::Array{T,3}, beta::Real) where {T<:Number}

Tilt a 3D scalar field around y-axis, following Euler angle

Parameters
------------------------------
obj: 3D array sized (N, N, N). 
beta: Euler angles in rad.

Outputs
----------------------------
new_obj: 3D array sized (N,N,N). Scalar field after rotation

"""
function tilt_y(obj::Array{T,3}, beta::Real) where {T<:Number}
    nx, ny, nz = size(obj)
    N = max(nx, ny, nz)
    new_obj = zeros(T, N, N, N)
    for j = 1:N
        new_obj[:, j, :] = warp(obj[:, j, :], -1*beta)
    end
    return new_obj
end




"""
    tilt_Euler(obj::Array{T,3}, alpha::Real, beta::Real, gamma::Real) where {T<:Number}

Tilt a 3D scalar field by Euler angle

Parameters
------------------------------
obj: 3D array sized (N, N, N). 
alpha,beta,gamma: Euler angles in rad around Z-axis, Y-axis and X-axis respectively.

Outputs
----------------------------
new_obj: 3D array sized (N,N,N). Scalar field after rotation

"""
function tilt_Euler(obj::Array{T,3}, alpha::Real, beta::Real, gamma::Real) where {T<:Number}
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



"""
    tilt_vecfld(v::Array{T,4}, alpha::Real, beta::Real, gamma::Real) where {T<:Number}

Tilt a 3D VECTOR field by Euler angle

Parameters
------------------------------
obj: 4D array sized (3,N, N, N). 
alpha,beta,gamma: Euler angles in rad around Z-axis, Y-axis and X-axis respectively.

Outputs
----------------------------
vnew: 4D array sized (3,N,N,N). Scalar field after rotation

"""
function tilt_vecfld(v::Array{T,4}, alpha::Real, beta::Real, gamma::Real) where {T<:Number}
    vx, vy, vz = v[1,:,:,:], v[2,:,:,:], v[3,:,:,:]
    vxn = tilt(vx, alpha, beta, gamma)
    vyn = tilt(vy, alpha, beta, gamma)
    vzn = tilt(vz, alpha, beta, gamma)
    
    vnew = zeros(T, size(v))
    M = Euler(alpha, beta, gamma)
    vnew[1,:,:,:] = M[1,1] .* vxn + M[1,2] .* vyn + M[1,3] .* vzn
    vnew[2,:,:,:] = M[2,1] .* vxn + M[2,2] .* vyn + M[2,3] .* vzn
    vnew[3,:,:,:] = M[3,1] .* vxn + M[3,2] .* vyn + M[3,3] .* vzn

    return vnew
end