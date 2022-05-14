"""
    radon(img::Array{T, 2}, theta::Real) where {T<:AbstractFloat}

Radon transform.

0 degree projection direction: x+
positive rotation direction: image towards clockwise

Parameters
------------------------------
img: 2D array sized (N, N). 
theta: Polar angle in rad.

Outputs
----------------------------
projection: 1D array sized (N)
"""
function radon(img::Array{T, 2}, theta::Real) where {T<:AbstractFloat}
    new_img = warp(img, -1*theta)
    projection = sum(new_img[:, :], dims=1)[1, :]

    return projection
end


"""
    radon_x(s::Array{T, 3}, gamma::Real) where T<:Number

3D projection of a SCALAR field rotated around X-axis.

0 degree projection direction: from z+ to z-

Parameters
------------------------------
s: 3D array sized (N, N, N). 
gamma: Euler angle around X-axis.

Outputs
----------------------------
projection: 2D array sized (N,N)
"""
function radon_x(s::Array{T, 3}, gamma::Real) where T<:Number
    new_obj = tilt_x(s, gamma)
    projection = sum(new_obj, dims=3)[:,:,1]
    return projection
end


"""
    radon_x(v::Array{T, 4}, gamma::Real) where T<:Number

3D projection of a VECTOR field rotated around X-axis.

0 degree projection direction: from z+ to z-

Parameters
------------------------------
v: 4D array sized (3,N, N, N). 
gamma: Euler angle around X-axis.

Outputs
----------------------------
projection: 2D array sized (N,N)
"""
function radon_x(v::Array{T,4}, gamma::Real) where {T<:Number}
    vy= radon_x(v[2,:,:,:], gamma)
    vz= radon_x(v[3,:,:,:], gamma)
    return -1 .* (cos(gamma) .* vz + sin(gamma) .* vy)
end



"""
    radon_y(s::Array{T, 3}, beta::Real) where T<:Number

3D projection of a SCALAR field rotated around Y-axis.

0 degree projection direction: from z+ to z-

Parameters
------------------------------
s: 3D array sized (N, N, N). 
beta: Euler angle around Y-axis.

Outputs
----------------------------
projection: 2D array sized (N,N)
"""
function radon_y(s::Array{T, 3}, beta::Real) where T<:Number
    new_obj = tilt_y(s, beta)
    projection = sum(new_obj, dims=3)[:,:,1]
    return projection
end


"""
    radon_y(v::Array{T, 4}, beta::Real) where T<:Number

3D projection of a VECTOR field rotated around X-axis.

0 degree projection direction: from z+ to z-

Parameters
------------------------------
v: 4D array sized (3,N, N, N). 
beta: Euler angle around Y-axis.

Outputs
----------------------------
projection: 2D array sized (N,N)
"""
function radon_y(v::Array{T,4}, beta::Real) where {T<:Number}
    vx= radon_y(v[1,:,:,:], beta)
    vz= radon_y(v[3,:,:,:], beta)
    return -1 .* (cos(beta) .* vz - sin(beta) .* vx)
end


function radon3d(v::Array{T}, theta::Real, axis::String) where {T<:Number}
    if axis == "X"
        return radon_x(v, theta)
    elseif axis == "Y"
        return radon_y(v, theta)
    else
        error("axis must be \"X\" or \"Y\" !")
    end
end

function radon3d_xyz(v::Array{T,4}, theta::Real, axis::String) where {T<:Number}
    vx = radon3d(v[1,:,:,:], theta, axis)
    vy = radon3d(v[2,:,:,:], theta, axis)
    vz = radon3d(v[3,:,:,:], theta, axis)
    st, ct = sin(theta), cos(theta)
    if axis == "X"
        vu = vx
        vv = ct .* vy - st .* vz
        vw = -1 .* (ct .* vz + st .* vy)
        return vu, vv, vw
    elseif axis == "Y"
        vu = ct .* vx + st .* vz
        vv = vy
        vw = -1 .* (ct .* vz - st .* vy)
        return vu, vv, vw
    else
        error("axis must be \"X\" or \"Y\" !")
    end
end

function scalar_projection_Euler(A::Array{T, 3}, 
        alpha::Float64, beta::Float64, gamma::Float64) where T<:Number
    new_obj = tilt(A, alpha, beta, gamma)
    projection = sum(new_obj, dims=3)[:,:,1]
    return projection
end


function vector_projection_Euler(v::Array{T,4}, alpha::Float64=0.0, beta::Float64=0.0, gamma::Float64=0.0) where {T<:Number}
    projection_x= radon_3d_object(v[1,:,:,:], alpha, beta, gamma)
    projection_y= radon_3d_object(v[2,:,:,:], alpha, beta, gamma)
    projection_z= radon_3d_object(v[3,:,:,:], alpha, beta, gamma)

    M = Euler(alpha, beta, gamma)
    phi = M[3,1] .* projection_x + M[3,2] .* projection_y + M[3,3] .* projection_z
    return phi
end

