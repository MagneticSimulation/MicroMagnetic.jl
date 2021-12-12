using PaddedViews

function vector_padding(v::Array{T,4}, Nx::Int, Ny::Int, Nz::Int) where {T<:Number}
    (dims,nx,ny,nz) = size(v)
    xp = floor(Int,(Nx-nx)/2)
    yp = floor(Int,(Ny-ny)/2)
    zp = floor(Int,(Nz-nz)/2)
    b_new = Array{T}(PaddedView(0,v,(1:dims,1:Nx,1:Ny,1:Nz),(1:dims,xp+1:xp+nx, yp+1:yp+ny, zp+1:zp+nz)))
    return b_new #return 4-d array
end

function vector_crop(v::Array{T,4}, nx::Int,ny::Int,nz::Int)  where {T<:Number}
    (dims,Nx,Ny,Nz) = size(v)
    xp = floor(Int,(Nx-nx)/2)
    yp = floor(Int,(Ny-ny)/2)
    zp = floor(Int,(Nz-nz)/2)
    b_new = v[:, xp+1:xp+nx, yp+1:yp+ny, zp+1:zp+nz]
    return b_new #return 4-d array
end

function radon(A::Array{T, 2}, theta::Float64) where {T<:AbstractFloat}
    """
    Perform radon transform.

    theta is polar angles in rad, with clockwise as positive direction
    """

    new_img = warp(A, theta)
    sinogram = sum(new_img[:, :], dims=2)[:, 1]

    return sinogram
end

function radon_3d_object(A::Array{T, 3}, alpha::Float64, beta::Float64, gamma::Float64) where T<:Number
    """
    Tilt a 3D object with Euler angles, and sum along z-direction
    """

    new_obj = tilt(A, alpha, beta, gamma)
    projection = sum(new_obj, dims=3)[:,:,1]

    return projection
end


function radon_vecfld(v::Array{T,4}, alpha::Float64=0.0, beta::Float64=0.0, gamma::Float64=0.0) where {T<:Number}
    """
    Return the projection component of a tilted vector field.
    """
    projection_x= radon_3d_object(v[1,:,:,:], alpha, beta, gamma)
    projection_y= radon_3d_object(v[2,:,:,:], alpha, beta, gamma)
    projection_z= radon_3d_object(v[3,:,:,:], alpha, beta, gamma)

    M = Euler(alpha, beta, gamma)
    phi = M[3,1] .* projection_x + M[3,2] .* projection_y + M[3,3] .* projection_z
    return phi
end

