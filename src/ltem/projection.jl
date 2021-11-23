using PaddedViews

"""
This part is based on radon from skimage.transform. So you need to install skikit-image.

The operation is :

```julia
using Conda
Conda.install("skikit-image")
```
"""


"""
Return a projection of a 3-d scalar field.
------------------------------------------
s:The scalar field, which should be 3-d padded. The shape should be (N,N,N).
angle: tilt_angle in degrees.
tilt_axis: Should be "x" or "y"
"""

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

function radon(A::Array{T, 2}, thetas) where {T<:AbstractFloat}
    """
    Perform radon transform.

    thetas are polar angles of projection backplane in degree
    """
    mrad = -1 * deg2rad.(thetas)
    nx,ny = size(A)
    N = max(nx,ny)
    L = length(thetas)
    sinogram = zeros(T, N, L)
    for a = 1:L
        new_img = warp(A, mrad[a])
        for i = 1:N
            sinogram[i, a] += sum(new_img[i, :])
        end
    end
    return sinogram
end

function radon_3d(A::Array{T,3}, angles, tilt_axis::String) where {T<:AbstractFloat}
    """
    Global coordinates: x y z . The project ray is always from z+ to z-.
    Local coordinates: u v w
    
    Definition of angles: When the sight is along rotate axis, the angle from z+ direction to w+ direction in counter-clockwise.
    """
    (nx,ny,nz) = size(A)
    (L,) = length(angles)
    sinogram_3d = zeros(T, nx, nx, L)

    if lowercase(tilt_axis) == "x"
        for i = 1:nx
            copyto!(view(sinogram_3d, i, :, :), radon(A[i, :, :], angles ) )
        end
    elseif lowercase(tilt_axis) == "y"
        for j = 1:ny
            copyto!(view(sinogram_3d, :, j, :), radon(A[:, j, :], -1 * angles) )
        end
    else
        @error("\"tilt_axis\" should be \"x\" or \"y\"")
    end
    return sinogram_3d
end

function radon_3d(s::Array{T,3}, angles, tilt_axis::String) where {T<:Complex}
    return radon_3d(real.(s),angles,tilt_axis) + 1im*radon_3d(imag.(s),angles,tilt_axis)
end

"""
    vector_field_projection(v::Array{T,4}, angle::Number, tilt_axis::String) where {T<:Number}

return a projection vector of the input array.

Input size: (3,N,N,N)
output size: (3, N, N)
angle: Number from -90 to 90
tilt_axis: "x" or "y"
"""

function vector_field_projection(v::Array{T,4}, angles, tilt_axis::String) where {T<:Number}
    (three, N, _, _) = size(v)
    (L,) = length(angles)
    phi = zeros(T,N,N,L)
    COS, SIN = cos.(deg2rad.(angles)), sin.(deg2rad.(angles))
    if lowercase(tilt_axis) == "x"
        sinogram_y = radon_3d(v[2,:,:,:], angles, tilt_axis)
        sinogram_z = radon_3d(v[3,:,:,:], angles, tilt_axis)
        for i = 1:L
            copyto!(view(phi,:,:,i), -COS[i] * sinogram_z[:,:,i] + SIN[i] * sinogram_y[:,:,i])
        end
    elseif lowercase(tilt_axis) == "y"
        sinogram_x = radon_3d(v[1,:,:,:], angles, tilt_axis)
        sinogram_z = radon_3d(v[3,:,:,:], angles, tilt_axis)
        for i = 1:L
            copyto!(view(phi,:,:,i), -COS[i] * sinogram_z[:,:,i] - SIN[i] * sinogram_x[:,:,i])
        end
    end
    return phi
end

