using PyCall
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

function radon_3d_scalar(s::Array{T,3}, angle::Number, tilt_axis::String) where {T<:AbstractFloat}
    # Note: Array s is expected to be odd-number-sized. For example: (3, 129,129,129) !
    np = pyimport("numpy")
    trans = pyimport("skimage.transform")
    radon = trans.radon
    (nx,ny,nz) = size(s)
    projection = zeros(T,nx,nx)
    angle = angle-90
    angle = np.array([angle])
    if lowercase(tilt_axis) == "x"
        for i = 1:nx
            copyto!(view(projection,i,:), radon(s[i,:,:], angle, circle=true))
        end
    elseif lowercase(tilt_axis) == "y"
        for i = 1:ny
            copyto!(view(projection,:,i), radon(s[:,i,:], angle, circle=true))
        end
    else
        @error("\"tilt_axis\" should be \"alpha\" or \"beta\"")
    end
    return projection
end

function radon_3d_scalar(s::Array{T,3}, angle::Number, tilt_axis::String) where {T<:Complex}
    np = pyimport("numpy")
    trans = pyimport("skimage.transform")
    radon = trans.radon
    (nx,ny,nz) = size(s)
    projection = zeros(T,nx,nx)
    angle = angle-90
    angle = np.array([angle])
    if lowercase(tilt_axis) == "x"
        for i = 1:nx
            real_part = radon(real.(s[i,:,:]), angle, circle=true)
            imag_part = radon(imag.(s[i,:,:]), angle, circle=true)
            copyto!(view(projection,i,:), real_part+im*imag_part)
        end
    elseif lowercase(tilt_axis) == "y"
        for i = 1:ny
            real_part = radon(real.(s[:,i,:]), angle, circle=true)
            imag_part = radon(imag.(s[:,i,:]), angle, circle=true)
            copyto!(view(projection,:,i), real_part+im*imag_part)
        end
    else
        @error("\"tilt_axis\" should be \"x\" or \"y\"")
    end
    return projection
end

"""
    vector_field_projection(v::Array{T,4}, angle::Number, tilt_axis::String) where {T<:Number}

return a projection vector of the input array.

Input size: (3,N,N,N)
output size: (3, N, N)
angle: Number from -90 to 90
tilt_axis: "x" or "y"
"""

function vector_field_projection(v::Array{T,4}, angle::Number, tilt_axis::String) where {T<:Number}
    (three, N, _, _) = size(v)
    vnew = zeros(T,3,N,N)
    xp = JuMag.radon_3d_scalar(v[1,:,:,:], angle, tilt_axis)
    yp = JuMag.radon_3d_scalar(v[2,:,:,:], angle, tilt_axis)
    zp = JuMag.radon_3d_scalar(v[3,:,:,:], angle, tilt_axis)
    COS, SIN = cos(deg2rad(angle)), sin(deg2rad(angle))
    if tilt_axis == "x"
        copyto!(view(vnew,1,:,:), xp)
        copyto!(view(vnew,2,:,:), COS * yp - SIN * zp)
        copyto!(view(vnew,3,:,:), -COS * zp - SIN * yp)
    elseif tilt_axis == "y"
        copyto!(view(vnew,1,:,:), COS * xp - SIN * zp)
        copyto!(view(vnew,2,:,:), yp)
        copyto!(view(vnew,3,:,:), -COS * zp - SIN * xp)
    end
    return vnew
end

function sum_ovf(ovf; axis=ez)
    if axis == ez
        return sum_z(ovf)
    end
    if axis == ey
        return sum_y(ovf)
    end
    if axis == ex
        return sum_x(ovf)
    end

    @error("axis should be one of ex ey ez!")
end

function sum_z(ovf)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    m = ovf.data
    b = reshape(m,(3,nx,ny,nz))
    mx, my, mz = b[1,:,:,:], b[2,:,:,:], b[3,:,:,:]
    mxp, myp, mzp = zeros(nx,ny),zeros(nx,ny),zeros(nx,ny)
    for k = 1:nz
        mxp .+= mx[:,:,k]
        myp .+= my[:,:,k]
        mzp .+= mz[:,:,k]
    end
    return mxp, myp, mzp
end

function sum_y(ovf)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    m = ovf.data
    b = reshape(m,(3,nx,ny,nz))
    mx, my, mz = b[1,:,:,:], b[2,:,:,:], b[3,:,:,:]
    mxp, myp, mzp = zeros(nx,nz),zeros(nx,nz),zeros(nx,nz)
    for k = 1:ny
        mxp .+= mx[:,k,:]
        myp .+= mz[:,k,:]
        mzp .-= my[:,k,:]
    end
    return mxp, myp, mzp
end

function sum_x(ovf)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    m = ovf.data
    b = reshape(m,(3,nx,ny,nz))
    mx, my, mz = b[1,:,:,:], b[2,:,:,:], b[3,:,:,:]
    mxp, myp, mzp = zeros(ny,nz),zeros(ny,nz),zeros(ny,nz)
    for k = 1:nx
        mxp .-= my[k,:,:]
        myp .+= mz[k,:,:]
        mzp .-= mx[k,:,:]
    end
    return mxp, myp, mzp
end