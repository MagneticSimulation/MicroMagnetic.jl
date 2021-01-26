using PyCall

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
tilt_axis: Should be "alpha" or "beta"
"""
function radon_3d_scalar(s, angle::Number, tilt_axis)
    np = pyimport("numpy")
    trans = pyimport("skimage.transform")
    radon = trans.radon
    (nx,ny,nz) = size(s)
    projection = zeros(nx,nx)
    angle = angle-90
    angle = np.array([angle])
    if lowercase(tilt_axis) == "alpha"
        for i = 1:nx
            projection[i,:] = radon(s[i,:,:], angle, circle=true) 
        end
    elseif lowercase(tilt_axis) == "beta"
        for i = 1:ny
            projection[:,i] = radon(s[:,i,:], angle, circle=true) 
        end
    else
        @error("\"tilt_axis\" should be \"alpha\" or \"beta\"")
    end
    return projection
end

function radon_transform_ovf(ovf::OVF2, angle::Number, tilt_axis; N=128)
    np = pyimport("numpy")
    if angle == 0
        smx,smy,smz = sum_z(ovf)
        new_smx = pad_2d_vec(smx, N)
        new_smy = pad_2d_vec(smy, N)
        new_smz = pad_2d_vec(smz, N)
        return new_smx,new_smy,new_smz
    end
    m = ovf.data
    nx,ny,nz = ovf.xnodes,ovf.ynodes,ovf.znodes
    b = reshape(m,(3,nx,ny,nz))    
    mag_pad = pad_3d_vec(b, N)
    mx,my,mz = mag_pad[1,:,:,:],mag_pad[2,:,:,:],mag_pad[3,:,:,:]
    mxp = radon_3d_scalar(mx, angle, tilt_axis)
    myp = radon_3d_scalar(my, angle, tilt_axis)
    mzp = radon_3d_scalar(mz, angle, tilt_axis)
    COS, SIN = cos(deg2rad(angle)), sin(deg2rad(angle))
    local_mx = mxp
    local_my = COS * myp - SIN * mzp
    local_mz = COS * mzp + SIN * myp

    return local_mx,local_my,local_mz
end

function pad_3d_vec(v, N)
    np = pyimport("numpy")
    l = length(size(v))
    if l == 4
        (dims,nx,ny,nz) = size(v)
        n1,n2,n3 = floor(int,(N-nx)/2), floor(int,(N-ny)/2), floor(int,(N-nz)/2)
        v_new = zeros(dims,N,N,N)
        for dim =1:dims
            v_new[dim,:,:,:] .= np.pad(v[dim,:,:,:], ((n1,n1),(n2,n2),(n3,n3)), "constant")
        end
    elseif l == 3
        (nx,ny,nz) = size(v)
        n1,n2,n3 = floor(int,(N-nx)/2), floor(int,(N-ny)/2), floor(int,(N-nz)/2)
        v_new= np.pad(v, ((n1,n1),(n2,n2),(n3,n3)), "constant")
    end
    return v_new
end

function pad_2d_vec(v, N)
    np = pyimport("numpy")
    l = length(size(v))
    if l == 3
        (dims,nx,ny) = size(v)
        n1,n2 = floor(int,(N-nx)/2), floor(int,(N-ny)/2)
        v_new = zeros(dims,N,N)
        for dim =1:dims
            v_new[dim,:,:] .= np.pad(v[dim,:,:], ((n1,n1),(n2,n2)), "constant")
        end
    elseif l == 2
        (nx,ny) = size(v)
        n1,n2 = floor(Int,(N-nx)/2), floor(Int,(N-ny)/2)
        v_new= np.pad(v, ((n1,n1),(n2,n2)), "constant")
    end
    return v_new
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