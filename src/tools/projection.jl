function ovf_projection(fname;axis=ez,slice=-1,style = "mz")
    ovf = read_ovf(fname)
    return ovf_projection(ovf,axis=axis,slice=slice)
end

function ovf_projection(ovf::OVF2;axis=ez,slice=-1)
    m = ovf.data
    nx,ny,nz = ovf.xnodes,ovf.ynodes,ovf.znodes
    b = reshape(m,(3,nx,ny,nz))
    mx,my,mz = b[1,:,:,:],b[2,:,:,:],b[3,:,:,:]

    if slice == -1
        mxp,myp,mzp = Normal_Projection(ovf,axis=axis)
    elseif axis==ez && slice in 1:nz
        mxp,myp,mzp = mx[:,:,slice],my[:,:,slice],mz[:,:,slice]
    elseif axis==ex && slice in 1:nx
        mxp,myp,mzp = -my[slice,:,:],mz[slice,:,:],-mx[slice,:,:]
    elseif axis==ey && slice in 1:ny
        mxp,myp,mzp = mx[:,:,slice],mz[:,:,slice],-my[:,:,slice]
    else
        @error("make slice error, slice should be in 1:layer")
    end

    return mxp,myp,mzp
end

"""
        Normal_Projection(ovf; Nx=-1, Ny=-1, Nz=-1, axis = ez)

return the projection magnetization to x/y/z plane

"""
function Normal_Projection(ovf; Nx=-1, Ny=-1, axis = ez)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    m = ovf.data
    b = reshape(m,(3,nx,ny,nz))
    mx, my, mz = b[1,:,:,:], b[2,:,:,:], b[3,:,:,:]

    if axis == ex
        return Normal_Projection_x(mx,my,mz;Nx=Nx,Ny=Ny)
    end

    if axis == ey
        return Normal_Projection_y(mx,my,mz;Nx=Nx,Ny=Ny)
    end

    if axis == ez
        return Normal_Projection_z(mx,my,mz;Nx=Nx,Ny=Ny)
    end

    @error("axis should be one of ex ey ez")
end

function Normal_Projection_z(mx::Array{Float64,2}, my::Array{Float64,2}, mz::Array{Float64,2}; Nx=-1, Ny=-1)
     (nx,ny) = size(mx)
     mx1 = reshape(mx,(nx,ny,1))
     my1 = reshape(my,(nx,ny,1))
     mz1 = reshape(mz,(nx,ny,1))

     return Normal_Projection_z(mx1,my1,mz1,Nx=Nx,Ny=Ny)
end

function Normal_Projection_y(mx::Array{Float64,2}, my::Array{Float64,2}, mz::Array{Float64,2}; Nx=-1, Ny=-1)
     (nx,nz) = size(mx)
     mx1 = reshape(mx,(nx,1,nz))
     my1 = reshape(my,(nx,1,nz))
     mz1 = reshape(mz,(nx,1,nz))

     return Normal_Projection_y(mx1,my1,mz1,Nx=Nx,Ny=Ny)
end

function Normal_Projection_x(mx::Array{Float64,2}, my::Array{Float64,2}, mz::Array{Float64,2}; Nx=-1, Ny=-1)
     (ny,nz) = size(mx)
     mx1 = reshape(mx,(1,ny,nz))
     my1 = reshape(my,(1,ny,nz))
     mz1 = reshape(mz,(1,ny,nz))

     return Normal_Projection_x(mx1,my1,mz1,Nx=Nx,Ny=Ny)
end

function Normal_Projection_z(mx::Array{Float64,3}, my::Array{Float64,3}, mz::Array{Float64,3}; Nx=-1, Ny=-1)
    (nx,ny,nz) = size(mx)
    Nx = Nx > nx ? Nx : nx
    Ny = Ny > ny ? Ny : ny

    smx,smy,smz = zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)
    startx,starty = Int(floor((Nx-nx)/2)),Int(floor((Ny-ny)/2))

    for i = 1:nx, j= 1:ny
        x,y = startx+i,starty+j
        for k = 1:nz
            smx[x,y] += mx[i,j,k]
            smy[x,y] += my[i,j,k]
            smz[x,y] += mz[i,j,k]
        end
    end

    return smx,smy,smz
end

function Normal_Projection_x(mx::Array{Float64,3}, my::Array{Float64,3}, mz::Array{Float64,3}; Nx=-1, Ny=-1)
    (nx,ny,nz) = size(mx)
    Nx = Nx > ny ? Nx : ny
    Ny = Ny > nz ? Ny : nz

    smx,smy,smz = zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)
    startx,starty = Int(floor((Nx-ny)/2)),Int(floor((Ny-nz)/2))

    for j = 1:ny, k= 1:nz
        x,y = startx+ny+1-j,starty+k
        for i = 1:nx
            smx[x,y] += -my[i,j,k]
            smy[x,y] += mz[i,j,k]
            smz[x,y] += -mx[i,j,k]
        end
    end

    return smx,smy,smz
end

function Normal_Projection_y(mx::Array{Float64,3}, my::Array{Float64,3}, mz::Array{Float64,3}; Nx=-1, Ny=-1)
    (nx,ny,nz) = size(mx)
    Nx = Nx > nx ? Nx : nx
    Ny = Ny > nz ? Ny : nz

    smx,smy,smz = zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)
    startx,starty = Int(floor((Nx-nx)/2)),Int(floor((Ny-nz)/2))

    for i = 1:nx, k= 1:nz
        x,y = startx+i,starty+k
        for j = 1:ny
            smx[x,y] += mx[i,j,k]
            smy[x,y] += mz[i,j,k]
            smz[x,y] += -my[i,j,k]
        end
    end

    return smx,smy,smz
end

"""
    Make_Projection(ovf::OVF2;beta=0,gamma=0,Nx=-1,Ny=-1)
    
nx ny nz is the shape of m, Nx Ny is the projection size

First step: rotate beta(rad) with axis (0,-1,0)

Second step : rotate gamma(rad) with axis (1,0,0)

positive degree for anti-clockwise,negative for clockwise

```julia
Make_Projection(m,nx,ny,nz,Nx=128,Ny=128,beta=pi/3)
```
"""
function Make_Projection(ovf::OVF2;beta=0,gamma=0,Nx=-1,Ny=-1)
    m = ovf.data
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes

    return Make_Projection(m,nx,ny,nz,beta=beta,gamma=gamma,Nx=Nx,Ny=Ny)
end

function Make_Projection(m::Array{T,1},nx::Int,ny::Int,nz::Int; Nx::Int=-1,Ny::Int=-1,
                        beta::Number=0,gamma::Number=0) where T<:AbstractFloat

    Nx = Nx > 0 ? Nx : nx
    Ny = Ny > 0 ? Ny : ny

    b = reshape(m,(3,nx,ny,nz))
    mx,my,mz = b[1,:,:,:],b[2,:,:,:],b[3,:,:,:]

    mxp,myp,mzp = Make_Projection(mx,my,mz,Nx=Nx,Ny=Ny,beta=beta,gamma=gamma)
    return mxp,myp,mzp
end

function Make_Projection(mx::Array{Float64,3},my::Array{Float64,3},mz::Array{Float64,3};
                          Nx::Int=-1,Ny::Int=-1,beta::Number=0,gamma::Number=0)

    if abs(beta) < 1e-6 && abs(gamma) < 1e-6
        return Normal_Projection_z(mx, my, mz, Nx=Nx, Ny=Ny)
    end

    Nx = Nx > 0 ? Nx : nx
    Ny = Ny > 0 ? Ny : ny

    (nx,ny,nz) = size(mx)
    mxp,myp,mzp = zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)
    for i=1:nx, j= 1:ny, k=1:nz
        x0,y0,z0= i-(nx+1)/2,j-(ny+1)/2,k-(nz+1)/2
        axis1 = [0,-1.0,0]
        x1,y1,z1 = rotation_operator([x0,y0,z0],axis1,float(beta))

        axis2 = rotation_operator([1.0,0,0],axis1,float(beta))
        x2,y2,z2 = rotation_operator([x1,y1,z1],axis2,float(gamma))
        x,y = x2+(Nx+1)/2,y2+(Ny+1)/2
        linear_separate(mxp,myp,mzp,mx[i,j,k],my[i,j,k],mz[i,j,k],x,y,Nx,Ny,beta,gamma)
    end
    return mxp,myp,mzp
end

function linear_separate(mxp,myp,mzp,mx,my,mz,x,y,Nx,Ny,beta,gamma)
    cb,sb,cg,sg = cos(beta),sin(beta),cos(gamma),sin(gamma)
    if 1 <= x <= Nx && 1 <= y <= Ny
        xf,yf = Int(floor(x)),Int(floor(y))
        xd,yd = x-xf,y-yf
        xp,yp = 1-xd,1-yd
        mx1,my1,mz1 = mx*cb-mz*sb, my, mz*cb+mx*sb
        mxx,myy,mzz = mx1, my1*cg-mz1*sg, my1*sg+mz1*cg

        mxp[xf,yf] += mxx*xp*yp
        myp[xf,yf] += myy*xp*yp
        mzp[xf,yf] += mzz*xp*yp

        if yf+1 <= Ny
            mxp[xf,yf+1] += mxx*xp*yd
            myp[xf,yf+1] += myy*xp*yd
            mzp[xf,yf+1] += mzz*xp*yd
        end

        if xf+1 <= Nx
            mxp[xf+1,yf] += mxx*xd*yp
            myp[xf+1,yf] += myy*xd*yp
            mzp[xf+1,yf] += mzz*xd*yp
        end
        
        if yf+1 <= Ny && xf+1 <= Nx
            mxp[xf+1,yf+1] += mxx*xd*yd
            myp[xf+1,yf+1] += myy*xd*yd
            mzp[xf+1,yf+1] += mzz*xd*yd
        end
    end
end