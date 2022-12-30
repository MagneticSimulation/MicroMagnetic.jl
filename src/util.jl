@inline function dot_product(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return x1*y1 + x2*y2 + x3*y3
end

@inline function cross_x(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return -x3*y2 + x2*y3
end

@inline function cross_y(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return x3*y1 - x1*y3
end

@inline function cross_z(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return -x2*y1 + x1*y2
end

@inline function cross_product(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return (-x3*y2 + x2*y3, x3*y1 - x1*y3, -x2*y1 + x1*y2)
end

@inline function cross_product(x::Tuple{Number, Number, Number}, y::Tuple{Number, Number, Number})
    return (-x[3]*y[2] + x[2]*y[3], x[3]*y[1] - x[1]*y[3], -x[2]*y[1] + x[1]*y[2])
end

#compute a.(bxc) = b.(cxa) = c.(axb)
@inline function volume(Sx::T, Sy::T, Sz::T, Six::T, Siy::T, Siz::T, Sjx::T, Sjy::T, Sjz::T) where {T<:AbstractFloat}
    tx = Sx * (-Siz * Sjy + Siy * Sjz);
    ty = Sy * (Siz * Sjx - Six * Sjz);
    tz = Sz * (-Siy * Sjx + Six * Sjy);
    return tx + ty + tz;
end

#compute the angle defined in equation (1) in paper [PRB 93, 174403 (2016)] or equation (3) in [New J. Phys. 20 (2018) 103014]
@inline function Berg_Omega(ux::T, uy::T, uz::T, vx::T, vy::T, vz::T, wx::T, wy::T, wz::T) where {T<:AbstractFloat}
    b = volume(ux, uy, uz, vx, vy, vz, wx, wy, wz)
    a = 1.0 + (ux*vx + uy*vy + uz*vz) + (ux*wx + uy*wy + uz*wz) + (vx*wx + vy*wy + vz*wz)
    return 2*atan(b, a)
end

function abs!(a::Array{T,1}, b::Array{T,1})  where {T<:AbstractFloat}
    a .= abs.(b)
    return nothing
end

function abs!(a::Array{T,1})  where {T<:AbstractFloat}
    for i in 1:length(a)
        if a[i] < 0
            @inbounds a[i] = -a[i]
        end
    end
    return nothing
end

#The frequency of discrete FFT
#f = [0, 1, ...,   N/2-1,     -N/2, ..., -1] / (d*N)   if n is even
#f = [0, 1, ..., (N-1)/2, -(N-1)/2, ..., -1] / (d*N)   if n is odd
#where d is the data interval.
function fftfreq(N; d=1)
    f = zeros(Float64, N)
    n = ceil(Int, N/2)
    for i = 1:n
        f[i] = (i-1)/(d*N)
    end
    for i = n+1:N
        f[i] = (i-1-N)/(d*N)
    end
    return f
end

"""
This is a helper function to convert m to cylindrical coordinates (mr, mt, mz) for spins 
uniformly located in CylindricalTubeMeshGPU
"""
function convert_m_to_cylindrical(m::Array{T, 1}, nr::Int64, nz::Int64) where {T<:AbstractFloat}
    N = nr*nz
    mc = zeros(T, 3, N)
    b = reshape(m, 3, N)
    for i = 1:N
        theta = 2*pi*(i-1)/nr
        mc[1, i] = b[1, i]*cos(theta) + b[2, i]*sin(theta)
        mc[2, i] = -b[1, i]*sin(theta) + b[2, i]*cos(theta)
        mc[3, i] = b[3, i]
    end
    return reshape(mc, 3*N)
end

function partial_xy(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    ngbs = mesh.ngbs
    nxyz = mesh.nxyz
    pxm = zeros(T, 3*nxyz)
    pym = zeros(T, 3*nxyz)
    for i = 1:nxyz
      j = 3*i-2
      #x-direction
      i1 = ngbs[1,i]
      i2 = ngbs[2,i]
      factor = i1*i2>0 ? 1/(2*dx) : 1/dx
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      j1 = 3*i1-2
      j2 = 3*i2-2
      pxm[j] = (m[j2] - m[j1]) * factor
      pxm[j+1]  = (m[j2+1] - m[j1+1]) * factor
      pxm[j+2]  = (m[j2+2] - m[j1+2]) * factor


      #y-direction
      i1 = ngbs[3,i]
      i2 = ngbs[4,i]
      factor = i1*i2>0 ? 1/(2*dy) : 1/dy
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      j1 = 3*i1-2
      j2 = 3*i2-2
      pym[j] = (m[j2] - m[j1]) * factor
      pym[j+1]  = (m[j2+1] - m[j1+1]) * factor
      pym[j+2]  = (m[j2+2] - m[j1+2]) * factor
  end
  return  pxm, pym
end


function partial_x(u::Array{T, 1}, mesh::Mesh, Ms::Array{T, 1}) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    ngbs = mesh.ngbs
    nxyz = mesh.nxyz
    px = zeros(T, nxyz)
    for i = 1:nxyz
      #x-direction
      i1 = ngbs[1,i]
      i2 = ngbs[2,i]
      i1>0 && Ms[i1]<0 && (i1 = -1)
      i2>0 && Ms[i2]<0 && (i2 = -1)
      factor = i1*i2>0 ? 1/(2*dx) : 1/dx
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      px[i] = (u[i2] - u[i1]) * factor
    end
  return px
end

function partial_y(u::Array{T, 1}, mesh::Mesh, Ms::Array{T, 1}) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    ngbs = mesh.ngbs
    nxyz = mesh.nxyz
    py = zeros(T, nxyz)
    for i = 1:nxyz
      #y-direction
      i1 = ngbs[3,i]
      i2 = ngbs[4,i]
      i1>0 && Ms[i1]<0 && (i1 = -1)
      i2>0 && Ms[i2]<0 && (i2 = -1)
      factor = i1*i2>0 ? 1/(2*dy) : 1/dy
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      py[i] = (u[i2] - u[i1]) * factor
    end
  return py
end

function partial_z(u::Array{T, 1}, mesh::Mesh, Ms::Array{T, 1}) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    ngbs = mesh.ngbs
    nxyz = mesh.nxyz
    pz = zeros(T, nxyz)
    for i = 1:nxyz
      #y-direction
      i1 = ngbs[5,i]
      i2 = ngbs[6,i]
      i1>0 && Ms[i1]<0 && (i1 = -1)
      i2>0 && Ms[i2]<0 && (i2 = -1)
      factor = i1*i2>0 ? 1/(2*dz) : 1/dz
      i1 < 0 && (i1 = i)
      i2 < 0 && (i2 = i)
      pz[i] = (u[i2] - u[i1]) * factor
    end
  return pz
end
