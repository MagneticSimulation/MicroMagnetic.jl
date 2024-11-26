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

@inline function cross_product(x::Array{T, 1}, y::Array{T, 1}) where {T<:AbstractFloat}
    return [-x[3]*y[2] + x[2]*y[3], x[3]*y[1] - x[1]*y[3], -x[2]*y[1] + x[1]*y[2]]
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

# compute A = A + B 
function vector_add(A::AbstractArray{T, 1}, B::AbstractArray{T, 1}) where {T<:AbstractFloat}
    @kernel function vector_add_kernel!(a, b)
        i = @index(Global)
        @inbounds a[i] = a[i] + b[i]
    end
    kernel! = vector_add_kernel!(default_backend[])
    kernel!(A, B; ndrange=length(A))
    return nothing
end

# compute a = a1 + c2*a2 
function vector_add2(a::A, a1::A, a2::A, c2::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add2_kernal!(a, @Const(a1), @Const(a2), c)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c * a2[i]
    end
    kernel! = vector_add2_kernal!(default_backend[], groupsize[])
    kernel!(a, a1, a2, c2; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3
function vector_add3(a::A, a1::A, a2::A, a3::A, c2::S, c3::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add3_kernal!(a, @Const(a1), @Const(a2), @Const(a3), c2, c3)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i]
    end
    kernel! = vector_add3_kernal!(default_backend[], groupsize[])
    kernel!(a, a1, a2, a3, c2, c3; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4
function vector_add4(a::A, a1::A, a2::A, a3::A, a4::A, c2::S, c3::S, c4::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add4_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4), c2, c3, c4)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i]
    end
    kernel! = vector_add4_kernel!(default_backend[], groupsize[])
    kernel!(a, a1, a2, a3, a4, c2, c3, c4; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5
function vector_add5(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, c2::S, c3::S, c4::S, c5::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add5_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4), @Const(a5), c2, c3, c4, c5)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i]
    end
    kernel! = vector_add5_kernel!(default_backend[], groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, c2, c3, c4, c5; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5 + c6*a6
function vector_add6(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, a6::A, c2::S, c3::S, c4::S, c5::S, c6::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add6_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4), @Const(a5), @Const(a6), c2, c3, c4, c5, c6)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i] + c6 * a6[i]
    end
    kernel! = vector_add6_kernel!(default_backend[], groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, a6, c2, c3, c4, c5, c6; ndrange=length(a))
    return nothing
end

# compute a = c1*a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5 + c6*a6
function vector_add6b(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, a6::A, c1::S, c2::S, c3::S, c4::S, c5::S, c6::S) where {T<:AbstractFloat, S<:AbstractFloat, A<:AbstractArray{T,1}}
    @kernel function vector_add6b_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4), @Const(a5), @Const(a6), c1, c2, c3, c4, c5, c6)
        i = @index(Global)
        @inbounds a[i] = c1 * a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i] + c6 * a6[i]
    end
    kernel! = vector_add6b_kernel!(default_backend[], groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, a6, c1, c2, c3, c4, c5, c6; ndrange=length(a))
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

function average(x::AbstractArray{T,1}) where {T<:AbstractFloat}
    return sum(x)/length(x)
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
    n_total = mesh.n_total
    pxm = zeros(T, 3*n_total)
    pym = zeros(T, 3*n_total)
    for i = 1:n_total
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
    n_total = mesh.n_total
    px = zeros(T, n_total)
    for i = 1:n_total
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
    n_total = mesh.n_total
    py = zeros(T, n_total)
    for i = 1:n_total
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
    n_total = mesh.n_total
    pz = zeros(T, n_total)
    for i = 1:n_total
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
