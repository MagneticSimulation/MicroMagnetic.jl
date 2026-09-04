@inline function dot_product(x1::T, x2::T, x3::T, y1::T, y2::T,
                             y3::T) where {T<:AbstractFloat}
    return x1*y1 + x2*y2 + x3*y3
end

"""
    cross_x(x1, x2, x3, y1, y2, y3)
    cross_y(x1, x2, x3, y1, y2, y3)
    cross_z(x1, x2, x3, y1, y2, y3)

Return the three components of  (x1,x2,x3) × (y1,y2,y3).

NOTE ON TYPE SIGNATURE (historical pitfall — DO NOT re-tighten to `where {T<:Number, x1::T,...}`):

  MicroMagnetic's `build_matrix` (see `src/eigen/eigen.jl`) performs *exact symbolic
  linearisation* of the LLG equation by feeding effective-field kernels a spin array
  whose elements are `Epsilon` / `AddExpr` duals, not plain `Float64`.  Many kernels
  multiply *scalar* coefficient tuples (e.g. the DMI direction masks of ±1 or the
  1/dx pre-factors, which are plain `Float64`) with *symbolic* neighbour-spin values
  inside  cross_*(ax,ay,az, m_kx,m_ky,m_kz).  In these calls the six arguments are
  NOT of the same Julia type:  (Float64,Float64,Float64, AddExpr,AddExpr,AddExpr).

  Hence  cross_x / cross_y / cross_z  MUST accept mixed numeric types.  The "all six
  identical" parametric method below is still valuable (it lets Julia fully inline and
  optimise the normal-Float64 production path), but the *first* three methods below
  are the generic, promotion-based fallbacks that work for any Number mixture.
"""
@inline cross_x(x1, x2, x3, y1, y2, y3) = -x3 * y2 + x2 * y3
@inline cross_y(x1, x2, x3, y1, y2, y3) =  x3 * y1 - x1 * y3
@inline cross_z(x1, x2, x3, y1, y2, y3) = -x2 * y1 + x1 * y2

# Homogeneous fast-path (same concrete numeric type across all 6 arguments):
@inline function cross_x(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return -x3*y2 + x2*y3
end

@inline function cross_y(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return x3*y1 - x1*y3
end

@inline function cross_z(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:Number}
    return -x2*y1 + x1*y2
end

@inline function cross_product(x1::T, x2::T, x3::T, y1::T, y2::T,
                               y3::T) where {T<:AbstractFloat}
    return (-x3*y2 + x2*y3, x3*y1 - x1*y3, -x2*y1 + x1*y2)
end

@inline function cross_product(x::Array{T,1}, y::Array{T,1}, n_total::Int) where {T<:AbstractFloat}
    x = reshape(x, 3, n_total)
    y = reshape(y, 3, n_total)
    result = similar(x)
    @inbounds for i in 1:n_total
        result[1, i] = -x[3, i] * y[2, i] + x[2, i] * y[3, i]
        result[2, i] =  x[3, i] * y[1, i] - x[1, i] * y[3, i]
        result[3, i] = -x[2, i] * y[1, i] + x[1, i] * y[2, i]
    end
    result = reshape(result, 3 * n_total)
    return result
end

@inline function cross_product(x::Array{T,1}, y::Array{T,1}) where {T<:AbstractFloat}
    return [-x[3]*y[2] + x[2]*y[3], x[3]*y[1] - x[1]*y[3], -x[2]*y[1] + x[1]*y[2]]
end

@inline function cross_product(x::Tuple{Number,Number,Number},
                               y::Tuple{Number,Number,Number})
    return (-x[3]*y[2] + x[2]*y[3], x[3]*y[1] - x[1]*y[3], -x[2]*y[1] + x[1]*y[2])
end

#compute a.(bxc) = b.(cxa) = c.(axb)
@inline function volume(Sx::T, Sy::T, Sz::T, Six::T, Siy::T, Siz::T, Sjx::T, Sjy::T,
                        Sjz::T) where {T<:AbstractFloat}
    tx = Sx * (-Siz * Sjy + Siy * Sjz);
    ty = Sy * (Siz * Sjx - Six * Sjz);
    tz = Sz * (-Siy * Sjx + Six * Sjy);
    return tx + ty + tz;
end

#compute the angle defined in equation (1) in paper [PRB 93, 174403 (2016)] or equation (3) in [New J. Phys. 20 (2018) 103014]
@inline function Berg_Omega(ux::T, uy::T, uz::T, vx::T, vy::T, vz::T, wx::T, wy::T,
                            wz::T) where {T<:AbstractFloat}
    b = volume(ux, uy, uz, vx, vy, vz, wx, wy, wz)
    a = 1.0 + (ux*vx + uy*vy + uz*vz) + (ux*wx + uy*wy + uz*wz) + (vx*wx + vy*wy + vz*wz)
    return 2*atan(b, a)
end

# compute A = A + B 
function vector_add(A::AbstractArray{T,1}, B::AbstractArray{T,1}) where {T<:AbstractFloat}
    @kernel function vector_add_kernel!(a, b)
        i = @index(Global)
        @inbounds a[i] = a[i] + b[i]
    end
    kernel! = vector_add_kernel!(get_backend(A))
    kernel!(A, B; ndrange=length(A))
    return nothing
end

# compute a = a1 + c2*a2 
function vector_add2(a::A, a1::A, a2::A,
                     c2::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add2_kernal!(a, @Const(a1), @Const(a2), c)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c * a2[i]
    end
    kernel! = vector_add2_kernal!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, c2; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3
function vector_add3(a::A, a1::A, a2::A, a3::A, c2::S,
                     c3::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add3_kernal!(a, @Const(a1), @Const(a2), @Const(a3), c2, c3)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i]
    end
    kernel! = vector_add3_kernal!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, a3, c2, c3; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4
function vector_add4(a::A, a1::A, a2::A, a3::A, a4::A, c2::S, c3::S,
                     c4::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add4_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4),
                                         c2, c3, c4)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i]
    end
    kernel! = vector_add4_kernel!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, a3, a4, c2, c3, c4; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5
function vector_add5(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, c2::S, c3::S, c4::S,
                     c5::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add5_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4),
                                         @Const(a5), c2, c3, c4, c5)
        i = @index(Global)
        @inbounds a[i] = a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i]
    end
    kernel! = vector_add5_kernel!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, c2, c3, c4, c5; ndrange=length(a))
    return nothing
end

# compute a = c1*a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5
function vector_add5b(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, c1::S, c2::S, c3::S, c4::S,
                      c5::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add5b_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4),
                                          @Const(a5), c1, c2, c3, c4, c5)
        i = @index(Global)
        @inbounds a[i] = c1 * a1[i] + c2 * a2[i] + c3 * a3[i] + c4 * a4[i] + c5 * a5[i]
    end
    kernel! = vector_add5b_kernel!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, c1, c2, c3, c4, c5; ndrange=length(a))
    return nothing
end

# compute a = a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5 + c6*a6
function vector_add6(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, a6::A, c2::S, c3::S, c4::S,
                     c5::S,
                     c6::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add6_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4),
                                         @Const(a5), @Const(a6), c2, c3, c4, c5, c6)
        i = @index(Global)
        @inbounds a[i] = a1[i] +
                         c2 * a2[i] +
                         c3 * a3[i] +
                         c4 * a4[i] +
                         c5 * a5[i] +
                         c6 * a6[i]
    end
    kernel! = vector_add6_kernel!(get_backend(a), groupsize[])
    kernel!(a, a1, a2, a3, a4, a5, a6, c2, c3, c4, c5, c6; ndrange=length(a))
    return nothing
end

# compute a = c1*a1 + c2*a2 + c3*a3 + c4*a4 + c5*a5 + c6*a6
function vector_add6b(a::A, a1::A, a2::A, a3::A, a4::A, a5::A, a6::A, c1::S, c2::S, c3::S,
                      c4::S, c5::S,
                      c6::S) where {T<:AbstractFloat,S<:AbstractFloat,A<:AbstractArray{T,1}}
    @kernel function vector_add6b_kernel!(a, @Const(a1), @Const(a2), @Const(a3), @Const(a4),
                                          @Const(a5), @Const(a6), c1, c2, c3, c4, c5, c6)
        i = @index(Global)
        @inbounds a[i] = c1 * a1[i] +
                         c2 * a2[i] +
                         c3 * a3[i] +
                         c4 * a4[i] +
                         c5 * a5[i] +
                         c6 * a6[i]
    end
    kernel! = vector_add6b_kernel!(get_backend(a), groupsize[])
    return kernel!(a, a1, a2, a3, a4, a5, a6, c1, c2, c3, c4, c5, c6; ndrange=length(a))
end

#The frequency of discrete FFT
#f = [0, 1, ...,   N/2-1,     -N/2, ..., -1] / (d*N)   if n is even
#f = [0, 1, ..., (N-1)/2, -(N-1)/2, ..., -1] / (d*N)   if n is odd
#where d is the data interval.
function fftfreq(N; d=1)
    f = zeros(Float64, N)
    n = ceil(Int, N/2)
    for i in 1:n
        f[i] = (i-1)/(d*N)
    end
    for i in (n + 1):N
        f[i] = (i-1-N)/(d*N)
    end
    return f
end

function average(x::AbstractArray{T,1}) where {T<:AbstractFloat}
    return sum(x)/length(x)
end

function normalise(a::AbstractArray{T,1}, N::Int64) where {T<:AbstractFloat}
    @kernel function local_kernel!(a)
        id = @index(Global)
        j = 3 * id - 2

        @inbounds m2 = a[j] * a[j] + a[j + 1] * a[j + 1] + a[j + 2] * a[j + 2]
        if m2 > 0
            length::T = 1 / sqrt(m2)
            @inbounds a[j] *= length
            @inbounds a[j + 1] *= length
            @inbounds a[j + 2] *= length
        end
    end

    local_kernel!(get_backend(a), groupsize[])(a; ndrange=N)
    return nothing
end

function compute_dm!(dm::AbstractArray{T,1}, m1::AbstractArray{T,1}, m2::AbstractArray{T,1},
                     N::Int64) where {T<:AbstractFloat}
    @kernel function local_kernel!(c, a, b)
        I = @index(Global)
        j = 3 * I - 2
        @inbounds mx = a[j] - b[j]
        @inbounds my = a[j + 1] - b[j + 1]
        @inbounds mz = a[j + 2] - b[j + 2]
        @inbounds c[I] = sqrt(mx * mx + my * my + mz * mz)
    end

    local_kernel!(get_backend(dm), groupsize[])(dm, m1, m2; ndrange=N)
    return nothing
end

function omega_to_spin(omega::AbstractArray{T,1}, spin::AbstractArray{T,1},
                       spin_next::AbstractArray{T,1}, N::Int64) where {T<:AbstractFloat}
    #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
    #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
    @kernel function local_kernel!(a, b, c)
        I = @index(Global)
        j = 3 * I - 2
        @inbounds w1 = a[j] * 0.5
        @inbounds w2 = a[j + 1] * 0.5
        @inbounds w3 = a[j + 2] * 0.5
        @inbounds m1 = b[j]
        @inbounds m2 = b[j + 1]
        @inbounds m3 = b[j + 2]
        r = 1 + w1 * w1 + w2 * w2 + w3 * w3
        a11 = 1 + w1 * w1 - w2 * w2 - w3 * w3
        a12 = 2 * (w1 * w2 - w3)
        a13 = 2 * (w2 + w1 * w3)
        a21 = 2 * (w1 * w2 + w3)
        a22 = 1 - w1 * w1 + w2 * w2 - w3 * w3
        a23 = -2 * (w1 - w2 * w3)
        a31 = 2 * (-w2 + w1 * w3)
        a32 = 2 * (w1 + w2 * w3)
        a33 = 1 - w1 * w1 - w2 * w2 + w3 * w3
        @inbounds c[j] = (a11 * m1 + a12 * m2 + a13 * m3) / r
        @inbounds c[j + 1] = (a21 * m1 + a22 * m2 + a23 * m3) / r
        @inbounds c[j + 2] = (a31 * m1 + a32 * m2 + a33 * m3) / r
    end

    local_kernel!(get_backend(omega), groupsize[])(omega, spin, spin_next; ndrange=N)
    return nothing
end

function omega_to_spin(omega::Array{Float64,1}, spin::Array{Float64,1},
                       spin_next::Array{Float64,1}, N::Int64)
    #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
    #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
    for i in 0:(N - 1)
        j = 3 * i + 1
        w1 = omega[j] * 0.5
        w2 = omega[j + 1] * 0.5
        w3 = omega[j + 2] * 0.5
        m1 = spin[j]
        m2 = spin[j + 1]
        m3 = spin[j + 2]
        r = 1 + w1 * w1 + w2 * w2 + w3 * w3
        a11 = 1 + w1 * w1 - w2 * w2 - w3 * w3
        a12 = 2 * (w1 * w2 - w3)
        a13 = 2 * (w2 + w1 * w3)
        a21 = 2 * (w1 * w2 + w3)
        a22 = 1 - w1 * w1 + w2 * w2 - w3 * w3
        a23 = -2 * (w1 - w2 * w3)
        a31 = 2 * (-w2 + w1 * w3)
        a32 = 2 * (w1 + w2 * w3)
        a33 = 1 - w1 * w1 - w2 * w2 + w3 * w3
        spin_next[j] = (a11 * m1 + a12 * m2 + a13 * m3) / r
        spin_next[j + 1] = (a21 * m1 + a22 * m2 + a23 * m3) / r
        spin_next[j + 2] = (a31 * m1 + a32 * m2 + a33 * m3) / r
    end
end

function compute_dm!(dm::Array{Float64,1}, m1::Array{Float64,1}, m2::Array{Float64,1},
                     N::Int64)
    for i in 1:N
        j = 3 * i - 2
        mx = m1[j] - m2[j]
        my = m1[j + 1] - m2[j + 1]
        mz = m1[j + 2] - m2[j + 2]
        dm[i] = sqrt(mx * mx + my * my + mz * mz)
    end
    return nothing
end

"""
This is a helper function to convert m to cylindrical coordinates (mr, mt, mz) for spins 
uniformly located in CylindricalTubeMeshGPU
"""
function convert_m_to_cylindrical(m::Array{T,1}, nr::Int64,
                                  nz::Int64) where {T<:AbstractFloat}
    N = nr*nz
    mc = zeros(T, 3, N)
    b = reshape(m, 3, N)
    for i in 1:N
        theta = 2*pi*(i-1)/nr
        mc[1, i] = b[1, i]*cos(theta) + b[2, i]*sin(theta)
        mc[2, i] = -b[1, i]*sin(theta) + b[2, i]*cos(theta)
        mc[3, i] = b[3, i]
    end
    return reshape(mc, 3*N)
end

function partial_xy(m::Array{T,1}, mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    ngbs = mesh.ngbs
    n_total = mesh.n_total
    pxm = zeros(T, 3*n_total)
    pym = zeros(T, 3*n_total)
    for i in 1:n_total
        j = 3*i-2
        #x-direction
        i1 = ngbs[1, i]
        i2 = ngbs[2, i]
        factor = i1*i2>0 ? 1/(2*dx) : 1/dx
        i1 < 0 && (i1 = i)
        i2 < 0 && (i2 = i)
        j1 = 3*i1-2
        j2 = 3*i2-2
        pxm[j] = (m[j2] - m[j1]) * factor
        pxm[j+1] = (m[j2+1] - m[j1+1]) * factor
        pxm[j+2] = (m[j2+2] - m[j1+2]) * factor

        #y-direction
        i1 = ngbs[3, i]
        i2 = ngbs[4, i]
        factor = i1*i2>0 ? 1/(2*dy) : 1/dy
        i1 < 0 && (i1 = i)
        i2 < 0 && (i2 = i)
        j1 = 3*i1-2
        j2 = 3*i2-2
        pym[j] = (m[j2] - m[j1]) * factor
        pym[j+1] = (m[j2+1] - m[j1+1]) * factor
        pym[j+2] = (m[j2+2] - m[j1+2]) * factor
    end
    return pxm, pym
end
