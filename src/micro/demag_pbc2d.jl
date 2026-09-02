#True 2D-periodic demag for arbitrary direction pairs (xy, xz, yz), the macro
#pipeline route (the real-space split of Lebecki et al., J. Phys. D 41 (2008)
#175005, cf. demag_pbc1d.jl): the quadrant kernel = the explicit image sum over
#the 2D image lattice (the standard tensors_kernel_*! with two nonzero image
#counts) plus the entry-wise analytic far-field tail -- the continuum integral
#of the dipolar kernel over the images beyond the half-integer critical
#rectangle (Wang et al., Comp. Mater. Sci. 49 (2010) 84, eq. 10-15) -- plus the
#periodic-DC column fix.  The solve reuses the macro padded-FFT pipeline
#unchanged (init_demag); the returned energy is the standard Demag.
#The superseded spectral solver (2D batched FFT + layer convolution, xy only)
#lives at commit 7d5e6295.
#
#Periodic-DC column: the state uniform in the periodic plane carries no bound
#charges, so its plane-summed kernel vanishes for every component except the
#open-axis diagonal, whose self-layer plane-sum equals 1 exactly (the field
#from the uniformly magnetized layer is -M between its faces, 0 outside; the
#package convention H = -sum N M makes the sum +1).  The fix forces these
#analytic targets, mirroring the spectral solver's k=0 override.

#the quarter-plane integrals of the dipolar kernel (Wang 2010, eq. 13);
#module-level functions with explicit arguments (KA-kernel safe, no closures)
@inline function f2d_xx(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    abc / (4pi * Tx * Ty) * p * q / ((p^2 + z^2) * sqrt(p^2 + q^2 + z^2))
end
@inline function f2d_yy(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    abc / (4pi * Tx * Ty) * p * q / ((q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
end
@inline function f2d_zz(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    -abc / (4pi * Tx * Ty) * p * q * (p^2 + q^2 + 2z^2) /
        ((p^2 + z^2) * (q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
end
@inline function f2d_xy(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    -abc / (4pi * Tx * Ty) / sqrt(p^2 + q^2 + z^2)
end
@inline function f2d_xz(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    abc / (4pi * Tx * Ty) * q * z / ((p^2 + z^2) * sqrt(p^2 + q^2 + z^2))
end
@inline function f2d_yz(p::T, q::T, z::T, abc::T, Tx::T, Ty::T) where {T}
    abc / (4pi * Tx * Ty) * p * z / ((q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
end

#the far-field tail of one canonical component: the continuum integral of the
#dipolar kernel over {x'>Xc} ∪ {x'<-Xc} ∪ {y'>Yc} ∪ {y'<-Yc} (the semi-infinite
#strips and the corners, inclusion-exclusion, Wang 2010 eq. 10-12); the
#convention matches the package kernel (H = -sum N M)
function pbc2d_tail(comp::Int, x::T, y::T, z::T, Xc::T, Yc::T, abc::T, Tx::T, Ty::T) where {T}
    p1 = x + Xc; p2 = x - Xc
    q1 = y - Yc; q2 = y + Yc
    if comp == 1
        return f2d_xx(p1, q1, z, abc, Tx, Ty) + f2d_xx(p2, q2, z, abc, Tx, Ty) -
               f2d_xx(p1, q2, z, abc, Tx, Ty) - f2d_xx(p2, q1, z, abc, Tx, Ty)
    elseif comp == 2
        return f2d_yy(p1, q1, z, abc, Tx, Ty) + f2d_yy(p2, q2, z, abc, Tx, Ty) -
               f2d_yy(p1, q2, z, abc, Tx, Ty) - f2d_yy(p2, q1, z, abc, Tx, Ty)
    elseif comp == 3
        return f2d_zz(p1, q1, z, abc, Tx, Ty) + f2d_zz(p2, q2, z, abc, Tx, Ty) -
               f2d_zz(p1, q2, z, abc, Tx, Ty) - f2d_zz(p2, q1, z, abc, Tx, Ty)
    elseif comp == 4
        return f2d_xy(p1, q1, z, abc, Tx, Ty) + f2d_xy(p2, q2, z, abc, Tx, Ty) -
               f2d_xy(p1, q2, z, abc, Tx, Ty) - f2d_xy(p2, q1, z, abc, Tx, Ty)
    elseif comp == 5
        return f2d_xz(p1, q1, z, abc, Tx, Ty) + f2d_xz(p2, q2, z, abc, Tx, Ty) -
               f2d_xz(p1, q2, z, abc, Tx, Ty) - f2d_xz(p2, q1, z, abc, Tx, Ty)
    else
        return f2d_yz(p1, q1, z, abc, Tx, Ty) + f2d_yz(p2, q2, z, abc, Tx, Ty) -
               f2d_yz(p1, q2, z, abc, Tx, Ty) - f2d_yz(p2, q1, z, abc, Tx, Ty)
    end
end

#add the tail of one canonical component to the quadrant tensor built
#`build_idx`-th by init_demag (order xx, yy, zz, xy, xz, yz).
#pair = 1 (xy), 2 (xz), 3 (yz): the canonical frame (X, Y, Z) with Z the open
#axis is (x,y,z) / (x,z,y) / (y,z,x); e.g. for the xz pair the actual xy
#coupling is H_X from M_Z (the canonical XZ), etc.
const PBC2D_PCOMP = ((1, 2, 3, 4, 5, 6), (1, 3, 2, 5, 4, 6), (3, 1, 2, 5, 6, 4))

@kernel function pbc2d_add_tail_kernel!(tensor, pair::Int64, pcomp::Int64,
                                        Xc::Float64, Yc::Float64, abc::Float64,
                                        Tx::Float64, Ty::Float64,
                                        dx::Float64, dy::Float64, dz::Float64)
    i, j, k = @index(Global, NTuple)
    if pair == 1
        x = (i - 1) * dx; y = (j - 1) * dy; z = (k - 1) * dz
    elseif pair == 2
        x = (i - 1) * dx; y = (k - 1) * dz; z = (j - 1) * dy
    else
        x = (j - 1) * dy; y = (k - 1) * dz; z = (i - 1) * dx
    end
    @inbounds tensor[i, j, k] += pbc2d_tail(pcomp, x, y, z, Xc, Yc, abc, Tx, Ty)
end

function pbc2d_add_tail!(tensor, pair::Int, build_idx::Int, Nx::Int, Ny::Int, Nz::Int,
                         dx::Float64, dy::Float64, dz::Float64, nx::Int, ny::Int, nz::Int)
    if pair == 1
        n1, n2 = Nx, Ny; l1, l2 = nx * dx, ny * dy
    elseif pair == 2
        n1, n2 = Nx, Nz; l1, l2 = nx * dx, nz * dz
    else
        n1, n2 = Ny, Nz; l1, l2 = ny * dy, nz * dz
    end
    Xc = (n1 + 0.5) * l1
    Yc = (n2 + 0.5) * l2
    abc = dx * dy * dz
    pcomp = PBC2D_PCOMP[pair][build_idx]
    kernel! = pbc2d_add_tail_kernel!(default_backend[], groupsize[])
    kernel!(tensor, pair, pcomp, Xc, Yc, abc, l1, l2, dx, dy, dz; ndrange=(nx, ny, nz))
    return nothing
end

#force the plane-summed kernel (the sum over the two periodic axes at each
#open-axis position) onto the analytic targets: 0 everywhere except +1 at the
#self layer of the open-axis component; the correction is uniform in the
#periodic plane, so it touches only the (fx=fy=0) mode
@kernel function pbc2d_dc_kernel!(tensor, pair::Int64, self::Int64,
                                  nx::Int64, ny::Int64, nz::Int64)
    t = @index(Global)
    target = (t == 1 && self == 1) ? 1.0 : 0.0
    if pair == 1
        cur = 0.0
        for j in 1:ny, i in 1:nx
            cur += tensor[i, j, t]
        end
        corr = (target - cur) / (nx * ny)
        for j in 1:ny, i in 1:nx
            tensor[i, j, t] += corr
        end
    elseif pair == 2
        cur = 0.0
        for k in 1:nz, i in 1:nx
            cur += tensor[i, t, k]
        end
        corr = (target - cur) / (nx * nz)
        for k in 1:nz, i in 1:nx
            tensor[i, t, k] += corr
        end
    else
        cur = 0.0
        for k in 1:nz, j in 1:ny
            cur += tensor[t, j, k]
        end
        corr = (target - cur) / (ny * nz)
        for k in 1:nz, j in 1:ny
            tensor[t, j, k] += corr
        end
    end
end

function pbc2d_dc_fix!(tensor, pair::Int, build_idx::Int, nx::Int, ny::Int, nz::Int)
    self = (pair == 1 ? (build_idx == 3) : pair == 2 ? (build_idx == 2) : (build_idx == 1)) ? 1 : 0
    ndrange = pair == 1 ? nz : pair == 2 ? ny : nx
    kernel! = pbc2d_dc_kernel!(default_backend[], groupsize[])
    kernel!(tensor, pair, self, nx, ny, nz; ndrange=ndrange)
    return nothing
end
