#True 1D-periodic demag (periodic along one axis, the other two open), following
#Lebecki et al., J. Phys. D 41 (2008) 175005: the kernel entry is the explicit
#image sum |p| <= Ic along the periodic axis plus the analytic far-field tail --
#the continuum integral of the dipolar kernel over the images beyond the
#half-integer boundary Xc = (Ic+0.5)*L (Euler-Maclaurin boundary error ~ Xc^-3,
#measured ~1e-4 in kernel units at Ic=2, see sandbox/pbc1d_tail_check.jl).
#The solve reuses the macro padded-FFT pipeline unchanged (init_demag): the
#periodic axis wraps circularly in the image-summed kernel, the open axes stay
#zero padded.  add_demag(; pbc1d=true, Nx=Ic) selects the axis via which one of
#Nx/Ny/Nz is nonzero.
#
#The tail is expressed in the axis-permuted frame (periodic axis = "x", the
#transverse offsets = (v, w)): with the antiderivative F_ab(u) = int_0^u k_ab of
#the dipolar kernel and its full-line integral Gtot_ab,
#    T_ab(x0, v, w) = (dv*dw/(4pi*n)) * [Gtot_ab - F_ab(x0+Xc) + F_ab(x0-Xc)],
#where x0 is the longitudinal displacement of the kernel entry, n the cell
#count along the axis and (dv, dw) the transverse cell dims (so dv*dw/n =
#V_cell/L, the image density).  The trace of T vanishes.

#F_ab(u): antiderivative of the dipolar kernel k_ab(u, v, w) from 0 (no 1/4pi;
#k_1 = (2u^2-rho^2)/R^5, k_2/3 = the transverse diagonals, k_4/5 = u*v, u*w,
#k_6 = v*w with R^2 = u^2+v^2+w^2); comp is in the permuted frame
function pbc1d_antider(comp::Int, u::Float64, v::Float64, w::Float64)
    rho2 = v^2 + w^2
    R = sqrt(u^2 + rho2)
    if comp == 1
        return -u / R^3
    end
    rho = sqrt(rho2)
    A3 = u / (rho2 * R)
    A5 = u * (2u^2 + 3rho2) / (3 * rho2^2 * R^3)
    if comp == 2
        return 3v^2 * A5 - A3
    elseif comp == 3
        return 3w^2 * A5 - A3
    elseif comp == 4
        return v * (1 / rho^3 - 1 / R^3)
    elseif comp == 5
        return w * (1 / rho^3 - 1 / R^3)
    else
        return 3v * w * A5
    end
end

#full-line integral of the dipolar kernel (the F(inf)-F(-inf) constant)
function pbc1d_gtot(comp::Int, v::Float64, w::Float64)
    rho2 = v^2 + w^2
    if comp == 1 || comp == 4 || comp == 5
        return 0.0
    elseif comp == 2
        return 2 * (v^2 - w^2) / rho2^2
    elseif comp == 3
        return 2 * (w^2 - v^2) / rho2^2
    else
        return 4v * w / rho2^2
    end
end

#the tail value of one permuted-frame component (scalar, used by the kernel);
#the package kernel convention is H = -sum N M, i.e. N = -(the field kernel),
#so the tail enters with the sign opposite to the field-kernel integral
#(the same convention as the N^cont of Lebecki et al., eq. B.6-B.7)
function pbc1d_tail(comp::Int, x0::Float64, v::Float64, w::Float64,
                    dv::Float64, dw::Float64, n::Int64, Xc::Float64)
    pref = dv * dw / (4pi * n)
    rho2 = v^2 + w^2
    if rho2 == 0  #on-axis limit (the antiderivative differences stay finite)
        s = 1 / (x0 + Xc)^2 + 1 / (Xc - x0)^2
        if comp == 1
            return -pref * s
        elseif comp == 2 || comp == 3
            return 0.5 * pref * s
        else
            return 0.0
        end
    end
    return -pref * (pbc1d_gtot(comp, v, w) - pbc1d_antider(comp, x0 + Xc, v, w) +
                    pbc1d_antider(comp, x0 - Xc, v, w))
end

#add the tail to the quadrant tensor built `build_idx`-th by init_demag
#(order xx, yy, zz, xy, xz, yz); axis = the periodic axis, Nimages = Ic.
#The permuted frame (X, Y, Z) = the axis and the two transverse dims:
#axis=1: (x,y,z); axis=2: (y,z,x); axis=3: (z,x,y) -- e.g. for axis=2 the
#actual xy coupling is H_X from M_Z (the permuted XZ), etc.
const PBC1D_PCOMP = ((1, 2, 3, 4, 5, 6), (3, 1, 2, 5, 6, 4), (2, 3, 1, 6, 4, 5))

@kernel function pbc1d_add_tail_kernel!(tensor, axis::Int64, pcomp::Int64, Xc::Float64,
                                        dx::Float64, dy::Float64, dz::Float64,
                                        n::Int64)
    i, j, k = @index(Global, NTuple)
    if axis == 1
        x0 = (i - 1) * dx; v = (j - 1) * dy; w = (k - 1) * dz
        dv = dy; dw = dz
    elseif axis == 2
        x0 = (j - 1) * dy; v = (k - 1) * dz; w = (i - 1) * dx
        dv = dz; dw = dx
    else
        x0 = (k - 1) * dz; v = (i - 1) * dx; w = (j - 1) * dy
        dv = dx; dw = dy
    end
    @inbounds tensor[i, j, k] += pbc1d_tail(pcomp, x0, v, w, dv, dw, n, Xc)
end

function pbc1d_add_tail!(tensor, axis::Int, build_idx::Int, Nimages::Int,
                         dx::Float64, dy::Float64, dz::Float64, nx::Int, ny::Int, nz::Int)
    n = axis == 1 ? nx : axis == 2 ? ny : nz
    daxis = axis == 1 ? dx : axis == 2 ? dy : dz
    Xc = (Nimages + 0.5) * n * daxis
    pcomp = PBC1D_PCOMP[axis][build_idx]
    kernel! = pbc1d_add_tail_kernel!(get_backend(tensor), groupsize[])
    kernel!(tensor, axis, pcomp, Xc, dx, dy, dz, n; ndrange=(nx, ny, nz))
    return nothing
end

#The periodic-DC column (the response to the magnetization uniform along the
#axis) is replaced by the direct axial sum of the raw kernel over |m| <= M cells
#plus the analytic dipolar tail beyond: the di-sum of images+tail carries the
#tail's boundary error, which the uniform mode picks up with full weight
#(a few % of Ms).  The direct sum has unit-cell spacing, so its own midpoint
#error at M ~ 64 cells is negligible; the odd-in-axis components have zero DC
#sum and the fix is a no-op for them.
const PBC1D_DC_M = 64

@kernel function pbc1d_dc_kernel!(tensor, axis::Int64, bcomp::Int64, pcomp::Int64,
                                  M::Int64, Xc::Float64, dx::Float64, dy::Float64,
                                  dz::Float64, nx::Int64, ny::Int64, nz::Int64)
    a2, a3 = @index(Global, NTuple)
    naxis = axis == 1 ? nx : (axis == 2 ? ny : nz)
    cur = zero(tensor[1, 1, 1])
    for t in 1:naxis
        if axis == 1
            cur += tensor[t, a2, a3]
        elseif axis == 2
            cur += tensor[a2, t, a3]
        else
            cur += tensor[a2, a3, t]
        end
    end
    dsum = 0.0
    for m in -M:M
        if axis == 1
            xa = m * dx; ya = (a2 - 1) * dy; za = (a3 - 1) * dz
        elseif axis == 2
            xa = (a2 - 1) * dx; ya = m * dy; za = (a3 - 1) * dz
        else
            xa = (a2 - 1) * dx; ya = (a3 - 1) * dy; za = m * dz
        end
        if bcomp == 1
            dsum += demag_tensor_xx(xa, ya, za, dx, dy, dz)
        elseif bcomp == 2
            dsum += demag_tensor_yy(xa, ya, za, dx, dy, dz)
        elseif bcomp == 3
            dsum += demag_tensor_zz(xa, ya, za, dx, dy, dz)
        elseif bcomp == 4
            dsum += demag_tensor_xy(xa, ya, za, dx, dy, dz)
        elseif bcomp == 5
            dsum += demag_tensor_xz(xa, ya, za, dx, dy, dz)
        else
            dsum += demag_tensor_yz(xa, ya, za, dx, dy, dz)
        end
    end
    if axis == 1
        v = (a2 - 1) * dy; w = (a3 - 1) * dz; dv = dy; dw = dz
    elseif axis == 2
        v = (a3 - 1) * dz; w = (a2 - 1) * dx; dv = dz; dw = dx
    else
        v = (a2 - 1) * dx; w = (a3 - 1) * dy; dv = dx; dw = dy
    end
    dsum += pbc1d_tail(pcomp, 0.0, v, w, dv, dw, 1, Xc)
    corr = (dsum - cur) / naxis
    for t in 1:naxis
        if axis == 1
            tensor[t, a2, a3] += corr
        elseif axis == 2
            tensor[a2, t, a3] += corr
        else
            tensor[a2, a3, t] += corr
        end
    end
end

function pbc1d_dc_fix!(tensor, axis::Int, build_idx::Int,
                       dx::Float64, dy::Float64, dz::Float64, nx::Int, ny::Int, nz::Int)
    daxis = axis == 1 ? dx : axis == 2 ? dy : dz
    Xc = (PBC1D_DC_M + 0.5) * daxis
    pcomp = PBC1D_PCOMP[axis][build_idx]
    kernel! = pbc1d_dc_kernel!(get_backend(tensor), groupsize[])
    ndrange = axis == 1 ? (ny, nz) : axis == 2 ? (nx, nz) : (nx, ny)
    kernel!(tensor, axis, build_idx, pcomp, PBC1D_DC_M, Xc, dx, dy, dz, nx, ny, nz;
            ndrange=ndrange)
    return nothing
end
