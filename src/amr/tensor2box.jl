# General two-box Newell demag tensor: a source box (sx,sy,sz) and a target
# box (tx,ty,tz) with center-to-center offset (x,y,z). Reduces to the
# package's 27-point equal-size formula when the sizes coincide (validated),
# and to the dipole far field at large separations.
#
# Corner structure per dimension: target corners ±tx/2 minus source corners
# ±sx/2 → four differences {(tx+sx)/2, (tx−sx)/2, −(tx−sx)/2, −(tx+sx)/2}
# with signs (+,−,−,+). The tensor component is the product corner sum of the
# Newell antiderivative (newell_f for the diagonal, newell_g for the
# off-diagonal components), normalized by the TARGET box volume (the field
# averaged over the target cell, per unit source magnetization -- the same
# convention as the package's equal-size formula, and energy-pairing
# consistent: N_ij*V_i = N_ji*V_j).

@inline function _t2b_diffs(s::Float64, t::Float64)
    return ((t + s) / 2, (t - s) / 2, -(t - s) / 2, -(t + s) / 2), (1.0, -1.0, -1.0, 1.0)
end

"""General two-box demag tensor, component `comp` (1=xx, 2=yy, 3=zz, 4=xy,
5=xz, 6=yz). (x,y,z): center-to-center offset; (sx,sy,sz): source box full
sizes; (tx,ty,tz): target box full sizes. Normalization follows the package
convention ( validated against `demag_tensor_xx` for equal sizes)."""
function demag_tensor_2box(comp::Int, x::Float64, y::Float64, z::Float64,
                           sx::Float64, sy::Float64, sz::Float64,
                           tx::Float64, ty::Float64, tz::Float64)
    (dxs, wx) = _t2b_diffs(sx, tx)
    (dys, wy) = _t2b_diffs(sy, ty)
    (dzs, wz) = _t2b_diffs(sz, tz)
    acc = 0.0
    @inbounds for kz in 1:4, kj in 1:4, ki in 1:4
        acc += wx[ki] * wy[kj] * wz[kz] *
               _t2b_F(comp, x + dxs[ki], y + dys[kj], z + dzs[kz])
    end
    return -acc / (4pi * tx * ty * tz)
end

# Per-component antiderivative and axis permutation, mirroring the package:
# xx = f(x,y,z); yy = f(y,x,z); zz = f(z,y,x); xy = g(x,y,z);
# xz = xy(x,z,y); yz = xy(y,z,x).
@inline function _t2b_F(comp::Int, x::Float64, y::Float64, z::Float64)
    if comp == 1
        return newell_f(x, y, z)
    elseif comp == 2
        return newell_f(y, x, z)
    elseif comp == 3
        return newell_f(z, y, x)
    elseif comp == 4
        return newell_g(x, y, z)
    elseif comp == 5
        return newell_g(x, z, y)
    else
        return newell_g(y, z, x)
    end
end

# ---------------------------------------------------------------- calibration

"""The package's equal-size tensor via my general formula (for validation)."""
function t2b_equal(comp::Int, x, y, z, h)
    return demag_tensor_2box(comp, x, y, z, h, h, h, h, h, h)
end
