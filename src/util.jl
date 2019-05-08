@inline function cross_x(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return -x3*y2 + x2*y3
end

@inline function cross_y(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return x3*y1 - x1*y3
end

@inline function cross_z(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return -x2*y1 + x1*y2
end

@inline function cross_product(x1::T, x2::T, x3::T, y1::T, y2::T, y3::T) where {T<:AbstractFloat}
    return (-x3*y2 + x2*y3, x3*y1 - x1*y3, -x2*y1 + x1*y2)
end

@inline function volume(Sx::T, Sy::T, Sz::T, Six::T, Siy::T, Siz::T, Sjx::T, Sjy::T, Sjz::T) where {T<:AbstractFloat}
    tx = Sx * (-Siz * Sjy + Siy * Sjz);
    ty = Sy * (Siz * Sjx - Six * Sjz);
    tz = Sz * (-Siy * Sjx + Six * Sjy);
    return tx + ty + tz;
end
