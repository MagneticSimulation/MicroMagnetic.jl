@kernel function zeeman_kernel!(@Const(m), @Const(h), energy, @Const(Ms), volume::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds energy[id] = -mu_0 * Ms[id] * volume * (m[j+1] * h[j+1] + m[j+2] * h[j+2] + m[j+3] * h[j+3])
end


@kernel function time_zeeman_kernel!(@Const(m), h, @Const(h_static), energy, @Const(Ms), volume::T, fx::T, fy::T, fz::T) where {T<:AbstractFloat}
    id = @index(Global)
    j = 3 * (id - 1)
    @inbounds h[j+1] = h_static[j+1] * fx
    @inbounds h[j+2] = h_static[j+2] * fy
    @inbounds h[j+3] = h_static[j+3] * fz
    @inbounds energy[id] = -mu_0 * Ms[id] * volume * (m[j+1] * h[j+1] + m[j+2] * h[j+2] + m[j+3] * h[j+3])
end