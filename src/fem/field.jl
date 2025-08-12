

function effective_field(zee::Zeeman, sim::MicroSimFE, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    v_coeff = sim.mesh.unit_length^3

    back = default_backend[]
    zeeman_kernel!(back, groupsize[])(spin, zee.field, zee.energy, sim.L_mu, T(v_coeff);
                                      ndrange=N)
    return nothing
end