module EnzymeExt

using Enzyme
using MicroMagnetic

function wrap_LLG(x::Array{Float64}, y::Array{Float64}, sim, gamma, alpha)
    N = sim.n_total
    Heff = zeros(Float64, 3*N)
    H = zeros(Float64, 3*N)

    # Compute the effective field due to all interactions
    for interaction in sim.interactions
        if !isa(interaction, MicroMagnetic.Zeeman)
            MicroMagnetic.effective_field(interaction, sim, x, 0.0, output=H)
            Heff .+= H
        else
            Heff .+= interaction.field
        end
    end
   
    for i = 1:N
        j = 3*(i-1) + 1
        Hx, Hy, Hz = Heff[j], Heff[j+1], Heff[j+2]
        mx, my, mz = x[j], x[j+1], x[j+2]
        y[j] = gamma*MicroMagnetic.cross_x(mx, my, mz, Hx, Hy, Hz)
        y[j+1] = gamma*MicroMagnetic.cross_y(mx, my, mz, Hx, Hy, Hz)
        y[j+2] = gamma*MicroMagnetic.cross_z(mx, my, mz, Hx, Hy, Hz)
    end
    return nothing
end

function MicroMagnetic.dynamic_matrix(sim; gamma=2.21e5, sparse=false, alpha=0.01)
    m0 = Array(sim.spin)
    N = length(m0)
    y  = similar(m0)
    dy = similar(m0)
    dx = similar(m0)
    A = zeros(N, N)

    for i = 1:N
        dx .= 0
        dx[i] = 1
        Enzyme.autodiff(set_runtime_activity(Forward), wrap_LLG, Duplicated(m0, dx), Duplicated(y, dy), Const(sim), Const(gamma), Const(alpha));
        A[:,i] = dy
    end

    N = sim.n_total
    C = zeros(2N, 2N)

    Rs = [MicroMagnetic.rotation_matrix(m0[3i-2], m0[3i-1], m0[3i]) for i=1:N]

    for i = 1:N
        for j = 1:N
            B = Rs[i]'*A[3i-2:3i, 3j-2:3j]*Rs[j]
            C[2i-1:2i, 2j-1:2j] .= B[1:2, 1:2]
        end
    end
    return C
end

function __init__()
end

end
