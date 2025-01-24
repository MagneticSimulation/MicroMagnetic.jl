using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using CairoMakie

MicroMagnetic.set_backend("cpu")

function setup(;H=(0,0,0), K=0.005, B=0)
    mesh = CubicMesh(nx=100, ny=1, nz=1)

    sim = Sim(mesh)
    set_mu_s(sim, 1.0)

    init_m0(sim, (0,0,1))

    J = 1
    add_exch(sim, J)
    add_anis(sim, K*J; axis=(0, 0, 1))
    add_zeeman(sim, H)
    add_exch_bq(sim, B)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency(;H0=0, K=0.005, B=0)
    H = (0, 0, H0)
    sim = setup(H=H, K=K, B=B)
    B = build_matrix(sim, gamma=1)
    all = imag(eigvals(B))
    sort!(all)
    return all[101:end]
end


function analytical(K, B, n)
    gamma = 1
    J = 1
    Ku = K*J
    mu_s = 1
    N=100
    freq = 2*gamma/mu_s*(K+(J+2B)*(1-cos(n*pi/N)))
    return freq
end

function analytical_freq(;K=0.005, B=0)
    N=100
    ns = [i for i in 1:N]
    freq = [analytical(K, B, n-1) for n in ns]
    return ns, freq
end

function plot_freq()

    fs = compute_frequency(B=0)
    fs2 = compute_frequency(B=0.1)
    
    fig = Figure(; size=(500, 360), fontsize=18)
    ax = Axis(fig[1, 1]; xlabel="Mode n", ylabel="Frequency")

    s1 = scatter!(ax, fs; marker=:rect, markersize=4, strokewidth=1, alpha=0, color=:white, label="B=0")
    s1 = scatter!(ax, fs2; markersize=4, strokewidth=1, alpha=0, color=:white, label="B=0.1")

    
    ns, freq = analytical_freq(K=0.005, B=0)
   
    l1 = lines!(ax, ns, freq; linestyle=:solid, color=:slateblue1, label="Analytical")

    ns, freq = analytical_freq(K=0.005, B=0.1)
    l1 = lines!(ax, ns, freq; linestyle=:solid, color=:slateblue1)
    
    axislegend(ax; position=(0.95, 0.05), labelsize=14)
    #axislegend(ax, [l1, l2, l3], [L"\beta=0", L"\beta=0.1", L"\beta=0.2"], position = :rt, orientation = :horizontal, labelsize=18)
    #axislegend(ax, [s1], [L"\text{JuMag}"], position=(0.78, 0.07), orientation = :horizontal, labelsize=16)

    save("fig3.pdf", fig)
    return nothing
end

plot_freq()