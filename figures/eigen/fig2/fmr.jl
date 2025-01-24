using MicroMagnetic
using DelimitedFiles
using LinearAlgebra
using CairoMakie

function setup(;m0=(0,0,1), H=(0,0,0))
    mesh = FDMesh(nx=1, ny=1, nz=1, dx=5e-9, dy=5e-9, dz=2e-9)
    sim = create_sim(mesh; H=H, m0=m0, Ms=8e5, Kc=2e4, demag=false)
    return sim
end

MicroMagnetic.set_precision(AbstractFloat)

function compute_frequency100(H0)
    H = (H0, 0, 0)
    sim = setup(m0=(1,0,0), H=H)
    B = build_matrix(sim, gamma=2.21e5)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function compute_frequency110(H0)
    H = (H0/sqrt(2), H0/sqrt(2),0)
    sim = setup(H=H, m0=(1,1,0))
    B = build_matrix(sim, gamma=2.21e5)
    println(B)
    return imag(eigvals(B)[2])/1e9/(2*pi)
end

function analytical100(Kc=2e4, Ms=8e5)
    H = range(0.8e5, 2e5; length=20)
    gamma = 2.21e5
    K = 4*Kc/(mu_0*Ms)
    return H*mu_0, gamma*(H .+ K)/1e9/(2*pi)
end

function analytical110(Kc=2e4, Ms=8e5)
    H = range(0.8e5, 2e5; length=20)
    gamma = 2.21e5
    K = Kc/(mu_0*Ms)
    freq = gamma*sqrt.(H.^2 .- 2K.*H .-8K^2)/1e9/(2*pi)
    return H*mu_0, freq
end

function plot_freq()
    freq110 = []
    freq100 = []

    Hs = range(0.8e5, 2e5; length=20)
    for H in Hs
        push!(freq100, compute_frequency100(H))
        push!(freq110, compute_frequency110(H))
    end

    fig = Figure(; size=(500, 360), fontsize=18)
    ax = Axis(fig[1, 1]; xlabel="H (T)", ylabel="Frequency (GHz)")

    s1 = scatter!(ax, Hs*mu_0, freq100; marker=:rect, markersize=10, strokewidth=1, alpha=0, color=:white, label="[100]")
    s1 = scatter!(ax, Hs*mu_0, freq110; markersize=10, strokewidth=1, alpha=0, color=:white, label="[110]")

    Hs, freq = analytical110()
    l1 = lines!(ax, Hs, freq; linestyle=:solid, color=:slateblue1, label="Analytical")

    Hs, freq = analytical100()
    l1 = lines!(ax, Hs, freq; linestyle=:solid, color=:slateblue1)

    axislegend(ax; position=(0.95, 0.05), labelsize=14)
    #axislegend(ax, [l1, l2, l3], [L"\beta=0", L"\beta=0.1", L"\beta=0.2"], position = :rt, orientation = :horizontal, labelsize=18)
    #axislegend(ax, [s1], [L"\text{JuMag}"], position=(0.78, 0.07), orientation = :horizontal, labelsize=16)

    save("fig2.pdf", fig)
    return nothing
end

plot_freq()