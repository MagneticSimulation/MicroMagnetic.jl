using MicroMagnetic

@using_gpu()

function relax_system_single(T; N=30, M=2000, m0=nothing)
    mesh = CubicMesh(; nx=N, ny=N, nz=N, pbc="xyz")
    sim = MonteCarlo(mesh; name="mc")
    
    add_exch(sim; J=300 * k_B)

    if m0 == nothing
        init_m0_random(sim)
        sim.T = 100000
        run_sim(sim; max_steps=10000, save_vtk_every=-1, save_m_every=-1)
    else
        init_m0(sim, m0)
    end

    sim.T = T
    run_sim(sim; max_steps=50000, save_vtk_every=-1, save_m_every=-1)

    ms = zeros(M)
    
    sim.T = T
    for i in 1:M
        run_sim(sim; max_steps=100, save_vtk_every=-1, save_m_every=-1)
        m = MicroMagnetic.average_m(sim)
        ms[i] = sqrt(m[1]^2 + m[2]^2 + m[3]^2)
    end
    
    m_average = sum(ms) / M
    m2_average = sum(ms.^2) / M

    return m_average, N^3*(m2_average-m_average^2)/T, Array(sim.spin)
end

function relax_system()
    f = open("assets/M_T.txt", "w")
    write(f, "#T(K)     m    chi\n")
    
    Ts = [600:-20:460; 455:-5:405; 400:-20:10]
    m0 = nothing
    for T in Ts
        println("Running for $T ...")
        m, chi, m0 = relax_system_single(T, m0=m0)
        write(f, "$T    $m   $chi\n")
    end
    return close(f)
end

relax_system()