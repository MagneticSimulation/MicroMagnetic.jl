using Printf

function set_Ms(sim::MicroSimFE, Ms::NumberOrArray)
    sim.mu0_Ms .= mu_0 * Ms
    compute_L_Ms!(sim.L_mu, sim.mesh, sim.mu0_Ms)
    return true
end

function init_m0(sim::MicroSimFE, m0::TupleOrArrayOrFunction; norm=true)
    init_vector!(sim.prespin, sim.mesh, m0)
    if norm
        normalise(sim.prespin, sim.n_total)
    end

    if any(isnan, sim.prespin)
        error("NaN is given by the input m0!")
    end

    sim.spin[:] .= sim.prespin[:]
end

function average_m(sim::MicroSimFE)
    b = sim.spin .* sim.L_mu

    b = reshape(b, 3, sim.n_total)

    mx = 3 * sum(b[1, :]) / sum(sim.L_mu)
    my = 3 * sum(b[2, :]) / sum(sim.L_mu)
    mz = 3 * sum(b[3, :]) / sum(sim.L_mu)

    return (mx, my, mz)
end

function save_inp(sim::MicroSimFE, fname::String)
    mesh = sim.mesh

    #create the parent folder if necessary. TODO: copy this line to other places.
    mkpath(dirname(fname))

    f = open(fname, "w")

    data_number = 1 # FIX ME
    write(f, @sprintf("%d %d %d 0 0\n", mesh.number_nodes, mesh.number_cells, data_number))# data number
    for n in 1:(mesh.number_nodes)
        write(f,
              @sprintf("%d %0.12g %0.12g %0.12g\n", n, mesh.coordinates[1, n],
                       mesh.coordinates[2, n], mesh.coordinates[3, n]))
    end

    for c in 1:(mesh.number_cells)
        write(f,
              @sprintf("%d %d tet %d %d %d %d\n", c, 1, #FIX ME: material id
                       mesh.cell_verts[1, c], mesh.cell_verts[2, c], mesh.cell_verts[3, c],
                       mesh.cell_verts[4, c]))
    end

    write(f, @sprintf("%d", data_number))
    for i in 1:data_number
        write(f, " 3")
    end
    write(f, "\nm, none")

    m = reshape(sim.spin, 3, mesh.number_nodes)
    # FIX ME: to add more data
    for n in 1:(mesh.number_nodes)
        write(f, @sprintf("\n%d %0.12g %0.12g %0.12g", n, m[1, n], m[2, n], m[3, n]))
    end

    return close(f)
end

function save_inp_field(sim::MicroSimFE, fname::String, field::Array)
    mesh = sim.mesh

    #create the parent folder if necessary. TODO: copy this line to other places.
    mkpath(dirname(fname))

    f = open(fname, "w")

    data_number = 1 # FIX ME
    write(f, @sprintf("%d %d %d 0 0\n", mesh.number_nodes, mesh.number_cells, data_number))# data number
    for n in 1:(mesh.number_nodes)
        write(f,
              @sprintf("%d %0.12g %0.12g %0.12g\n", n, mesh.coordinates[1, n],
                       mesh.coordinates[2, n], mesh.coordinates[3, n]))
    end

    for c in 1:(mesh.number_cells)
        write(f,
              @sprintf("%d %d tet %d %d %d %d\n", c, 1, #FIX ME: material id
                       mesh.cell_verts[1, c], mesh.cell_verts[2, c], mesh.cell_verts[3, c],
                       mesh.cell_verts[4, c]))
    end

    write(f, @sprintf("%d", data_number))
    for i in 1:data_number
        write(f, " 3")
    end
    write(f, "\nm, none")

    m = reshape(field, 3, mesh.number_nodes)
    # FIX ME: to add more data
    for n in 1:(mesh.number_nodes)
        write(f, @sprintf("\n%d %0.12g %0.12g %0.12g", n, m[1, n], m[2, n], m[3, n]))
    end

    return close(f)
end