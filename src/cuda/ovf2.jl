function write_OVF2_Text(io::IOStream, sim::AbstractSimGPU)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz

    spin = zeros(Float32, 3, nx, ny, nz)
    copyto!(spin, sim.spin)

    b = reshape(spin, (3, nx, ny, nz))

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, string(b[1, i, j, k], " ", b[2, i, j, k], " ", b[3, i, j, k], "\n"))
    end

end