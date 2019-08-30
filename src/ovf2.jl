"""
  save_ovf(sim::AbstractSim, fname::String; dataformat::String = "text")

Save spins in format of ovf, which can be viewed by Muview. Dataformat can be "text" or "Binary4"
For example:

```julia
   save_ovf(sim, "m0")
```

"""
function save_ovf(sim::AbstractSim, fname::String; dataformat::String = "binary")
    io = open(fname * ".ovf", "w")
    write_OVF2_Header(io, sim)
    write_OVF2_Data(io, sim, dataformat)
    hdr(io, "End", "Segment")
    close(io)
end

function hdr(io::IOStream, string1::Any, string2::Any)
    write(io, string("# ", string(string1), ": ", string(string2), "\n"))
end

function write_OVF2_Header(io::IOStream, sim::AbstractSim)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    xyz = zeros(Float32, 3, nx+1, ny+1, nz+1)
    dx, dy, dz = 1.0,1.0,1.0

    if isa(mesh, CubicMesh)
      dx, dy, dz=mesh.a, mesh.a, mesh.a
    else
      dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
    end

    write(io, "# OOMMF OVF 2.0\n")
    hdr(io, "Segment count", "1")
    hdr(io, "Begin", "Segment")
    hdr(io, "Begin", "Header")

    hdr(io, "Title", sim.name)
    hdr(io, "meshtype", "rectangular")
    hdr(io, "meshunit", "m")

	hdr(io, "xmin", 0)
	hdr(io, "ymin", 0)
	hdr(io, "zmin", 0)

	hdr(io, "xmax", nx*dx)
	hdr(io, "ymax", ny*dy)
	hdr(io, "zmax", nz*dz)

    ##hdr(io, "valuedim", length(sim.interactions))
    hdr(io, "valuedim", 3)

    labels, units = "m_x m_y m_z", "1 1 1"
    #= for i in sim.interactions
        labels = labels * " " * "i.name"
        units = units * " " * "J"
    end =#
    hdr(io, "valuelabels", labels)
    hdr(io, "valueunits", units)

	## We don't really have stages
	##fmt.Fprintln(io, "# Desc: Stage simulation time: ", meta.TimeStep, " s") // TODO
	hdr(io, "Desc",string("Total simulation time: ", sim.saver.t, " s"))

    hdr(io, "xbase", dx/2)
    hdr(io, "ybase", dy/2)
    hdr(io, "zbase", dz/2)
    hdr(io, "xnodes", nx)
    hdr(io, "ynodes", ny)
    hdr(io, "znodes", nz)
    hdr(io, "xstepsize", dx)
    hdr(io, "ystepsize", dy)
    hdr(io, "zstepsize", dz)
    hdr(io, "End", "Header")
end

function write_OVF2_Data(io::IOStream, sim::AbstractSim, dataformat::String)
    if dataformat == "text"
        hdr(io, "Begin", "Data " * dataformat)
        write_OVF2_Text(io, sim)
        hdr(io, "End", "Data " * dataformat)
    elseif dataformat == "binary4"
        hdr(io, "Begin", "Data Binary 4")
        write_OVF2_Binary4(io, sim)
        write(io, "\n")
        hdr(io, "End", "Data Binary 4")
    elseif dataformat == "binary8" || dataformat == "binary"
        hdr(io, "Begin", "Data Binary 8")
        write_OVF2_Binary8(io, sim)
        write(io, "\n")
        hdr(io, "End", "Data Binary 8")
    else
        @info "Data format error!"
    end
end

function write_OVF2_Text(io::IOStream, sim::AbstractSim)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz

    b = reshape(sim.spin, (3, nx, ny, nz))

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, string(b[1, i, j, k], " ", b[2, i, j, k], " ", b[3, i, j, k], "\n"))
    end

end

function write_OVF2_Binary4(io::IOStream, sim::AbstractSim)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    nxyz = nx*ny*nz
    spin = zeros(Float32, 3*nxyz)
    copyto!(spin, sim.spin)

    b = reshape(spin, (3, nx, ny, nz))

    write(io, Float32(1234567.0))   ##OOMMF requires this number to be first to check the format

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, b[1, i, j, k])
        write(io, b[2, i, j, k])
        write(io, b[3, i, j, k])
    end

end

function write_OVF2_Binary8(io::IOStream, sim::AbstractSim)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz

    b = reshape(sim.spin, (3, nx, ny, nz))

    write(io, Float64(123456789012345.0))   ##OOMMF requires this number to be first to check the format

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, Float64(b[1, i, j, k]))
        write(io, Float64(b[2, i, j, k]))
        write(io, Float64(b[3, i, j, k]))
    end

end
"""
read_ovf(sim. fname)
Initialize sim with an ovf file named of "fname.ovf".

"""

function read_ovf(sim::AbstractSim,fname::String)
    nxyz = sim.nxyz

    io = open(fname * ".ovf", "r")
    n = countlines(io)
    seek(io,0)
    for i = 1:n
        f = readline(io)
        if startswith(f, "# Begin: Data")
            if f[15:end] == "Binary 8"
                read_OVF2_Binary8(io,sim)
                close(io)
                return nothing
            elseif  f[15:end] == "Binary 4"
                read_OVF2_Binary4(io,sim)
                close(io)
                return nothing
            elseif  f[15:end] == "text"
                read_OVF2_Text(io,sim)
                close(io)
                return nothing
            else
                @info "Data format error!"
                return nothing
            end
        end
    end

end

function read_OVF2_Binary4(io::IOStream, sim::AbstractSim)
    nxyz = sim.nxyz
    spin = zeros(Float64, 3*nxyz)
    if read(io,Float32) == 1234567.0
      for i = 1:3*nxyz
        spin[i] = Float64(read(io, Float32))
      end
    else
        @info "Data format error!"
    end

    copyto!(sim.prespin, spin)
    copyto!(sim.spin, spin)
end

function read_OVF2_Binary8(io::IOStream, sim::AbstractSim)
    nxyz = sim.nxyz
    spin = zeros(Float64, 3*nxyz)
    if read(io, Float64) == 123456789012345.0
      for i = 1:3*nxyz
        spin[i] = read(io,Float64)
      end
    else
        @info "Data format error! ppp"
    end

    copyto!(sim.prespin, spin)
    copyto!(sim.spin, spin)
end

function read_OVF2_Text(io::IOStream, sim::AbstractSim)
    nxyz = sim.nxyz
    spin = zeros(Float64, 3*nxyz)
    for i = 1:3*nxyz
      spin[i] = Float64(read(io,Float64))
    end

    copyto!(sim.prespin, spin)
    copyto!(sim.spin, spin)
end
