mutable struct OVF2
    xnodes::Int64
    ynodes::Int64
    znodes::Int64
    xstepsize::Float64
    ystepsize::Float64
    zstepsize::Float64
    data_type::String
    data::Array{Float64, 1}
    OVF2() = new()
end


"""
  save_ovf(sim::AbstractSim, fname::String; dataformat::String = "text")

Save spins in format of ovf, which can be viewed by Muview. Dataformat could be "text", "binary4" or "binary8".

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

    dx, dy, dz=mesh.dx, mesh.dy, mesh.dz

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
    #hdr(io, "Desc", string("Total simulation time: ", sim.saver.t, " s")) // sim.saver may not be initialized

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
    elseif dataformat == "binary4" || dataformat == "binary"
        hdr(io, "Begin", "Data Binary 4")
        write_OVF2_Binary4(io, sim)
        write(io, "\n")
        hdr(io, "End", "Data Binary 4")
    elseif dataformat == "binary8"
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
    read_ovf(sim, fname)

    Initialize sim with an ovf file named of "fname.ovf".
"""
function read_ovf(sim::AbstractSim, fname::String)
    ovf  = read_ovf(fname)
    nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
    if nxyz != sim.nxyz
        error("The ovf does not match the sim.mesh")
    end
    copyto!(sim.prespin, ovf.data)
    copyto!(sim.spin, ovf.data)
end

function read_OVF2_Binary4(io::IOStream, ovf::OVF2)
    ovf.data_type = "Binary 4"
    nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*nxyz)
    if read(io,Float32) == 1234567.0
      for i = 1:3*nxyz
        spin[i] = Float64(read(io, Float32))
      end
    end
    ovf.data = spin
    return nothing
end

function read_OVF2_Binary8(io::IOStream, ovf::OVF2)
    ovf.data_type = "Binary 8"
    nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*nxyz)
    if read(io, Float64) == 123456789012345.0
      for i = 1:3*nxyz
        spin[i] = read(io, Float64)
      end
    end
    ovf.data = spin
    return nothing
end

function read_OVF2_Text(io::IOStream, ovf::OVF2)
    ovf.data_type = "text"
    nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*nxyz)
    i = 0
    for line in eachline(io)
        if startswith(line, "#")
            break
        end
        for s in split(line)
            i += 1
            spin[i] = parse(Float64, s)
        end
    end
    ovf.data = spin
    return nothing
end


"""
    read_ovf(fname)

    Load ovf file as OVF2 Type, where spin is stored in OVF2.data
"""
function read_ovf(fname::String)
    if !endswith(fname, ".ovf")
        fname = fname*".ovf"
    end

    ovf = OVF2()

    io = open(fname, "r")
    if !startswith(readline(io), "# OOMMF OVF 2")
        error("Input is not a OVF2 file.")
    end

    for line in eachline(io)
        if startswith(line, "# xnodes:")
            ovf.xnodes = parse(Int64, line[10:end])
        elseif startswith(line, "# ynodes:")
            ovf.ynodes = parse(Int64, line[10:end])
        elseif startswith(line, "# znodes:")
            ovf.znodes = parse(Int64, line[10:end])
        elseif startswith(line, "# xstepsize:")
            ovf.xstepsize = parse(Float64, line[13:end])
        elseif startswith(line, "# ystepsize:")
            ovf.ystepsize = parse(Float64, line[13:end])
        elseif startswith(line, "# zstepsize:")
            ovf.zstepsize = parse(Float64, line[13:end])
        elseif startswith(line, "# Begin: Data")
            if line[15:end] == "Binary 8"
                read_OVF2_Binary8(io, ovf)
                close(io)
            elseif  line[15:end] == "Binary 4"
                read_OVF2_Binary4(io, ovf)
                close(io)
            elseif lowercase(line[15:end]) == "text"
                read_OVF2_Text(io, ovf)
                close(io)
            else
                error("Data format error!")
            end
        end
    end
    return ovf
end
