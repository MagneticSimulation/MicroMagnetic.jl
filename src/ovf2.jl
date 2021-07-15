mutable struct OVF2{T<:AbstractFloat}
    xnodes::Int64
    ynodes::Int64
    znodes::Int64
    xstepsize::Float64
    ystepsize::Float64
    zstepsize::Float64
    data_type::String
    name::String
    data::Array{T, 1}
    OVF2{T}() where {T<:AbstractFloat} = new()
end

"""
  save_ovf(sim::AbstractSim, fname::String; dataformat::String = "text")

Save spins in format of ovf, which can be viewed by Muview. Dataformat could be "text", "binary4" or "binary8".

For example:

    ```julia
        save_ovf(sim, "m0")
    ```
"""
function save_ovf(sim::AbstractSim, fname::String; dataformat::String = "Binary 8")
    mesh = sim.mesh
    nxyz = mesh.nx*mesh.ny*mesh.nz

    ovf = OVF2{Float64}()
    ovf.xnodes = mesh.nx
    ovf.ynodes = mesh.ny
    ovf.znodes = mesh.nz
    ovf.xstepsize = mesh.dx
    ovf.ystepsize = mesh.dy
    ovf.zstepsize = mesh.dz
    ovf.data_type = dataformat
    ovf.name = sim.name
    ovf.data = zeros(Float64,3*nxyz)
    copyto!(ovf.data, sim.spin)

    save_ovf(ovf, fname, dataformat=dataformat)
end

function save_ovf(ovf, fname; dataformat::String = "Binary 8")
    if !endswith(fname,".ovf")
        fname = fname* ".ovf"
    end

    io = open(fname, "w")
    ovf.data_type = dataformat
    write_OVF2_Header(io, ovf)
    write_OVF2_Data(io, ovf)
    hdr(io, "End", "Segment")
    close(io)
end

function hdr(io::IOStream, string1::Any, string2::Any)
    write(io, string("# ", string(string1), ": ", string(string2), "\n"))
end

function write_OVF2_Header(io::IOStream, ovf::OVF2)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    xyz = zeros(Float32, 3, nx, ny, nz)

    dx, dy, dz=ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    write(io, "# OOMMF OVF 2.0\n")
    hdr(io, "Segment count", "1")
    hdr(io, "Begin", "Segment")
    hdr(io, "Begin", "Header")

    hdr(io, "Title", ovf.name)
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

function write_OVF2_Data(io::IOStream, ovf::OVF2)
    dataformat = lowercase(ovf.data_type)

    if dataformat == "text"
        hdr(io, "Begin", "Data Text")
        write_OVF2_Text(io, ovf)
        hdr(io, "End", "Data Text")
    elseif dataformat == "binary 4" || dataformat == "binary"
        hdr(io, "Begin", "Data Binary 4")
        write_OVF2_Binary4(io, ovf)
        write(io, "\n")
        hdr(io, "End", "Data Binary 4")
    elseif dataformat == "binary 8"
        hdr(io, "Begin", "Data Binary 8")
        write_OVF2_Binary8(io, ovf)
        write(io, "\n")
        hdr(io, "End", "Data Binary 8")
    else
        @info "Data format error!"
    end
end

function write_OVF2_Text(io::IOStream, ovf::OVF2)

    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes

    b = reshape(ovf.data, (3, nx, ny, nz))

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, string(b[1, i, j, k], " ", b[2, i, j, k], " ", b[3, i, j, k], "\n"))
    end

end

function write_OVF2_Binary4(io::IOStream, ovf::OVF2)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    nxyz = nx*ny*nz
    spin = zeros(Float32, 3*nxyz)
    copyto!(spin, ovf.data)

    b = reshape(spin, (3, nx, ny, nz))

    write(io, Float32(1234567.0))   ##OOMMF requires this number to be first to check the format

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, b[1, i, j, k])
        write(io, b[2, i, j, k])
        write(io, b[3, i, j, k])
    end

end

function write_OVF2_Binary8(io::IOStream, ovf::OVF2)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes

    b = reshape(ovf.data, (3, nx, ny, nz))

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
    spin = zeros(Float32, 3*nxyz)
    if read(io,Float32) == 1234567.0
	  read!(io, spin)
    end
    ovf.data = convert(Array{Float64,1}, spin)
    return nothing
end

function read_OVF2_Binary8(io::IOStream, ovf::OVF2)
    ovf.data_type = "Binary 8"
    nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*nxyz)
    if read(io, Float64) == 123456789012345.0
      read!(io, spin)
    end
    ovf.data = spin
    return nothing
end

function read_OVF2_Text(io::IOStream, ovf::OVF2)
    ovf.data_type = "Text"
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
function read_ovf(fname::String; T::DataType=Float64)
    if !endswith(fname, ".ovf")
        fname = fname*".ovf"
    end

    ovf = OVF2{T}()

    io = open(fname, "r")
    if !startswith(readline(io), "# OOMMF OVF 2")
        error("Input is not a OVF2 file.")
    end

    for line in eachline(io)
        if startswith(line, "# Title:")
            ovf.name = line[10:end]
        elseif startswith(line, "# xnodes:")
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


function ovf2npy(ovf_name, npy_name=nothing; box=nothing)
    if npy_name == nothing
        npy_name = endswith(ovf_name, ".ovf") ? ovf_name[1:end-4]*".npy" : ovf_name*".npy"
    end
    ovf = read_ovf(ovf_name)
    nx,ny,nz = ovf.xnodes,ovf.ynodes,ovf.znodes
    if box != nothing
        nx = box[2]-box[1]+1
        ny = box[4]-box[3]+1
        nz = box[6]-box[5]+1
    end

    data = reshape(ovf.data, (3, ovf.xnodes, ovf.ynodes, ovf.znodes))
    b = data
    if box != nothing
        b = data[:, box[1]:box[2], box[3]:box[4], box[5]:box[6]]
    end

    b = reshape(b, 3*nx*ny*nz)

    npzwrite(npy_name, b)
end