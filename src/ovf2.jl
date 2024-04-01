mutable struct OVF2{T<:AbstractFloat}
    xnodes::Int64
    ynodes::Int64
    znodes::Int64
    xstepsize::Float64
    ystepsize::Float64
    zstepsize::Float64
    type::DataType
    name::String
    data::Array{T, 1}
    OVF2{T}() where {T<:AbstractFloat} = new()
end

"""
    save_ovf(sim::AbstractSim, fname::String; type::DataType = Float64)

Save spins by ovf, which can be viewed by Muview. 

Parameters:

    Sim : Sim struct whose spin to be saved.

    fname : Save file name.

Optional:

    type : Data type of ovf2 file. Can be chosen from Float32, Float64 or String.

For example:

    ```julia
        save_ovf(sim, "ovf_example")
    ```

Or to specify a certain data type:

    ```julia
        save_ovf(sim, "ovf_example", type = String)
    ```
"""
function save_ovf(sim::AbstractSim, fname::String; type::DataType = Float64)
    mesh = sim.mesh
    n_total = mesh.nx*mesh.ny*mesh.nz

    ovf = OVF2{Float64}()
    ovf.xnodes = mesh.nx
    ovf.ynodes = mesh.ny
    ovf.znodes = mesh.nz
    ovf.xstepsize = mesh.dx
    ovf.ystepsize = mesh.dy
    ovf.zstepsize = mesh.dz
    ovf.type = type
    ovf.name = sim.name
    ovf.data = zeros(Float64,3*n_total)
    copyto!(ovf.data, sim.spin)

    save_ovf(ovf, fname)
end

"""
    save_ovf(ovf::OVF2, fname::String; type::DataType = ovf.type)

Save an ovf2 file by an OVF2 struct.
"""
#FIXME: type is not used.
function save_ovf(ovf::OVF2, fname::String; type::DataType = ovf.type)
    if !endswith(fname,".ovf")
        fname = fname* ".ovf"
    end
    io = open(fname, "w")
    write_ovf2_header(io, ovf)
    write_ovf2_data(io, ovf)
    hdr(io, "End", "Segment")
    close(io)
end

function hdr(io::IOStream, string1::Any, string2::Any)
    write(io, string("# ", string(string1), ": ", string(string2), "\n"))
end

function write_ovf2_header(io::IOStream, ovf::OVF2)
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

function write_ovf2_data(io::IOStream, ovf::OVF2)
    type = ovf.type

    if type == String
        hdr(io, "Begin", "Data Text")
        write_ovf2_string(io, ovf)
        hdr(io, "End", "Data Text")
    elseif type == Float32
        hdr(io, "Begin", "Data Binary 4")
        write_ovf2_float32(io, ovf)
        write(io, "\n")
        hdr(io, "End", "Data Binary 4")
    elseif type == Float64
        hdr(io, "Begin", "Data Binary 8")
        write_ovf2_float64(io, ovf)
        write(io, "\n")
        hdr(io, "End", "Data Binary 8")
    else
        @info "OVF type error!"
    end
end

function write_ovf2_string(io::IOStream, ovf::OVF2)

    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes

    b = reshape(ovf.data, (3, nx, ny, nz))

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, string(b[1, i, j, k], " ", b[2, i, j, k], " ", b[3, i, j, k], "\n"))
    end

end

function write_ovf2_float32(io::IOStream, ovf::OVF2)
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    n_total = nx*ny*nz
    spin = zeros(Float32, 3*n_total)
    copyto!(spin, ovf.data)

    b = reshape(spin, (3, nx, ny, nz))

    write(io, Float32(1234567.0))   ##OOMMF requires this number to be first to check the format

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, b[1, i, j, k])
        write(io, b[2, i, j, k])
        write(io, b[3, i, j, k])
    end

end

function write_ovf2_float64(io::IOStream, ovf::OVF2)
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
    n_total = ovf.xnodes*ovf.ynodes*ovf.znodes
    if n_total != sim.n_total
        error("The ovf does not match the sim.mesh")
    end
    copyto!(sim.prespin, ovf.data)
    copyto!(sim.spin, ovf.data)
end

function read_ovf2_float32(io::IOStream, ovf::OVF2)
    ovf.type = Float32
    n_total = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float32, 3*n_total)
    if read(io,Float32) == 1234567.0
	  read!(io, spin)
    end
    ovf.data = convert(Array{Float64,1}, spin)
    return nothing
end

function read_ovf2_float64(io::IOStream, ovf::OVF2)
    ovf.type = Float64
    n_total = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*n_total)
    if read(io, Float64) == 123456789012345.0
      read!(io, spin)
    end
    ovf.data = spin
    return nothing
end

function read_ovf2_string(io::IOStream, ovf::OVF2)
    ovf.type = String
    n_total = ovf.xnodes*ovf.ynodes*ovf.znodes
    spin = zeros(Float64, 3*n_total)
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
                read_ovf2_float64(io, ovf)
                close(io)
            elseif  line[15:end] == "Binary 4"
                read_ovf2_float32(io, ovf)
                close(io)
            elseif lowercase(line[15:end]) == "text"
                read_ovf2_string(io, ovf)
                close(io)
            else
                error("Data format error!")
            end
        end
    end
    return ovf
end

"""
    mag2ovf(m::Array{T, 1}, nx::Int, ny::Int, nz::Int; 
    dx::Float64 = 1e-9, dy::Float64=1e-9, dz::Float64=1e-9, 
    name::String = "mag", type::DataType = Float64) where T <: AbstractFloat


    Save a magnetization array into an ovf file.
    
Required Parameters:
        m : 1d m array to be converted
        nx, ny, nz: shape of mesh

Optional Parameters:
        fname : output file name
        type : output data type which can be chosen from Float64, Float32 or String
"""

function mag2ovf(m::Array{T, 1}, nx::Int, ny::Int, nz::Int; 
    dx::Float64 = 1e-9, dy::Float64 = 1e-9, dz::Float64 = 1e-9, 
    fname::String = "mag", type::DataType = Float64) where T <: AbstractFloat

    ovf = JuMag.OVF2{Float64}()
    ovf.xnodes = nx
    ovf.ynodes = ny
    ovf.znodes = nz
    ovf.xstepsize = dx
    ovf.ystepsize = dy
    ovf.zstepsize = dz
    ovf.type = type
    ovf.name = fname

    n_total = nx*ny*nz
    ovf.data = zeros(Float64, 3*n_total)
    copyto!(ovf.data, m)

    save_ovf(ovf, fname)
end