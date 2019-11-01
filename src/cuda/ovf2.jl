function write_OVF2_Data(io::IOStream, sim::AbstractSimGPU, dataformat::String)
    if dataformat == "text"
        hdr(io, "Begin", "Data " * dataformat)
        write_OVF2_Text(io, sim)
    elseif dataformat == "binary" || dataformat == "binary4" || dataformat == "binary8"
        write_OVF2_Binary(io, sim, dataformat)
    else
        @info "Data format error!"
    end
end

function write_OVF2_Text(io::IOStream, sim::AbstractSimGPU)
    Float = _cuda_using_double.x ? Float64 : Float32
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    nxyz = nx*ny*nz
    spin = zeros(Float, 3*nxyz)

    copyto!(spin, sim.spin)

    b = reshape(spin, (3, nx, ny, nz))

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, string(b[1, i, j, k], " ", b[2, i, j, k], " ", b[3, i, j, k], "\n"))
    end

end

function write_OVF2_Binary(io::IOStream, sim::AbstractSimGPU,dataformat)
    Float = _cuda_using_double.x ? Float64 : Float32
    binary = _cuda_using_double.x ? "Binary 8" : "Binary 4"
    hdr(io, "Begin", "Data "  * binary)

    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    nxyz = nx*ny*nz
    spin = zeros(Float, 3*nxyz)

    copyto!(spin, sim.spin)

    b = reshape(spin, (3, nx, ny, nz))
    checknumber = _cuda_using_double.x ? Float(123456789012345.0) : Float(1234567.0)

    write(io, checknumber)   ##OOMMF requires this number to be first to check the format

    for k = 1:nz, j = 1:ny, i = 1:nx
        write(io, Float(b[1, i, j, k]))
        write(io, Float(b[2, i, j, k]))
        write(io, Float(b[3, i, j, k]))
    end
    write(io, "\n")
    hdr(io, "End", "Data " * binary)

end


function read_OVF2_Binary4(io::IOStream, sim::AbstractSimGPU)
    Float = _cuda_using_double.x ? Float64 : Float32
    nxyz = sim.nxyz
    spin = zeros(Float, 3*nxyz)
    if read(io,Float32) == Float32(1234567.0)
      for i = 1:3*nxyz
        spin[i] = Float(read(io,Float32))
      end
    else
        @info "Data format error!"
    end
    copyto!(sim.prespin, spin)
    copyto!(sim.spin, spin)
end

function read_OVF2_Binary8(io::IOStream, sim::AbstractSimGPU)
    Float = _cuda_using_double.x ? Float64 : Float32
    nxyz = sim.nxyz
    spin = zeros(Float, 3*nxyz)
    if read(io, Float64) == Float64(123456789012345.0)
      for i = 1:3*nxyz
        spin[i] = Float(read(io,Float64))
      end
    else
        @info "Data format error in read_OVF2_Binary8"
    end
    copyto!(sim.prespin, spin)
    copyto!(sim.spin, spin)
end

function read_OVF2_Binary4(io::IOStream, sim::MonteCarlo)
    Float = _cuda_using_double.x ? Float64 : Float32
    nxyz = sim.nxyz
    spin = zeros(Float, 3*nxyz)
    if read(io,Float32) == Float32(1234567.0)
      for i = 1:3*nxyz
        spin[i] = Float(read(io,Float32))
      end
    else
        @info "Data format error!"
    end

    copyto!(sim.spin, spin)
end

function read_OVF2_Binary8(io::IOStream, sim::MonteCarlo)
    Float = _cuda_using_double.x ? Float64 : Float32
    nxyz = sim.nxyz
    spin = zeros(Float, 3*nxyz)
    if read(io, Float64) == Float64(123456789012345.0)
      for i = 1:3*nxyz
        spin[i] = Float(read(io,Float64))
      end
    else
        @info "Data format error in read_OVF2_Binary8"
    end

    copyto!(sim.spin, spin)
end