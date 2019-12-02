function save_ovf(sim::AbstractSimGPU, fname::String; dataformat::String = "auto")
  T = _cuda_using_double.x ? Float64 : Float32

  if dataformat == "auto"
    dataformat = _cuda_using_double.x ? "Binary 8" : "Binary 4"
  end

  mesh = sim.mesh
  nx,ny,nz = mesh.nx,mesh.ny,mesh.nz
  nxyz = nx*ny*nz

  ovf = OVF2{T}()
  ovf.xnodes = mesh.nx
  ovf.ynodes = mesh.ny
  ovf.znodes = mesh.nz
  ovf.xstepsize = mesh.dx
  ovf.ystepsize = mesh.dy
  ovf.zstepsize = mesh.dz
  ovf.data_type = dataformat
  ovf.name = sim.name
  ovf.data = zeros(T,3*nxyz)
  copyto!(ovf.data, sim.spin)

  if !endswith(fname,".ovf")
      fname = fname* ".ovf"
  end
  io = open(fname, "w")
  write_OVF2_Header(io, ovf)
  write_OVF2_Data(io, ovf)
  hdr(io, "End", "Segment")
  close(io)
end


function read_ovf(sim::AbstractSimGPU, fname::String)
  T = _cuda_using_double.x ? Float64 : Float32
  ovf  = read_ovf(fname, T=T)
  nxyz = ovf.xnodes*ovf.ynodes*ovf.znodes
  if nxyz != sim.nxyz
      error("The ovf does not match the sim.mesh")
  end
  copyto!(sim.prespin, ovf.data)
  copyto!(sim.spin, ovf.data)
end