using JuMag
using Test

Nx = 50

function m0_fun(i, j, k, dx, dy, dz)
    L = Nx * dx
    x = i * dx
    return sin(2 * pi * x / L), sin(2 * pi * x / L + 1.2), sin(2 * pi * x / L + 2.3)
end

function test_exch_scalar()
    
    mesh = FDMesh(dx=2e-9, nx=Nx, ny=1, nz=1, pbc="x")
    @test mesh.dx == 2e-9
    @test mesh.nx == Nx
    ngbs = Array(mesh.ngbs)
    @test ngbs[1, 1] == Nx
    @test ngbs[1, Nx] == Nx - 1
    @test ngbs[2, 1] == 2
    @test ngbs[2, Nx] == 1

    Ms = 8.6e5
    A = 1.3e-11

    sim = Sim(mesh)
    set_Ms(sim, Ms)
    init_m0(sim, m0_fun, norm=false)
    exch = add_exch(sim, A)

    JuMag.effective_field(sim, sim.spin, 0.0)

    if isa(sim.spin, Array)

      f1 = Array(exch.field)

      JuMag.effective_field_debug(exch, sim, sim.spin, 0.0)

      @test isapprox(f1, exch.field, atol=1e-10)
    end
    
    xs = (1:Nx) * 2e-9
    mu0 = 4 * pi * 1e-7
    L = Nx * 2e-9
    expected_x = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs)
    expected_y = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs .+ 1.2)
    expected_z = -2 * A / (mu0 * Ms) * (2 * pi / L)^2 * sin.((2 * pi / L) .* xs .+ 2.3)

    b = reshape(Array(exch.field), 3, sim.n_total)
    @test (maximum(b[1, :] .- expected_x) ./ Ms < 2e-4)
    @test (maximum(b[2, :] .- expected_y) ./ Ms < 2e-4)
    @test (maximum(b[3, :] .- expected_z) ./ Ms < 2e-4)


  Delta = 2e-9
  mesh = FDMesh(dx=2e-9, dz=Delta, nx=1, ny=1, nz=3, pbc="x")
  sim = Sim(mesh)
  set_Ms(sim, Ms)

  sigma = 1e-5
  init_m0(sim, (0.6, 0.8, 0))
  #r = add_exch_rkky(sim, sigma)

  #JuMag.effective_field(sim, sim.spin, 0.0)
  #b = reshape(sim.field, 3, sim.n_total)
  
  #fx = sigma / Delta / (mu0 * Ms) * 0.6
  #fy = sigma / Delta / (mu0 * Ms) * 0.8
  #println(fx - b[1, 1])
  #@test fx - b[1, 1] == 0.0


end










mesh = FDMesh(nx=3, ny=3, nz=3)
function m(i, j, k, dx, dy, dz)
  return (i^2, j + 1.0, k * j)
end

function test_exch_vector_x(mesh)
  function Ms_x(i, j, k, dx, dy, dz)
    if j == 2 && k == 2
      return 1e5
    else
      return 0
    end
  end
  sim1 = Sim(mesh)
  sim2 = Sim(mesh)
  set_Ms(sim1, Ms_x)
  set_Ms(sim2, Ms_x)
  add_exch(sim1, 1e-12)
  add_exch_vector(sim2, (1e-12, 0, 0))
  init_m0(sim1, m)
  init_m0(sim2, m)
  JuMag.effective_field(sim1, sim1.spin, 0.0)
  JuMag.effective_field(sim2, sim2.spin, 0.0)
  for i = 1:3*27
    @test sim1.field[i] - sim2.field[i] == 0.0
  end
end
function test_exch_vector_x(mesh)
  function Ms_x(i, j, k, dx, dy, dz)
    if j == 2 && k == 2
      return 1e5
    else
      return 0
    end
  end
  sim1 = Sim(mesh)
  sim2 = Sim(mesh)
  set_Ms(sim1, Ms_x)
  set_Ms(sim2, Ms_x)
  add_exch(sim1, 1e-12)
  add_exch_vector(sim2, (1e-12, 0, 0))
  init_m0(sim1, m)
  init_m0(sim2, m)
  JuMag.effective_field(sim1, sim1.spin, 0.0)
  JuMag.effective_field(sim2, sim2.spin, 0.0)
  for i = 1:3*27
    @test sim1.field[i] - sim2.field[i] == 0.0
  end
end
function test_exch_vector_y(mesh)
  function Ms_y(i, j, k, dx, dy, dz)
    if i == 2 && k == 2
      return 1e5
    else
      return 0
    end
  end
  sim1 = Sim(mesh)
  sim2 = Sim(mesh)
  set_Ms(sim1, Ms_y)
  set_Ms(sim2, Ms_y)
  add_exch(sim1, 1e-12)
  add_exch_vector(sim2, (0, 1e-12, 0))
  init_m0(sim1, m)
  init_m0(sim2, m)
  JuMag.effective_field(sim1, sim1.spin, 0.0)
  JuMag.effective_field(sim2, sim2.spin, 0.0)
  for i = 1:3*27
    @test sim1.field[i] - sim2.field[i] == 0.0
  end
end
function test_exch_vector_z(mesh)
  function Ms_z(i, j, k, dx, dy, dz)
    if i == 2 && j == 2
      return 1e5
    else
      return 0
    end
  end
  sim1 = Sim(mesh)
  sim2 = Sim(mesh)
  set_Ms(sim1, Ms_z)
  set_Ms(sim2, Ms_z)
  add_exch(sim1, 1e-12)
  add_exch_vector(sim2, (0, 0, 1e-12))
  init_m0(sim1, m)
  init_m0(sim2, m)
  JuMag.effective_field(sim1, sim1.spin, 0.0)
  JuMag.effective_field(sim2, sim2.spin, 0.0)
  for i = 1:3*27
    @test sim1.field[i] - sim2.field[i] == 0.0
  end
end
function test_exch_vector_all(mesh)
  function Ms_all(i, j, k, dx, dy, dz)
    return i^2 + j * 1e3 + k
  end
  sim1 = Sim(mesh)
  sim2 = Sim(mesh)
  set_Ms(sim1, Ms_all)
  set_Ms(sim2, Ms_all)
  add_exch(sim1, 1e-12)
  add_exch_vector(sim2, (1e-12, 1e-12, 1e-12))
  init_m0(sim1, m)
  init_m0(sim2, m)
  JuMag.effective_field(sim1, sim1.spin, 0.0)
  JuMag.effective_field(sim2, sim2.spin, 0.0)
  for i = 1:3*27
    @test sim1.field[i] - sim2.field[i] == 0.0
  end
end

#test_exch_vector_x(mesh)
#test_exch_vector_y(mesh)
#test_exch_vector_z(mesh)
#test_exch_vector_all(mesh)

test_exch_scalar()