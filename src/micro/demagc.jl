function F(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    return atan(y * z / (x * r))
end

function G(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    return -log(abs(r + z))
end

function Nxx(x, y, z, b, c)
    D = F(x, y + b, z + c) - F(x, y - b, z + c) - F(x, y + b, z - c) + F(x, y - b, z - c)
    return D / (4 * pi)
end

function Nyx(x, y, z, b, c)
    D = G(x, y + b, z + c) - G(x, y - b, z + c) - G(x, y + b, z - c) + G(x, y - b, z - c)
    return D / (4 * pi)
end

function field_x_surface(x, y, z, b, c, sigma)
    Hx = -Nxx(x, y, z, b, c)*sigma
    Hy = -Nyx(x, y, z, b, c)*sigma
    Hz = -Nyx(x, z, y, c, b)*sigma
    return SVector(Hx, Hy, Hz)
end

function field_y_surface(x, y, z, a, c, sigma)
    Hy, Hz, Hx = field_x_surface(y, z, x, c, a, sigma)
    return SVector(Hx, Hy, Hz)
end

function field_z_surface(x, y, z, a, b, sigma)
    Hz, Hx, Hy = field_x_surface(z, x, y, a, b, sigma)
    return SVector(Hx, Hy, Hz)
end

function compute_demag_charge_field(sim, mesh::FDMesh, Ms=8e5, xminus=0, xplus=0, yminus=0, yplus=0, zminus=0, zplus=0)
    dx, dy, dz = mesh.dx*1e9, mesh.dy*1e9, mesh.dz*1e9
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    field = zeros(3, nx, ny, nz)
    Lx, Ly, Lz = nx*dx, ny*dy, nz*dz

    a, b, c = Lx/2, Ly/2, Lz/2

    xminus = xminus*Ms
    xplus = xplus*Ms
    yminus = yminus*Ms
    zminus = zminus*Ms
    
    yplus = yplus*Ms
    zplus = zplus*Ms

    if xminus != 0 || xplus != 0
        for k in 1:nz, j in 1:ny, i in 1:nx
            x = (i-0.5)*dx
            y = (j-0.5)*dy - b
            z = (k-0.5)*dz - c
            field[:, i, j, k] .+= field_x_surface(x, y, z, b, c, xminus) + field_x_surface(x-Lx, y, z, b, c, xplus)
        end
    end

    if yminus != 0 || yplus != 0
        for k in 1:nz, j in 1:ny, i in 1:nx
            x = (i-0.5)*dx - a
            y = (j-0.5)*dy
            z = (k-0.5)*dz - c
            field[:, i, j, k] .+= field_y_surface(x, y, z, a, c, yminus) + field_y_surface(x, y-Ly, z, a, c, yplus)
        end
    end
    
    if zminus != 0 || zplus != 0
        for k in 1:nz, j in 1:ny, i in 1:nx
            x = (i-0.5)*dx - a
            y = (j-0.5)*dy - b
            z = (k-0.5)*dz 
            field[:, i, j, k] .+= field_z_surface(x, y, z, a, b, zminus) + field_z_surface(x, y, z-Lz, a, b, zplus)
        end
    end
    
    field = reshape(field, 3*sim.n_total)
    return field
end