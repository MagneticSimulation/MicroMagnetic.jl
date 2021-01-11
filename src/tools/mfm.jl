#using PyCall

function divergence(vector_field::Array{T,1}, mesh::Mesh) where T<:AbstractFloat
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    divM = zeros(T, nx, ny, nz)
    b = reshape(vector_field,(3,:))

    divM .+= divergence_x(b[1,:],mesh)
    divM .+= divergence_y(b[2,:],mesh)
    divM .+= divergence_z(b[3,:],mesh)

    return divM
end

function divergence_x(scalar_field::Array{T,1}, mesh) where T<:AbstractFloat
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx = mesh.dx
    div_x = zeros(T, nx, ny, nz)

    for i=1:nx,j=1:ny,k=1:nz
        id = index(i, j, k, nx, ny, nz)
        id1,id2 = index(i-1,j,k, nx,ny,nz),index(i+1,j,k, nx,ny,nz)
        if id1 > 0 && id2 >0
            div_x[i,j,k] = (scalar_field[id2] - scalar_field[id1])/(2*dx)
        end
        if id1 < 0
            div_x[i,j,k] = (scalar_field[id2] - scalar_field[id])/(dx)
        end

        if id2 < 0
            div_x[i,j,k] = (scalar_field[id] - scalar_field[id1])/(dx)
        end

    end

    return div_x
end

function divergence_y(scalar_field::Array{T,1}, mesh) where T<:AbstractFloat
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dy = mesh.dy
    div_y = zeros(T, nx, ny, nz)

    for i=1:nx,j=1:ny,k=1:nz
        id = index(i, j, k, nx, ny, nz)
        id1,id2 = index(i,j-1,k, nx,ny,nz),index(i,j+1,k, nx,ny,nz)
        if id1 > 0 && id2 >0
            div_y[i,j,k] = (scalar_field[id2] - scalar_field[id1])/(2*dy)
        end
        if id1 < 0
            div_y[i,j,k] = (scalar_field[id2] - scalar_field[id])/(dy)
        end

        if id2 < 0
            div_y[i,j,k] = (scalar_field[id] - scalar_field[id1])/(dy)
        end

    end

    return div_y
end

function divergence_z(scalar_field::Array{T,1}, mesh) where T<:AbstractFloat
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dz = mesh.dz
    div_z = zeros(T, nx, ny, nz)

    for i=1:nx,j=1:ny,k=1:nz
        id = index(i, j, k, nx, ny, nz)
        id1,id2 = index(i,j,k-1, nx,ny,nz),index(i,j,k+1, nx,ny,nz)
        if id1 > 0 && id2 >0
            div_z[i,j,k] = (scalar_field[id2] - scalar_field[id1])/(2*dz)
        end
        if id1 < 0
            div_z[i,j,k] = (scalar_field[id2] - scalar_field[id])/(dz)
        end

        if id2 < 0
            div_z[i,j,k] = (scalar_field[id] - scalar_field[id1])/(dz)
        end

    end

    return div_z
end

function second_divergence_z(scalar_field::Array{T,1}, mesh) where T<:AbstractFloat
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dz = mesh.dz
    div2_z = zeros(T, nx, ny, nz)

    for i=1:nx,j=1:ny,k=1:nz
        id = index(i, j, k, nx, ny, nz)
        id1,id2 = index(i,j,k-1, nx,ny,nz),index(i,j,k+1, nx,ny,nz)
        if id1 > 0
            div2_z[i,j,k] += (scalar_field[id1] - scalar_field[id])/dz^2
        end

        if id2 > 0
            div2_z[i,j,k] += (scalar_field[id2] - scalar_field[id])/dz^2
        end

    end

    return div2_z
end


function OVF2MFM(m,nx,ny,nz,dx,dy,dz,Nx,Ny,height)

    hz = floor(Int, height/dz)

    mesh = FDMesh(nx=nx,ny=ny,nz=nz+hz+1,dx=dx,dy=dy,dz=dz)

    println(@sprintf("mesh size is : %d %d %d \n",nx,ny,nz+hz+1))

    println("Creating sim...")

    sim = Sim(mesh,driver = "none")

    println("Initializing magnetization...")

    b = reshape(m,(3,nx,ny,nz))

    mm = zeros(3*nx*ny*(nz+hz+1))

    bb = reshape(mm,(3,nx,ny,nz+hz+1))

    Ms_array=zeros(nx*ny*(nz+hz+1))

    Ms = reshape(Ms_array,(nx,ny,nz+hz+1))

    # find the boundry of magnetized material

    for i = 1:nx,j=1:ny,k=1:nz

        if abs(b[1,i,j,k])+abs(b[2,i,j,k])+abs(b[3,i,j,k]) > 1e-1

           Ms[i,j,k] = 1e5

           bb[1,i,j,k] = b[1,i,j,k]

           bb[2,i,j,k] = b[2,i,j,k]

           bb[3,i,j,k] = b[3,i,j,k]

        end

    end

    set_Ms(sim,Ms_array)

    init_m0(sim,mm)


    println("Initializing demag kernel...")

    add_demag(sim)

    println("Computing demag field...")


    effective_field(sim,sim.spin,0.0)

    println("Computing divergence of demag field...")

    fieldz = reshape(sim.field,(3,:))[3,:]

    dFdz = second_divergence_z(fieldz, sim.mesh)

    println("All finished!")

    return dFdz[:,:,nz+hz]
end


"""
        OVF2MFM(fname::String;Nx=-1,Ny=-1,height=100e-9)

Get the simulated magnetic force microscope contrast from an ovf file.

And return the contrast matrix.

Method: Calculate the divergence of demagnetization on the platform of probe height.



                   probe |___________________________
                                                     |
                                                     |
                                                     |
                                                     |
             _____________________________________   | height
           /                                     /|  |
          /                                     / |  |
         /                                     / /   |
        /                ____________________________|
       /                                     / /
      /      material                       / /
     /                                     / /
    /_____________________________________/ /
    |_____________________________________|/

The tip is considered as a single domain along z-axis.

E = Integral(-dM.H) where dM is the local spin on the tip ( Mtip = Integral(dM) ), so dM ~ dz*(0,0,1).
                    H is the demagnetization field

So, phase ~ d2E/dz2 ~ dHz/dz, where Hz is the z-component of demagnetization field

We just need to calculate the z-ladder of Hz at the experimental height
"""
function OVF2MFM(fname::String;Nx=-1,Ny=-1,height=100e-9)
    np =  pyimport("numpy")
    mpl =  pyimport("matplotlib")
    mpl.use("Agg")
    mcolors = pyimport("matplotlib.colors")
    plt = pyimport("matplotlib.pyplot")
    colorsys = pyimport("colorsys")
    ag = pyimport("mpl_toolkits.axes_grid1")


    if endswith(fname,".ovf")
        fname= fname[1:end-4]
    end

    path=joinpath(dirname(fname),"MFM")

    if !isdir(path)
        mkpath(path)
    end

    ovf = read_ovf(fname)
    m = ovf.data
    nx = ovf.xnodes
    ny = ovf.ynodes
    nz = ovf.znodes
    dx = ovf.xstepsize
    dy = ovf.ystepsize
    dz = ovf.zstepsize

    println(@sprintf("data size is : %d %d %d, cell size is: %1.2ge-9 %1.2ge-9 %1.2ge-9 \n",nx,ny,nz,dx/1e-9,dy/1e-9,dz/1e-9))

    mfm = OVF2MFM(m,nx,ny,nz,dx,dy,dz,Nx,Ny,height)

    plt.imshow(np.transpose(mfm),cmap="gray",origin="lower")
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()
    figname = @sprintf("%s_MFM_h_%g.png",joinpath(path,basename(fname)),height/1e-9)
    plt.savefig(figname,dpi=300)
    plt.close()

    return mfm
end
