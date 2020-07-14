using LinearAlgebra
using CUDA

function compute_distance(neb::NEB_GPU_MPI)
    nxyz = neb.sim.nxyz
    N = neb.N
    dof = 3*nxyz
    blk, thr = cudims(nxyz)
    ds = neb.sim.energy #we borrow the sim.energy
    for n=0:N
       m1 = n==0 ? neb.image_l : view(neb.spin, (n-1)*dof+1:n*dof)
       m2 = n==N ? neb.image_r : view(neb.spin, n*dof+1:(n+1)*dof)
       @cuda blocks=blk threads=thr compute_distance_kernel!(ds, m1, m2, nxyz)
       neb.distance[n+1] = LinearAlgebra.norm(ds)
   end
   return nothing
end

function compute_tangents(neb::NEB_GPU_MPI)
    N = neb.N
    nxyz = neb.sim.nxyz
    blk, thr = cudims(3*nxyz)
    @cuda blocks=blk threads=thr compute_tangents_kernel!(neb.tangent, neb.spin,
                                                           neb.image_l, neb.image_r,
                                                           neb.energy, N, nxyz)

    blk, thr = cudims(N*nxyz)
    @cuda blocks=blk threads=thr neb_projection_kernel!(neb.tangent, neb.spin, N*nxyz)

    dof = 3*nxyz
    for n=1:N
        t = view(neb.tangent, dof*(n-1)+1:dof*n)
        norm_t = LinearAlgebra.norm(t)
        t .= t/norm_t
    end

   return nothing
end
