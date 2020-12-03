using JuMC3DTri
using Printf
using Statistics 
using LinearAlgebra
using CUDA
CUDA.seed!(1000)

function rand_m(i, j, k, nx, ny, nz)
  return 2*rand(3).-1
end

function rand_m_uniform(i, j, k, nx, ny, nz)
  return (0,0,1)
end

# JuMag.cuda_using_double(true)

function TriAngle3D()
  mesh =  TriMesh3DGPU(nx=16, ny=16, nz=2,pbc="xyz")
  sim = MonteCarlo(mesh, name="mc")
  init_m0(sim, rand_m)

  #add_exch中有除以k_B
  set_exch(sim,Jxy=10*k_B,Jz=10*k_B)
  set_anis_6fold2d(sim,K1=40*k_B,K2=80*k_B)
  set_zeeman(sim,Hx=0.1*k_B,Hy=0,Hz=1.5*k_B)

  sim.T = 100000
  # run_sim(sim, maxsteps=20000, save_vtk_every=-1, save_m_every=-1)
  run_sim(sim, maxsteps=20, save_vtk_every=-1, save_m_every=-1)
  

  # NSample=5000
  NSample=5
	mx=zeros(NSample)
	my=zeros(NSample)
	mz=zeros(NSample)
	ms=zeros(NSample)
	Qs=zeros(NSample)
	EnMean =zeros(NSample)
	EnStd  =zeros(NSample)
	mxStd  =zeros(NSample)
	myStd  =zeros(NSample)
	mzStd  =zeros(NSample)
  f = open("table.txt", "w")
  write(f, "#T(K)     m      Q    EnPerSite(k_B)    Cv    ChiX ChiY ChiZ   mx   my   mz\n")

  # TArray=[ [i for i in 300.1:-50:100.1];[i for i in 90.1:-10:60.1];[i for i in 59.1:-1:0.1] ]
  TArray=[ 300,100 ]

  for T in TArray
      sim.T = T
      println("Running for $T ...")
      run_sim(sim, maxsteps=50, save_vtk_every=-1, save_m_every=-1)
      # run_sim(sim, maxsteps=50000, save_vtk_every=-1, save_m_every=-1)

      for i = 1:NSample
          run_sim(sim, maxsteps=50, save_vtk_every=-1, save_m_every=-1)
          calcEnAndStd(sim)
          mx[i]=sim.mxMean
          my[i]=sim.myMean
          mz[i]=sim.mzMean
          ms[i] = sqrt(sim.mxMean^2+sim.myMean^2+sim.mzMean^2)
          Qs[i]=sim.Q
          EnMean[i] = sim.EnMean
          EnStd[i]=sim.En2Mean - sim.EnMean^2

          mxStd[i]=sim.mx2Mean-sim.mxMean^2
          myStd[i]=sim.my2Mean-sim.myMean^2
          mzStd[i]=sim.mz2Mean-sim.mzMean^2

      end

      m = sum(ms)/NSample
      mxout= sum(mx)/NSample
      myout= sum(my)/NSample
      mzout= sum(mz)/NSample
      Q  = sum(Qs)/NSample
      En = sum(EnMean)/NSample
      Cv = sum(EnStd)/NSample/T/T
      ChiX=sum(mxStd)/NSample/T
      ChiY=sum(myStd)/NSample/T
      ChiZ=sum(mzStd)/NSample/T

      write(f, "$T    $m    $Q    $En   $Cv  $ChiX $ChiY $ChiZ  $mxout $myout $mzout \n")
  save_vtk(sim, @sprintf("T%dK.vts",T))
  save_ovf(sim, @sprintf("T%dK.ovf",T))

  end
  close(f)
  # save_vtk(sim, "T0.1.vts")

  return
end

outdir=replace(PROGRAM_FILE,".jl"=>".out")
# rm(outdir,force=true,recursive=true)
mkpath(outdir)
cd(outdir)
rm("table.txt",force=true)
for i in readdir()
  if length(i)>3 && i[length(i)-3:length(i)]==".vts"
  rm(i,force=true)
  end
end
TriAngle3D()
