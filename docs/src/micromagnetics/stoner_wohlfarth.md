```@meta
ShareDefaultModule = true
```

# Stoner–Wohlfarth model

link: <https://en.wikipedia.org/wiki/Stoner%E2%80%93Wohlfarth_model>

When the size of the studied system is below the exchange length, the magnetization of the system is uniform
and thus can be described using a single 3d vector. In this situation, the demagnetization field can be calculated
by simply multiplying the demagnetization tensor and the magnetization. Therefore, the demagnetization energy is
equivalent to an effective anisotropy.
Assuming that the external field $\vec{H}=(H, 0, 0)$, easy axis $\hat{u}=(\cos\theta, \sin\theta, 0)$ and unit
magnetization vector $\hat{m}= (\cos\phi, \sin\phi, 0)$, we have
```math
\begin{align}
E &= - K (\vec{m} \cdot \hat{u})^2 - \mu_0 M_s H \cos \phi \\
&= -\frac{K}{2}\left[ 1 + \cos (2(\theta-\phi))+ 4 h \cos\phi \right ]
\end{align}
```
where $h=H/H_k$ and $H_k = 2K/(\mu_0 M_s)$. In the equilibrium state, the first derivative of the energy with respect to the
magnetization direction is zero, i.e.,
```math
\frac{\partial E}{\partial \phi} = K [-\sin (2(\theta-\phi)) + 2 h \sin \phi] = 0
```
In principle, there will be a solution of $\phi$ for each given $h$ and $\theta$. The hysteresis loops can be constructed
by plotting $\cos\phi$ as a function of  $h$. The switching field can be obtained by extra setting the second derivative
of the energy with respect to the magnetization direction to zero, i.e.,
```math
\frac{\partial^2 E}{\partial \phi^2} = 2 K [\cos (2(\theta-\phi)) +  h \cos \phi] = 0
```
The obtained switching field is
```math
h_s=\frac{\left(1-t^2+t^4\right)^{1 / 2}}{1+t^2}
```
where $t=\tan ^{1 / 3} \theta$. Specifically, $h_s=1/2$ if $\theta=\pi/4$. In this example, we use MicroMagnetic to demostrate
this result.

We chose the system to be a cubic sample, so the demagnetization tensor is $N_x=N_y=N_z=1/3$.
That is, the demagnetization itself does not contribute to the effective anisotropy. So in this simulation,
we have ignored the demagnetization field. The MicroMagnetic script is shown below:

````@example
using MicroMagnetic
````

We create a mesh for a cubic geometry 4nm x 4nm x 4nm

````@example
mesh = FDMesh(; nx=4, ny=4, nz=4, dx=1e-9, dy=1e-9, dz=1e-9);
nothing #hide
````
We create a simulation instance with the LLG driver and the BS23 integrator.
````@example
sim = Sim(mesh, driver="LLG", integrator="BS23", name="SW")

sim.driver.alpha = 0.5
sim.driver.integrator.tol = 1e-7

set_Ms(sim, 1e6)
init_m0(sim, (-1,1,0))
    
add_exch(sim, 1.3e-11)
add_anis(sim, 5e4, axis=(1, 1, 0))
add_zeeman(sim, (0,0,0))
nothing #hide
````

We use the `hysteresis` function to simulate the hysteresis loops.
````@example
Hs = [i*mT for i=-100:5:100]
hysteresis(sim, Hs, direction=(1,0,0), full_loop=false, stopping_dmdt=0.05, output="vts")
nothing #hide
````

For the used anisotropy $K_u=5\times 10^4$ A/m$^3$, the expected switch field is $H_c = (1/2) H_K = 39788.7$ A/m.
We plot the hysteresis loops using the `plot_ts` function

````@example
using CairoMakie
fig = plot_ts("SW_llg.txt", x_key="Hx", ["m_x", "m_y"], x_unit=1/mT, xlabel="H (mT)", ylabel="m", mirror_loop=true);
````

```@setup
save("../public/sw.png", fig)
```