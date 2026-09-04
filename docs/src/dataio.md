# Data Input/Output

MicroMagnetic writes magnetization snapshots in the **OVF** format (the vector-field
format used by OOMMF and mumax3), meshes and fields in **VTK** (for ParaView), and
scalar time series in plain **text tables**. This page lists the available readers
and writers.

## OVF files

OVF is the recommended exchange format. It round-trips with OOMMF and mumax3
datasets and is what the [Standard Problem](micromagnetics/std4.md) tutorials use.

```julia
# write: snapshot of the current magnetization
save_ovf(sim, "m0.ovf")                     # data stored as Float64
save_ovf(sim, "m0.ovf"; type=Float32)       # smaller files

# read: either straight into a simulation, or into a container
read_ovf(sim, "m0.ovf")                     # loads the magnetization into sim.spin
ovf = read_ovf("m0.ovf")                    # returns an OVF2 struct
nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
m = ovf.data                                # flat 3N vector, interleaved [mx,my,mz,...]
```

`OVF2` objects can be written back with `save_ovf(ovf, "copy.ovf")`. To convert a
raw flat array (for example one produced by your own post-processing) into an OVF
file, use [`mag2ovf`](@ref):

```julia
mag2ovf(m, nx, ny, nz; dx=2.5e-9, dy=2.5e-9, dz=3e-9, fname="processed")
```

The `read_ovf` container is also what the analysis tools accept — see
[Analysis Tools](tools.md) for `compute_skyrmion_number("m.ovf")` and friends.

## VTK files

VTK is for ParaView-style visualization of meshes, magnetization and fields:

```julia
save_vtk(sim, "config")            # mesh + m + Ms (+ demag charges if enabled)
save_vtk(mesh, shape, "shape1")    # visualize a Shape on the mesh (see basics)
```

Existing OVF snapshots can be batch-converted for ParaView:

```julia
ovf2vtk("skx_LLG")                 # converts every OVF in the folder
```

## Time-series tables

Every run appends to a plain-text table `"<name>_<driver>.txt"` with time, average
magnetization and per-term energies (a sample is shown on the
[basics page](basics.md#data-tables)). Read it back as dictionaries:

```julia
data, units = read_table("std4_llg.txt")
data["time"]; data["m_x"]; data["E_total"]
```

Additional columns during a run are added with a custom `SaverItem` —
the recipe (guiding center example) is on the [basics page](basics.md#data-tables).

## Movies and plots

```julia
ovf2movie("std4_LLG"; output="std4.mp4", component='x')  # folder of OVF -> mp4
ovf2png("std4_LLG")                                      # folder of OVF -> pngs
plot_m(sim; component='z')                               # live view of sim
plot_ts("std4_llg.txt")                                  # quick table plot
```

## Which format when?

| Situation | Format |
| :--- | :--- |
| Exchanging snapshots with OOMMF/mumax3, or feeding analysis tools | OVF |
| 3D visualization in ParaView, meshes, vector fields | VTK |
| Time series you will read in Julia/Python | text table |
| Movies for papers/talks | mp4 via `ovf2movie` |
