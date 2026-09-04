# Bitwise regression for the arithmetic-addressing stencil kernels:
# new kernels (ngbs-free) must reproduce the master ngbs-table kernels
# (snapshot in stencil_ref_kernels.jl) exactly — including periodic boundaries,
# vacuum regions, spatial (inline) parameters and odd mesh dimensions.
# Context: NGBS_ARITH_TASK.md; design+experiments in EXCH_DMI_OPT.md §7.

using MicroMagnetic
using MicroMagnetic: Fill
using Test
using KernelAbstractions
using Printf

module RefKernels
using KernelAbstractions
import MicroMagnetic
const cross_x = MicroMagnetic.cross_x
const cross_y = MicroMagnetic.cross_y
const cross_z = MicroMagnetic.cross_z
include(joinpath(@__DIR__, "stencil_ref_kernels.jl"))
end

const _ref_ngbs(mesh) = MicroMagnetic.kernel_array(Array(mesh.ngbs))
const _ref_mat2(v::T) where {T} = MicroMagnetic.kernel_array(fill(v, 2, 2))

function ref_exch!(h2, e2, sim, exch, cls, ok)
    T = eltype(sim.spin); mesh = sim.mesh; back = get_backend(sim.spin)
    ngbs = _ref_ngbs(mesh)
    if ok
        pAx = exch.pair_Ax isa Fill ? _ref_mat2(exch.pair_Ax.value) : exch.pair_Ax
        pAy = exch.pair_Ay isa Fill ? _ref_mat2(exch.pair_Ay.value) : exch.pair_Ay
        pAz = exch.pair_Az isa Fill ? _ref_mat2(exch.pair_Az.value) : exch.pair_Az
        RefKernels.exchange_partition_kernel!(back, 256)(sim.spin, h2, e2, cls, pAx, pAy,
            pAz, sim.inv_ms, T(mesh.dx), T(mesh.dy), T(mesh.dz), ngbs, T(mesh.volume);
            ndrange=sim.n_total)
    else
        RefKernels.exchange_kernel!(back, 256)(sim.spin, h2, e2, sim.mu0_Ms, exch.Ax,
            exch.Ay, exch.Az, T(mesh.dx), T(mesh.dy), T(mesh.dz), ngbs, T(mesh.volume);
            ndrange=sim.n_total)
    end
end

function ref_dmi!(h2, e2, sim, dmi, cls, ok)
    T = eltype(sim.spin); mesh = sim.mesh; back = get_backend(sim.spin)
    ngbs = _ref_ngbs(mesh)
    tfac = T(dmi.ft(0.0))
    if dmi.type === :bulk
        if ok
            pDx = dmi.pair_Dx isa Fill ? _ref_mat2(dmi.pair_Dx.value) : dmi.pair_Dx
            pDy = dmi.pair_Dy isa Fill ? _ref_mat2(dmi.pair_Dy.value) : dmi.pair_Dy
            pDz = dmi.pair_Dz isa Fill ? _ref_mat2(dmi.pair_Dz.value) : dmi.pair_Dz
            RefKernels.bulkdmi_partition_kernel!(back, 256)(sim.spin, h2, e2, cls, pDx,
                pDy, pDz, sim.inv_ms, T(mesh.dx), T(mesh.dy), T(mesh.dz), ngbs,
                T(mesh.volume); ndrange=sim.n_total)
        else
            RefKernels.bulkdmi_kernel!(back, 256)(sim.spin, h2, e2, sim.mu0_Ms, dmi.Dx,
                dmi.Dy, dmi.Dz, T(mesh.dx), T(mesh.dy), T(mesh.dz), ngbs,
                T(mesh.volume), tfac; ndrange=sim.n_total)
        end
    else
        if ok
            pD = dmi.pair_Dx isa Fill ? _ref_mat2(dmi.pair_Dx.value) : dmi.pair_Dx
            RefKernels.interfacial_dmi_partition_kernel!(back, 256)(sim.spin, h2, e2, cls,
                pD, dmi.Dcls, sim.inv_ms, T(mesh.dx), T(mesh.dy), ngbs, T(mesh.volume);
                ndrange=sim.n_total)
        else
            RefKernels.interfacial_dmi_kernel!(back, 256)(sim.spin, h2, e2, sim.mu0_Ms,
                dmi.Dx, T(mesh.dx), T(mesh.dy), T(mesh.dz), ngbs, T(mesh.volume), tfac;
                ndrange=sim.n_total)
        end
    end
end

function ref_zhangli!(h2, e2, sim, stt, cls, ok)
    T = eltype(sim.spin); mesh = sim.mesh; back = get_backend(sim.spin)
    ut = T(stt.ufun(0.0) / sim.driver.gamma)
    RefKernels.zhangli_torque_kernel!(back, 256)(sim.spin, h2, stt.bJx, stt.bJy, stt.bJz,
        _ref_ngbs(mesh), stt.xi, ut, T(mesh.dx), T(mesh.dy), T(mesh.dz);
        ndrange=sim.n_total)
end

function build_case(nx, ny, nz, pbc, pattern)
    mesh = FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=nx, ny=ny, nz=nz, pbc=pbc)
    sim = Sim(mesh)
    exch = nothing; dmi_b = nothing; dmi_i = nothing; stt = nothing
    if pattern == :uniform
        set_Ms(sim, 8.6e5)
        exch = add_exch(sim, 1.3e-11; name="exch_u")
        dmi_b = add_dmi(sim, 3.0e-3; type="bulk", name="dmiB_u")
        dmi_i = add_dmi(sim, 3.0e-3; type="interfacial", name="dmiI_u")
        stt = add_stt(sim; model=:zhang_li, J=(1e12, 0.0, 0.0), P=0.95, Ms=8.6e5, xi=0.05, name="stt_u")
    elseif pattern == :piecewise
        set_region(mesh, 1; i=1:div(nx, 2))
        set_region(mesh, 2; i=div(nx, 2)+1:nx)
        set_Ms(sim, 8.6e5)
        exch = add_exch(sim, region_map(1 => 1.3e-11, 2 => 2.6e-11; default=1.3e-11); name="exch_p")
        dmi_b = add_dmi(sim, region_map(1 => 3.0e-3, 2 => 1.5e-3; default=3.0e-3); type="bulk", name="dmiB_p")
        dmi_i = add_dmi(sim, region_map(1 => 3.0e-3, 2 => 1.5e-3; default=3.0e-3); type="interfacial", name="dmiI_p")
    elseif pattern == :vacuum
        set_region(mesh, 1; i=1:div(nx, 2))
        set_Ms(sim, region_map(1 => 8.6e5, 2 => 0.0; default=0.0))
        exch = add_exch(sim, 1.3e-11; name="exch_v")
        dmi_b = add_dmi(sim, 3.0e-3; type="bulk", name="dmiB_v")
        dmi_i = add_dmi(sim, 3.0e-3; type="interfacial", name="dmiI_v")
    elseif pattern == :spatial
        set_Ms(sim, 8.6e5)
        exch = add_exch(sim, (x, y, z) -> 1.3e-11 * (1 + 0.5 * sin(x * 2e9) * cos(y * 2e9)); name="exch_s")
        dmi_b = add_dmi(sim, (x, y, z) -> 3.0e-3 * (1 + 0.3 * sin(x * 1e9)); type="bulk", name="dmiB_s")
        dmi_i = add_dmi(sim, (x, y, z) -> 3.0e-3 * (1 + 0.3 * cos(y * 1e9)); type="interfacial", name="dmiI_s")
    end
    init_m0_random(sim; seed=42)
    return (; mesh, sim, exch, dmi_b, dmi_i, stt)
end

function check_terms!(fails, tag, c)
    sim = c.sim
    for (name, term, ref!, tables!) in (
        ("exch", c.exch, ref_exch!, MicroMagnetic._exchange_tables!),
        ("dmiB", c.dmi_b, ref_dmi!, MicroMagnetic._dmi_tables!),
        ("dmiI", c.dmi_i, ref_dmi!, MicroMagnetic._dmi_tables!),
        ("zhangli", c.stt, ref_zhangli!, nothing),
    )
        term === nothing && continue
        cls = nothing; ok = false
        if tables! !== nothing
            cls, ok = tables!(sim, term)
        end
        MicroMagnetic.effective_field(term, sim, sim.spin, 0.0)
        h2 = zero(term.field)
        has_energy = hasproperty(term, :energy)
        e2 = has_energy ? zero(term.energy) : nothing
        ref!(h2, e2, sim, term, cls, ok)
        same_h = isequal(Array(term.field), Array(h2))
        same_e = !has_energy || isequal(Array(term.energy), Array(e2))
        if !(same_h && same_e)
            push!(fails, @sprintf("%s %s (%s)", tag, name, ok ? "partition" : "inline"))
        end
    end
end

const PBC_MESHES = [
    (17, 13, 5, "open"), (16, 12, 8, "xyz"), (16, 12, 8, "xz"), (16, 12, 8, "y"),
    (24, 18, 1, "x"), (24, 18, 1, "xyz"), (33, 21, 9, "open"), (3, 2, 2, "open"),
]

@testset "stencil arith vs ngbs reference (bitwise, PBC)" begin
    fails = String[]
    for pattern in (:uniform, :piecewise, :vacuum, :spatial), (nx, ny, nz, pbc) in PBC_MESHES
        c = build_case(nx, ny, nz, pbc, pattern)
        tag = @sprintf("%s %dx%dx%d pbc=%s", pattern, nx, ny, nz, pbc)
        check_terms!(fails, tag, c)
    end
    @test isempty(fails)
    isempty(fails) || println(join(fails, "\n"))
end
