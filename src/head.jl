abstract type Driver end
abstract type AbstractSim end
abstract type MicroEnergy end
abstract type Integrator end
abstract type IntegratorCayley <: Integrator end

"""
    NumberOrArrayOrFunction

In Micromagnetics, typical parameters such as saturation magnetization `Ms` and exchange stiffness constant `A` are constant.
However, there are many cases that a spatial `Ms` is needed. For instance, if the simulated system is a circular disk the Ms
in the corners should be set to zero. If the simulated system contains mutiple materials, the exchange constant `A` should be
spatial as well. The Union `NumberOrArrayOrFunction` is designed to deal with such situations. As indicated from its name, it
means that the input parameter could be a number or an array or a function:
  - Number: should be Real.
  - Array: the length of the array should be `N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` where `i,j,k` is the cell index and `dx,dy,dz` is the cellsize.
    The return value of the function should be a real number. For example,

    ```julia
    function circular_Ms(i,j,k,dx,dy,dz)
        if (i-50.5)^2 + (j-50.5)^2 <= 50^2
            return 8.6e5
        end
        return 0.0
    end
    ```
"""
NumberOrArrayOrFunction = Union{Number,Array,Function}

# the aliases below appear in public method signatures, so export them (§8.5)
export NumberOrArrayOrFunction, NumberOrTupleOrArrayOrFunction, TupleOrArrayOrFunction

"""
    NumberOrTupleOrArrayOrFunction

In micromagnetics, there are also cases where the input parameters can be either scalars or vectors and vary with space. For example,
the parameters for the DMI could be a const for bulk DMI or interfacial DMI. In some materials, the DMI const may differ in different 
directions and thus a tuple with three numbers is required. In MicroMagnetic, the union `NumberOrTupleOrArrayOrFunction` is designed to deal 
with such situations. Similar to `NumberOrArrayOrFunction`, `NumberOrTupleOrArrayOrFunction` means that the input parameter could be 
a number, a tuple, an array or a function:
  - Number: should be Real.
  - Tuple: should be Real with length 3. For example, `(1,2e-5,0)`.
  - Array: the length of the array should be `N` or `3N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` and the return value should be a tuple with length 3.

    For example,
    ```julia
    function uniform_m0(i,j,k,dx,dy,dz)
        return (0,0,1)
    end
    ```
"""
NumberOrTupleOrArrayOrFunction = Union{Number,Tuple,Array,Function}

"""
    NumberOrArray

Similar to Union `NumberOrArrayOrFunction`, the Union `NumberOrArray` is designed to deal with cases that a number
or an array is needed.
"""
NumberOrArray = Union{Number,Array}

"""
    ArrayOrFunction

Similar to Union `NumberOrArrayOrFunction`, the Union `ArrayOrFunction` is designed to deal with cases that
an array or a function is needed.
"""
ArrayOrFunction = Union{Array,Function}

"""
    TupleOrArrayOrFunction

Similar to `NumberOrArrayOrFunction`, `TupleOrArrayOrFunction` means that the input parameter could be a tuple or
an array or a function:
  - Tuple: should be Real with length 3. For example, `(0,0,1e5)`.
  - Array: the length of the array should be `3N` where `N` is the total spin number of the system.
  - Function: the parameter of the function should be `(i,j,k,dx,dy,dz)` and the return value should be a tuple with length 3.
    For example,
    ```julia
    function uniform_m0(i,j,k,dx,dy,dz)
        return (0,0,1)
    end
    ```
"""
TupleOrArrayOrFunction = Union{Tuple,Array,Function}

"""
    TupleOrArray

The Union `TupleOrArray` is designed to deal with cases that a Tuple or an array is needed.
"""
TupleOrArray = Union{Tuple,Array}

mutable struct DataSaver
    name::String
    header_saved::Bool
    t::Float64
    nsteps::Int
    items::Array
end

mutable struct SaverItem
    name::Union{String,Tuple}
    unit::Union{String,Tuple}
    result::Function
end

"""
    MicroSim{T<:AbstractFloat} <: AbstractSim
"""
mutable struct MicroSim{T<:AbstractFloat} <: AbstractSim
    time::Float64
    mesh::FDMesh
    driver::Driver
    saver::DataSaver
    spin::AbstractArray{T,1}
    prespin::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    mu0_Ms::AbstractArray{T,1}
    pins::AbstractArray{Bool,1}
    n_total::Int64
    name::String
    driver_name::String
    interactions::Array{MicroEnergy}
    save_data::Bool
    # stencil-layer caches (see micro/stencil.jl): mat_class freshness keys on
    # mesh.layout_version (set_region bumps it; set_Ms drops caches directly);
    # inv_ms is derived at set_Ms. Both are caches — the parameter arrays stay
    # the permanent source.
    mat_class::Union{Nothing,AbstractArray}  # class per cell; 0 = vacuum; may be Fill
    n_classes::Int                           # material classes are 1..R
    mat_class_layout::Int
    inv_ms::Union{Nothing,AbstractArray}     # 1/mu0_Ms in T; 0 for vacuum; may be Fill
    MicroSim{T}() where {T<:AbstractFloat} = new()
end

mutable struct MicroSimFE{T<:AbstractFloat} <: AbstractSim
    mesh::FEMesh
    driver::Driver
    saver::DataSaver
    spin::AbstractArray{T,1}
    prespin::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    mu0_Ms::Array{T,1}
    L_mu::AbstractArray{T,1}
    pins::AbstractArray{Bool,1}
    n_total::Int64  # n_nodes
    n_cells::Int64
    name::String
    driver_name::String
    interactions::Array{MicroEnergy}
    save_data::Bool
    kwargs::Any
    MicroSimFE{T}() where {T<:AbstractFloat} = new()
end

mutable struct Zeeman{T<:AbstractFloat} <: MicroEnergy
    Hx::AbstractArray{T,1}  # per-spin static field components; uniform = O(1) Fill
    Hy::AbstractArray{T,1}
    Hz::AbstractArray{T,1}
    ft::Function            # time modulation; the static term uses `_static_time`
    H_repr::Any             # original input, kept for saver display / H_output
    time_fx::T              # current time factors (saver display)
    time_fy::T
    time_fz::T
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

# Magnetoelastic energy term (unified interface)
mutable struct Magnetoelastic{T<:AbstractFloat} <: MicroEnergy
    model::Symbol                 # :tensor, :cubic
    lambda_s::T                 # saturation magnetostriction (dimensionless, for :tensor)
    B1::T                      # magnetoelastic coupling B1 (J/m^3, for :cubic)
    B2::T                      # magnetoelastic coupling B2 (J/m^3, for :cubic)
    field_data::Any
    field::AbstractArray{T,1}   # effective field
    energy::AbstractArray{T,1}     # energy density
    name::String
end

mutable struct Anisotropy{T<:AbstractFloat} <: MicroEnergy
    Ku::AbstractArray{T,1}
    axis_x::AbstractArray{T,1}  # per-spin easy-axis components; uniform axes are Fills
    axis_y::AbstractArray{T,1}
    axis_z::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct CubicAnisotropy{T<:AbstractFloat} <: MicroEnergy
    Kc::AbstractArray{T,1}
    axis1x::AbstractArray{T,1}  # per-spin easy-axis components; uniform = O(1) Fill
    axis1y::AbstractArray{T,1}
    axis1z::AbstractArray{T,1}
    axis2x::AbstractArray{T,1}
    axis2y::AbstractArray{T,1}
    axis2z::AbstractArray{T,1}
    axis3x::AbstractArray{T,1}
    axis3y::AbstractArray{T,1}
    axis3z::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct HexagonalAnisotropy{T<:AbstractFloat} <: MicroEnergy
    K1::T
    K2::T
    K3::T
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct TwinMonoclinicAnisotropy{T<:AbstractFloat} <: MicroEnergy
    Ka::T
    Kb::T
    Kaa::T
    Kbb::T
    Kab::T
    Ku::T
    axis_a::AbstractArray{T,1}
    axis_b::AbstractArray{T,1}
    axis_u111::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct Exchange{T<:AbstractFloat} <: MicroEnergy
    Ax::AbstractArray{T,1}  # bond stiffness along x; uniform values are O(1) Fills
    Ay::AbstractArray{T,1}  # bond stiffness along y
    Az::AbstractArray{T,1}  # bond stiffness along z
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
    # stencil-layer caches (see micro/stencil.jl): pair_A* are harmonic-mean
    # pair tables built when Ax/Ay/Az are per-class uniform; stencil_layout is
    # the mesh.layout_version stamp of the last build (-1 = needs rebuild).
    pair_Ax::Union{Nothing,AbstractArray{T,2}}
    pair_Ay::Union{Nothing,AbstractArray{T,2}}
    pair_Az::Union{Nothing,AbstractArray{T,2}}
    stencil_layout::Int
end

mutable struct InterlayerExchange{T<:AbstractFloat} <: MicroEnergy
    Js::AbstractArray{T,1}
    k1::Int32
    k2::Int32
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

mutable struct InterlayerDMI{T<:AbstractFloat} <: MicroEnergy
    Dx::T
    Dy::T
    Dz::T
    k1::Int32
    k2::Int32
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

# The micromagnetic (continuum) DMI energy, serving the bulk and interfacial
# forms through `type`; atomistic DMI lives in the Heisenberg* interaction types.
mutable struct DMI{T<:AbstractFloat} <: MicroEnergy
    Dx::AbstractArray{T,1}  # per-spin D components; uniform values are Fills
    Dy::AbstractArray{T,1}
    Dz::AbstractArray{T,1}
    ft::Function            # time modulation; the static term uses `_static_time`
    type::Symbol            # :bulk or :interfacial
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
    # stencil-layer caches (see micro/stencil.jl): pair_D* are pair tables
    # built when Dx/Dy/Dz are per-class uniform; Dcls is the per-class D used
    # for the interfacial D_I==0 neighbour guard. stencil_layout is the
    # mesh.layout_version stamp of the last build (-1 = needs rebuild).
    pair_Dx::Union{Nothing,AbstractArray{T,2}}
    pair_Dy::Union{Nothing,AbstractArray{T,2}}
    pair_Dz::Union{Nothing,AbstractArray{T,2}}
    Dcls::AbstractArray{T,1}
    stencil_layout::Int
end

mutable struct StochasticField{T<:AbstractFloat} <: MicroEnergy
    temperature::AbstractArray{T,1}
    offset_temp::T
    eta::AbstractArray{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    nsteps::Int64
    name::String
    k_B::Float64
    scaling_fun::Function
    mean_temperature::T
    scaling_factor::T
    spatiotemporal_mode::Bool
end

mutable struct TorqueField{T<:AbstractFloat} <: MicroEnergy
    torque_fun::Function
    field::AbstractArray{T,1}
    name::String
end

mutable struct DFTorqueField{T<:AbstractFloat} <: MicroEnergy
    px::AbstractArray{T,1}  # per-spin polarization components; uniform = O(1) Fill
    py::AbstractArray{T,1}
    pz::AbstractArray{T,1}
    aj::AbstractArray{T,1}
    bj::T
    field::AbstractArray{T,1}
    name::String
end

mutable struct ZhangLiTorque{T<:AbstractFloat} <: MicroEnergy
    xi::AbstractArray{T,1}
    bJx::AbstractArray{T,1} # b*J along x; uniform values are O(1) Fills
    bJy::AbstractArray{T,1} # b*J along y
    bJz::AbstractArray{T,1} # b*J along z
    field::AbstractArray{T,1}
    ufun::Function
    name::String
end

mutable struct SlonczewskiTorque{T<:AbstractFloat} <: MicroEnergy
    px::AbstractArray{T,1}  # per-spin polarization components; uniform = O(1) Fill
    py::AbstractArray{T,1}
    pz::AbstractArray{T,1}
    beta::T
    Lambda::T
    xi::T
    P::T
    J::AbstractArray{T,1} # bJ = b*J
    field::AbstractArray{T,1}
    ufun::Function
    name::String
end

mutable struct SAHETorqueField{T<:AbstractFloat} <: MicroEnergy
    sigma_s_a1::AbstractArray{T,1} # = sigma_s x a1
    sigma_sa1_a1::AbstractArray{T,1} # = sigma_sa1 x a1 
    sigma_sa2::AbstractArray{T,1}
    a2::AbstractArray{T,1}  
    a3::AbstractArray{T,1}
    beta::T
    field::AbstractArray{T,1}
    name::String
end

mutable struct ExchangeFE{T<:AbstractFloat} <: MicroEnergy
    A::Array{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    K_matrix::AbstractSparseMatrix
    Mass_matrix::AbstractSparseMatrix
    method::Any
    name::String
end

mutable struct AnisotropyFE{T<:AbstractFloat} <: MicroEnergy
    Ku::Array{T,1}
    axis::Array{T,1}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    K_matrix::AbstractSparseMatrix
    name::String
end

mutable struct InterlayerExchangeFE{T<:AbstractFloat} <: MicroEnergy
    J::T
    s1::Int32 # surface id
    s2::Int32 # surface id
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    K_matrix::AbstractSparseMatrix
    name::String
end
