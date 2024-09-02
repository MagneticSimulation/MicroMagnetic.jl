function formatstring(s::String)
    return @sprintf("%18s ", s)
end

function formatstring(x::Int)
    return @sprintf("%18d ", x)
end

function formatstring(x::Float64)
    return @sprintf("%18.12e ", x)
end

function formatstring(x::Float32)
    return @sprintf("%18.12e ", x)
end

function formatstring(data::Tuple)
    # TODO: ("A","B","C") -> ("<A>", "<B>", "<C>")
    s = ""
    for x in data
        s = string(s, formatstring(x))
    end
    return s
end

"""
    init_saver(name::String, driver::String)

Use this function if a custom call back function is required. 

-`name`: output file name.
-`driver`: the driver name of your simulation. If is "LLG", the simulation time will be stored.
        
For example:

```julia
    saver = init_saver("output.txt", driver="LLG")
    skyrmion_center = SaverItem("center", "m", compute_guiding_center)
    run_sim(sim, steps=100, dt=1e-10, save_m_every=1, saver=saver)
```
"""
function init_saver(name::String, driver::String)
    saver = DataSaver(name, false, 0.0, 0, [])

    step = SaverItem("step", "<unitless>", o::AbstractSim -> o.saver.nsteps)
    push!(saver.items, step)

    if driver in ("LLG", "LLG_STT", "LLG_STT_CPP")
        time = SaverItem("time", "<s>", o::AbstractSim -> o.saver.t)
        push!(saver.items, time)
    end
    return saver
end

function create_saver(name::String, driver::String)
    saver = init_saver(name, driver)

    total_energy = SaverItem("E_total", "<J>", o::AbstractSim -> sum(o.energy))
    push!(saver.items, total_energy)

    m_all = SaverItem(("m_x", "m_y", "m_z"), ("<unitless>", "<unitless>", "<unitless>"), average_m)
    push!(saver.items, m_all)

    return saver
end

function write_data(sim::AbstractSim)
    return write_data(sim, sim.saver)
end

function save_sim_data(sim::AbstractSim)
    return write_data(sim, sim.saver)
end

function write_data(sim::AbstractSim, saver::DataSaver)
    if !saver.header_saved
        io = open(saver.name, "w")
        write(io, "#")
        for item in saver.items
            write(io, formatstring(item.name))
        end
        write(io, "\n#")

        for item in saver.items
            write(io, formatstring(item.unit))
        end
        write(io, "\n")
        saver.header_saved = true
    else
        io = open(saver.name, "a")
    end

    for item in saver.items
        write(io, formatstring(item.result(sim)))
    end

    write(io, "\n")
    return close(io)
end


"""
    read_table(filename::String)

Reads a data table from the specified file `filename`. The file should contain comment lines starting with `#`, where the first comment line 
includes the physical quantity names and the second comment line includes the units. The data starts from the `start_idx` line, which is the first line of actual data.

### Parameters
- `filename::String` : Path to the data file.

### Returns
- Two dictionaries:
  - `data` : A dictionary where keys are the physical quantity names and values are vectors of corresponding floating-point data.
  - `units` : A dictionary where keys are the physical quantity names and values are the corresponding unit strings.

### Example
```julia
data, units = read_table("data.txt")
"""
function read_table(filename::String)
 
    data_lines = open(filename) do file
        readlines(file)
    end

    start_idx = findfirst(line -> !startswith(line, "#"), data_lines)
    
    headers = split(strip(replace(data_lines[start_idx - 2], r"[#]" => "")), r"\s+")
    units = split(strip(replace(data_lines[start_idx - 1], r"[#<>]" => "")), r"\s+")

    for i in 1:length(units)
        if units[i] == "1" || isempty(units[i])
            units[i] = "unitless"
        end
    end

    data_dict = Dict{String, Vector{Float64}}()
    units_dict = Dict{String, String}()

    for line in data_lines[start_idx:end]
        values = split(strip(line), r"\s+")
        for (i, header) in enumerate(headers)
            if !haskey(data_dict, header)
                data_dict[header] = Float64[]
            end
            push!(data_dict[header], parse(Float64, values[i]))
        end
    end

    for (i, header) in enumerate(headers)
        units_dict[header] = units[i]
    end

    return (data_dict, units_dict)
end

export read_table