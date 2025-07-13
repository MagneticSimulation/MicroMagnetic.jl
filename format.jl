using JuliaFormatter

function find_newest_modified_files(folder_path::String, n::Int)
    # This will hold tuples of (modification_time, file_path)
    files_heap = []

    for (root, dirs, files) in walkdir(folder_path)
        if occursin(r"node_modules($|[/\\])", root)
            empty!(dirs)
            continue
        end

        for file in files
            file_path = joinpath(root, file)
            stats = stat(file_path)
            modification_time = round(Int64, stats.mtime)

            if lowercase(last(file_path, 3)) in [".jl", ".md"]
                push!(files_heap, (modification_time, file_path))
            end
        end
    end

    # Sort the heap by modification time in descending order
    sort!(files_heap; by=x -> x[1], rev=true)

    return [file[2] for file in files_heap[1:n]]
end

folder_path = "."
files = find_newest_modified_files(folder_path, 5)
for f in files
    format_file(f, YASStyle(); join_lines_based_on_source=false)
    println(f, " formated!")
end
