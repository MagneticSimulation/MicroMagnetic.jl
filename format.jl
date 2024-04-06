using JuliaFormatter

function find_newest_modified_file(folder_path::String)
    newest_file = ""
    newest_time = 0

    for (root, dirs, files) in walkdir(folder_path)
        for file in files
            file_path = joinpath(root, file)

            stats = stat(file_path)
            modification_time = round(Int64, stats.mtime)

            if modification_time > newest_time
                newest_time = modification_time
                newest_file = file_path
            end
        end
    end

    return newest_file
end

folder_path = "."
newest_modified_file = find_newest_modified_file(folder_path)
println("Newest modified file: ", newest_modified_file)

if lowercase(last(newest_modified_file, 3)) in [".jl", ".md"]
    format_file(newest_modified_file, YASStyle(); join_lines_based_on_source=false)
    println("Formated!")
end
