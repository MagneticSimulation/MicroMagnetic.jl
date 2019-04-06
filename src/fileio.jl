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
  s = ""
  for x in data
    s = string(s, formatstring(x))
  end
  return s
end

function write_data(sim::AbstractSim)
  saver = sim.saver
  if !saver.header_saved
    io = open(saver.name, "w");
    write(io, "#")
    for header in saver.headers
      write(io, formatstring(header));
    end
    write(io, "\n#");

    for s in saver.units
      write(io, formatstring(s));
    end
    write(io, "\n")
    saver.header_saved = true
  else
    io = open(saver.name, "a");
  end

  for fun in saver.results
    write(io, formatstring(fun(sim)));
  end
  write(io, "\n")
  close(io);
end
