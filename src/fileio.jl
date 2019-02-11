function formatstring(s::String)
  return @sprintf("%18s", s)
end

function formatstring(x::Number)
  return @sprintf("%18.12g ", x)
end

function formatstring(data::Tuple)
  s = ""
  for x in data
    s = string(s, formatstring(x))
  end
  return s
end

function write_data(sim::SimData)
  saver = sim.saver
  io = open(saver.name, "a");
  if !saver.header_saved
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
  end

  for fun in saver.results
    write(io, formatstring(fun(sim)));
  end
  write(io, "\n")
  close(io);
end
