
using Printf
using JLD2
using PNGFiles
using Colors, ColorSchemes
using FFMPEG

"""
    ovf2png(ovf_name; k = 1, component='z')

Create a png from the given ovf file.
`k` indicates the layer index (starting from 1) 
`component` refers to the used component ('x', 'y' or 'z') for plotting. 

"""
function ovf2png(ovf_name, output=nothing; k = 1, component='z')
  if output == nothing
    output = endswith(ovf_name, ".ovf") ? ovf_name[1:end-4] : ovf_name
  end

  ovf = read_ovf(ovf_name)

  spin = reshape(ovf.data, 3, ovf.xnodes, ovf.ynodes, ovf.znodes)

  m_index = component - 'x' + 1
  m_t = transpose(spin[m_index, :, :, k]) # transpose the magnetization
  m_f = reverse(m_t, dims = 1) # flip the m_t vertically
  PNGFiles.save(output*".png", get(ColorSchemes.coolwarm, m_f))
end



"""
    jdl2png(jdl_file; k = 1, component='z')

Create pngs from the given jdl2 file where `k` indicates the layer index (starting from 1) 
and `component` refers to the used component ('x', 'y' or 'z') for plotting. 

"""
function jdl2png(jdl_file; k = 1, component='z')
    data = load(jdl_file)
    steps = data["steps"]
    save_m_every = data["save_m_every"]
    nx, ny, nz = data["mesh/nx"], data["mesh/ny"], data["mesh/nz"]
    (k < 1) && (k = 1)
    (k > nz) && (k = nz)
    m_index = component - 'x' + 1
    if save_m_every < 0
      @info @sprintf("save_m_every is %d, which is negative, exiting~", save_m_every)
      return
    end

    png_folder = jdl_file[1:length(jdl_file)-5]*"_pngs"
    isdir(png_folder) || mkdir(png_folder)
    for i=1:save_m_every:steps
      index = @sprintf("m/%d", i)
      m = reshape(data[index], 3, nx, ny, nz)
      m_t = transpose(m[m_index, :, :, k]) # transpose the magnetization
      m_f = reverse(m_t, dims = 1) # flip the m_t vertically
      png_name = @sprintf("%s/img_%d.png", png_folder, i) 
      PNGFiles.save(png_name, get(ColorSchemes.coolwarm, m_f))
    end

    return png_folder
end

"""
  jdl2movie(jdl_file; k = 1, component='z', r=12, rm_png=false, output=nothing)

Create a moive from the given jdl2 file. `k` indicates the layer index (starting from 1), 
`component` refers to the used component ('x', 'y' or 'z') for plotting, `r` is frame rate, 
output is the filename of the video and the support formats are 'mp4', 'avi' and 'gif'.
"""
function jdl2movie(jdl_file; k = 1, component='z', r=12, rm_png=false, output=nothing)
  if output==nothing
    base_name = jdl_file[1:length(jdl_file)-5]
    output = @sprintf("%s.mp4", base_name)
  end

  isfile(output) && rm(output, force=true)

  png_folder = jdl2png(jdl_file, k=k, component=component)

  ffmpeg_exe( "-loglevel", "fatal", "-nostats", "-r", @sprintf("%d", r), "-i", png_folder*"/img_%d.png", output)

  rm_png && rm(png_folder, force=true, recursive=true)
end





