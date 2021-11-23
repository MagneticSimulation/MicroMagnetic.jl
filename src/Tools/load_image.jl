using PyCall

function rgb_to_mag(fname; nx=-1, ny=-1, nz=1)
    rgb = load(fname)
    (Nx,Ny) = size(rgb)
    nx = nx > 0 ? nx : Nx
    ny = ny > 0 ? ny : Ny
    
    rgb = imresize(rgb,(nx,ny))
    
    mag_array = zeros(3,nx,ny,nz)
    
    for i = 1:nx, j=1:ny
        h,s,v = colorsys.hsv_to_rgb(rgb[i,j].r, rgb[i,j].g, rgb[i,j].b)
        sh, ch = sin(2*pi*h), cos(2*pi*h)
        mag_array[1,i,j,:] .= ch * sqrt(s)
        mag_array[2,i,j,:] .= sh/(ch^2) * sqrt(s)
        mag_array[3,i,j,:] .= 2*v -1
    end
    return mag_array
end