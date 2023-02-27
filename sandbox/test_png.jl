using PNGFiles
using Colors, ColorSchemes

M, N = 500, 200
data = zeros(M, N)

for i = 1:M, j=1:N
    x = cos(i/M*2*pi)*sin(j/N*2*pi)*0.5+0.5
    data[i, j] = x
end

PNGFiles.save("img.png", get(ColorSchemes.coolwarm, data))

PNGFiles.save("img2.png", rand(RGB, M, N))
