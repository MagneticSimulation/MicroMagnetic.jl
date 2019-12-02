using JuMag
using PyCall

kx, ky, intensity = JuMag.fft_m("skx.ovf", axis='z')

np =  pyimport("numpy")
mpl =  pyimport("matplotlib")
mpl.use("Agg")
plt = pyimport("matplotlib.pyplot")

kx = kx*1e-9
ky = ky*1e-9

cs = plt.imshow(np.transpose(intensity), interpolation = "none", origin="lower",
                extent=(kx[1], kx[end], ky[1], ky[end]), cmap=mpl.cm.PuBu_r)

plt.xlim([-0.2, 0.2])
plt.ylim([-0.2, 0.2])
plt.xlabel("kx")
plt.ylabel("ky")
cbar = plt.colorbar(cs)
plt.savefig("fftz.png")
