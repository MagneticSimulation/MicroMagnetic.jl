using PyCall
using Printf

np =  pyimport("numpy")
mpl =  pyimport("matplotlib")
mpl.use("Agg")
plt = pyimport("matplotlib.pyplot")


function plot_m() 
    data = np.loadtxt("std4_dyn.txt")
    oommf = np.loadtxt("oommf.txt")

    plt.plot(data[:,2],data[:,4],"-",label="JuMag_mx")
    plt.plot(data[:,2],data[:,5],"-",label="JuMag_my")
    plt.plot(data[:,2],data[:,6],"-",label="JuMag_mz")

    plt.plot(oommf[:,1],oommf[:,2],"-",dashes=(1,1),label="oommf_mx")
    plt.plot(oommf[:,1],oommf[:,3],"-",dashes=(1,1),label="oommf_my")
    plt.plot(oommf[:,1],oommf[:,4],"-",dashes=(1,1),label="oommf_mz")
    plt.ylabel("m")
    plt.xlabel("time")
    plt.legend()

    plt.savefig("std4.png")
    plt.close()
end

plot_m()