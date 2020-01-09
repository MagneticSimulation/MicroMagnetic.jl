using PyCall
using Printf

np =  pyimport("numpy")
mpl =  pyimport("matplotlib")
mpl.use("Agg")
plt = pyimport("matplotlib.pyplot")


function plot_m() 
    data_cpu = np.loadtxt("std5_dyn_cpu.txt")
    data_gpu = np.loadtxt("std5_dyn_gpu.txt")
    fidimag = np.loadtxt("fidimag.txt")

    plt.plot(data_cpu[:,2],data_cpu[:,4],"-",label="JuMag_mx")
    plt.plot(data_cpu[:,2],data_cpu[:,5],"-",label="JuMag_my")
    plt.plot(data_cpu[:,2],data_cpu[:,6],"-",label="JuMag_mz")

    plt.plot(data_gpu[:,2],data_gpu[:,4],"-",dashes=(2,2),label="JuMag_mx_gpu")
    plt.plot(data_gpu[:,2],data_gpu[:,5],"-",dashes=(2,2),label="JuMag_my_gpu")
    plt.plot(data_gpu[:,2],data_gpu[:,6],"-",dashes=(2,2),label="JuMag_mz_gpu")

    plt.plot(fidimag[:,1],fidimag[:,5],"-",dashes=(1,1),label="fidimag_mx")
    plt.plot(fidimag[:,1],fidimag[:,6],"-",dashes=(1,1),label="fidimag_my")
    plt.plot(fidimag[:,1],fidimag[:,7],"-",dashes=(1,1),label="fidimag_mz")
    plt.ylabel("m")
    plt.xlabel("time")
    plt.legend()

    plt.savefig("std5.png")
    plt.close()
end

function plot_centre()
    data = np.loadtxt("std5_dyn_cpu_centre.txt")
    plt.plot(data[:,1],data[:,2])

    plt.xlim = [0,1e-7]
    plt.ylim = [0,1e-7]
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("std5_centre.png")
    plt.close()
end

plot_m()
plot_centre()