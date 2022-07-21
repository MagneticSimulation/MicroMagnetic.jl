
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def plot_loop():
    data = np.loadtxt('sw.txt')
    Hx, Hy = data[:, 7], data[:, 8]
    H = Hx/np.sqrt(2) + Hy/np.sqrt(2)

    H = H*1e-3 #convert H to kA

    mx, my = data[:, 2], data[:, 3] 
    m = mx/np.sqrt(2) + my/np.sqrt(2)

    for i in range(len(H)-1):
        if m[i]>0 and m[i+1]<0:
            print("Hc=", (H[i]+H[i+1])/2.0)

    fig = plt.figure(facecolor='w', figsize=(4, 3))
    plt.plot(H, m, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='0K', markersize=1.2)
    plt.plot(-H, -m, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='', markersize=1.2)

    plt.grid(c='lightgrey')
    plt.ylabel("m")
    plt.xlabel("H (kA)")
    
    plt.tight_layout()
    fig.savefig('loop.png', dpi=300)

if __name__ == "__main__":
    plot_loop()