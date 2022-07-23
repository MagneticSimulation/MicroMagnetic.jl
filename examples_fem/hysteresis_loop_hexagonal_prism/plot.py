
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12,
        }

def plot_loop():
    #kOe = 79577.47154594767
    data = np.loadtxt('hexagonal.txt')
    Hx, Hy = data[:, 7], data[:, 8]
    
    Hx = Hx*1e-3 #convert Hx to kA

    mx, my = data[:, 2], data[:, 3] 

    plt.rc('font', **font)
    fig = plt.figure(facecolor='w', figsize=(4, 3))
   
    plt.plot(Hx, mx, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='m_x', markersize=2)
    plt.plot(-Hx, -mx, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='', markersize=2)

    plt.grid(c='lightgrey')

    plt.ylabel(r"$m$")
    plt.xlabel("H (kA)")

    plt.tight_layout()
    plt.savefig('loop.png', dpi=300)

if __name__ == "__main__":
    plot_loop()