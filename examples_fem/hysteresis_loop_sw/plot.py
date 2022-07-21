
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
#plt.style.use(['science','no-latex', 'retro'])
#plt.style.use('seaborn-deep')
#mpl.rcParams['font.size'] = 13
#plt.style.use(['science', 'muted'])
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import numpy as np

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12,
        }

def custom_legend(legend):
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize(9)

    for label in legend.get_lines():
        label.set_linewidth(1)  # the legend line width

def plot_loop():

    #kOe = 79577.47154594767
    data = np.loadtxt('sw.txt')
    Hx, Hy = data[:, 7], data[:, 8]
    H = Hx/np.sqrt(2) + Hy/np.sqrt(2)

    H = H*1e-3 #convert H to kA

    mx, my = data[:, 2], data[:, 3] 

    m = mx/np.sqrt(2) + my/np.sqrt(2)

    for i in range(len(H)-1):
        if m[i]>0 and m[i+1]<0:
            print("Hc=", (H[i]+H[i+1])/2.0)

    plt.rc('font', **font)
    fig = plt.figure(facecolor='w', figsize=(4, 3))
    plt.plot(H, m, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='0K', markersize=1.2)
    plt.plot(-H, -m, linestyle='-', color='slateblue', marker='o', linewidth=1.2, label='', markersize=1.2)

    plt.grid(c='lightgrey')

    plt.ylabel(r"$m$")
    plt.xlabel("H (kA)")
    
    #l = axes.legend(loc='lower right')
    l = plt.legend()
    custom_legend(l)
    plt.tight_layout()
    fig.savefig('loop.png', dpi=600)

if __name__ == "__main__":
    plot_loop()