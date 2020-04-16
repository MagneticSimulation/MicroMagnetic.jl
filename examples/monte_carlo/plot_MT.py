import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


if __name__ == '__main__':

    fig, ax = plt.subplots(figsize=(4, 3))

    data = np.loadtxt("M_H.txt")

    ax.plot(data[:, 0], data[:, 1],'o-', markersize=5, markeredgecolor='k', markerfacecolor='white')

    plt.xlabel("Temperature (K)")
    plt.ylabel("m")
    fig.tight_layout()

    fig.savefig('mz_T2.png')
