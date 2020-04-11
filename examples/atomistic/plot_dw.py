import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

if __name__ == '__main__':

    fig, ax = plt.subplots(figsize=(4, 3))

    data = np.load("dw.npy")
    data.shape = (-1, 3)


    ax.plot(data[:, 0], '-', markersize=3, markeredgecolor='k', markerfacecolor='white', label="mx")
    ax.plot(data[:, 1], '-', markersize=3, markeredgecolor='k', markerfacecolor='white', label="my")
    ax.plot(data[:, 2], '-', markersize=3, markeredgecolor='k', markerfacecolor='white', label="mz")

    plt.legend()
    plt.xlabel("xs")
    plt.ylabel("m")
    fig.tight_layout()

    fig.savefig('dw.png')
