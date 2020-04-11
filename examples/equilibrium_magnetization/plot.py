import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

def coth(x):
    return (np.exp(2*x)+1)/(np.exp(2*x)-1)

def theory():
    Ms = 8e5
    V = 8e-27
    k_B = 1.380649e-23
    T = 100
    mu0 = 4*np.pi*1e-7
    mT = 0.001 / mu0
    c = mu0*Ms*V/(k_B*T)
    print("c=", c)
    H = np.linspace(1, 2000, 100)
    xi = H*c*mT

    #xi = np.linspace(1, 10, 10)
    print(xi)
    L = coth(xi) - 1/xi
    return H, L


if __name__ == '__main__':

    fig, ax = plt.subplots(figsize=(4, 3))

    data = np.loadtxt("M_H.txt")

    H, L = theory()
    ax.plot(H, L, "-", label="theory")
    ax.plot(data[:, 0], data[:, 1],'o', markersize=5, markeredgecolor='k', markerfacecolor='white', label="sim")

    plt.legend()
    plt.xlabel("H (mT)")
    plt.ylabel("m")
    fig.tight_layout()

    fig.savefig('m_H.png')
