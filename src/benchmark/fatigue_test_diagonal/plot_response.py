#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def plot_response():
    names = ["horizontal elongation", "vertical squeez", "vertical force"]
    fig = plt.figure(figsize=(7, 5))
    data = np.genfromtxt("results/LD.out", skip_header=1)
    timo = data[:, 1]
    for i in range(3):
        ax = fig.add_subplot(3, 1, i+1)
        values = data[:, i+2]
        ax.plot(timo, values, 'k-')
        ax.set_ylabel(names[i])
        if i != 2:
            ax.set_xticks([])
        else:
            ax.set_xlabel("time")
    fig.savefig("response.pdf")
    return


if __name__ == '__main__':
    plot_response()
    print("DONE")
