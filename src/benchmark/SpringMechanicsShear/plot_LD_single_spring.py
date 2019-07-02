#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from distutils.spawn import find_executable


if find_executable('latex'):
    print("latex installed")
    plt.rcParams.update({'text.usetex': True})

plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.linewidth': 2})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.serif': 'Times New Roman'})

# mm to inch
MTI = 0.0393700787

def plot_results_single_spring():
    names = ["displY [m]", "displX [m]", "loadY [N]", "loadX [N]"]
    # width = 147 * MTI  # textwidth on A4 
    # height = 200 * MTI
    # fig = plt.figure(figsize=(width, height))
    f, ((ax1, ax2, ax3, ax4)) = plt.subplots(4, sharex='col')  # , sharey='row')
    axs = [ax1, ax2, ax3, ax4]
    data = np.genfromtxt("LD.out", skip_header=1)
    steps = data[:, 0]
    for i, nam in enumerate(names):
        values = data[:, i+1]
        axs[i].plot(steps, values, 'k-', marker='+', markersize=3)
        # axs[i].axvline(x=40)
        axs[i].set_ylabel(nam)
        axs[i].yaxis.set_label_position("right")
    axs[-1].set_xlabel("time")
    f.subplots_adjust(hspace=.3)  # , wspace=.4)
    # plt.show()
    f.savefig("LD_single_spring.pdf")
    # fig.savefig("LD_single_spring.pdf")
    return


if __name__ == "__main__":
    plot_results_single_spring()
    print("DONE")
