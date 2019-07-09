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


def get_damage_data():
    damage = []
    cumSlip = []
    slips = []
    multip = []
    i = 1
    while True:
        file_name = "results/damage_%05d.out" % i
        if not os.path.isfile(file_name):
            break
        # print("exporting file %s" % file_name, end='')
        data = np.genfromtxt(file_name, usecols=(2, 3, 4, 5))
        damage.append(data[0, 0])
        cumSlip.append(data[0, 1])
        slips.append(data[0, 2])
        multip.append(data[0, 3])
        i += 1
        # print("damage = %lg, cumSlip = %lg, slip = %lg, multip = %lg" % (damage[-1], cumSlip[-1], slips[-1], multip[-1]))

    return [damage, cumSlip]


def plot_results_single_spring():
    names = ["displY [m]", "displX [m]", "loadY [N]", "loadX [N]"]
    data_values = []
    # width = 147 * MTI  # textwidth on A4
    # height = 200 * MTI
    # fig = plt.figure(figsize=(width, height))
    f, ((ax1, ax3, ax2, ax4, axD)) = plt.subplots(5, sharex='col')  # , sharey='row')
    axs = [ax1, ax2, ax3, ax4]
    data = np.genfromtxt("LD.out", skip_header=1)
    steps = data[:, 0]
    timo = data[:, 1]
    try:
        [damage, cumSlip] = get_damage_data()
    except Exception as e:
        print(e)
        print("but continue")

    for i, nam in enumerate(names):
        values = data[:, i+2]
        data_values.append(values)
        axs[i].plot(timo, values, 'k-', marker='o', markersize=0.5, mfc='r', mec='r', mew=0.5)
        # axs[i].axvline(x=40)
        axs[i].set_ylabel(nam)
        axs[i].yaxis.set_label_position("right")
    try:
        axD.plot(timo, damage, 'm-', marker='o', markersize=0.5, mfc='r', mec='r', mew=0.5, label="damage")
        # axD.plot(timo, cumSlip, 'c-', marker='o', markersize=0.5, mfc='g', mec='g', mew=0.5, label="cumSlip")
        axD.set_xlabel("time")
        axD.legend(loc='best', fontsize=6)
    except Exception as e:
        print(e)
        print("but continue")
    
    f.subplots_adjust(hspace=.3)  # , wspace=.4)
    # plt.show()
    f.savefig("LD_single_spring.pdf")
    plt.close()
    # fig.savefig("LD_single_spring.pdf")

    f, (ax5, ax6) = plt.subplots(2)

    ax5.plot(data_values[0], data_values[2], 'k-', marker='o', markersize=0.5, mfc='r', mec='r', mew=0.5)
    ax5.set_xlabel(names[0])
    ax5.set_ylabel(names[2])

    try:
        ax6.plot(timo, damage, 'm-', marker='o', markersize=0.5, mfc='r', mec='r', mew=0.5, label="damage")
        # ax6.plot(timo, cumSlip, 'c-', marker='o', markersize=0.5, mfc='g', mec='g', mew=0.5, label="cumSlip")
        ax6.set_xlabel("time")
        ax6.set_ylabel("damage")

        ax6.legend(loc='best')
    except Exception as e:
        print(e)
        print("but continue")

    f.subplots_adjust(hspace=.3)  # , wspace=.4)
    # plt.show()
    f.savefig("slipLoadDamage_single_spring.pdf")
    plt.close()

    f = plt.figure(figsize=(7, 5))
    ax = f.add_axes([0.1, 0.1, 0.85, 0.85])
    ax.plot(timo, data[:, 2], 'r-', lw=0.5)
    ax.plot([0, timo[-1]], [0, data[-1, 2]], 'k--')
    f.savefig("check_slip.pdf")
    plt.close()

    return


if __name__ == "__main__":
    plot_results_single_spring()
    print("DONE")
