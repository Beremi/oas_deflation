import numpy as np
import matplotlib.pyplot as plt
import os


def plot_response():
    names = ["damage", "sPi", "slip", "EAlg", "alphaKin", "zIso", "tauTildaPiTrial", "tauPiTrial", "ftrial"]
    all_data = [[[], [], [], [], [], [], [], [], []], 
                [[], [], [], [], [], [], [], [], []], 
                [[], [], [], [], [], [], [], [], []], 
                [[], [], [], [], [], [], [], [], []], 
                [[], [], [], [], [], [], [], [], []], ]
    fig = plt.figure(figsize=(7, 12))
    step = 1
    steps = []
    while True:
        filename = "results/damage_%05d.out" % step
        if not os.path.exists(filename):
            break
        data = np.genfromtxt(filename, usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10))
        for j in range(5):
            for i in range(9):
                all_data[j][i].append(data[j, i])
        steps.append(step)
        step += 1

    colors = ['k', 'b', 'c', 'r', 'm']
    for i in range(9):
        ax = fig.add_subplot(9, 1, i+1)
        for j in range(5):
            values = all_data[j][i]
            ax.plot(steps, values, color=colors[j], ls='-', lw=0.1, label="eleme No.%d" % j)
        ax.set_ylabel(names[i])
        ax.yaxis.set_label_position("right")
        if i != 8:
            ax.set_xticks([])
        else:
            ax.set_xlabel("time")
    ax.legend(loc='best', fontsize=6)
    fig.savefig("all_elems.pdf")
    return


if __name__ == '__main__':
    plot_response()
    print("DONE")
