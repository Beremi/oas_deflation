#!/usr/bin/env python3
import os
import sys
import pathlib
import argparse
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import arrow
from distutils.spawn import find_executable

if find_executable('latex'):
    print("latex installed")
    plt.rcParams.update({'text.usetex': True})

plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'axes.linewidth': 2})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.serif': 'Times New Roman'})

BASEDIR = pathlib.Path('.').resolve()

def loadNodes(filename):
    filename = filename
    nodes = np.loadtxt(BASEDIR.joinpath(filename), usecols=[1, 2], skiprows=1)
    return nodes


def loadVariables(filename, step):
    filename = BASEDIR.joinpath(filename + "_%05d.out" % step)
    return np.loadtxt(filename)


def plotData(ax, cax, nodes, values):

    cmap = matplotlib.cm.get_cmap('jet')
    norm = matplotlib.colors.Normalize(vmin=min(values), vmax=max(values))
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    ax.scatter(nodes[:, 0], nodes[:, 1], color=cmap(norm(values)))\
    #vmin = min(values)
    #vmax = max(values)
    #c = ax.scatter(nodes[:, 0], nodes[:, 1], cmap=cmap, vmin=vmin, vmax=vmax)

    #cbar = plt.colorbar(c, cax=cax, orientation="horizontal")

    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                                            orientation="horizontal")
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")

    tick_locator = matplotlib.ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()

    # ax.axis("equal")
    ax.axis("off")


######################################################################
def masterPlot(step, nodes, values, labels, xylim):

    print ("plotting step %d" % step)

    N = len(labels)
    axwidth = 8
    axdepth = 8/(xylim[0, 1]-xylim[0, 0])*(xylim[1, 1]-xylim[1, 0])
    textdepth = 0.5
    cbardepth = 0.5
    figwidth = N*axwidth
    figdepth = axdepth+textdepth+cbardepth
    fig = plt.figure(figsize=(figwidth, figdepth))
    for i in range(N):
        ax = fig.add_axes([i*axwidth/figwidth,
                           cbardepth/figdepth,
                           axwidth/figwidth,
                           axdepth/figdepth])
        ax.set_xlim(xylim[0, :])
        ax.set_ylim(xylim[1, :])
        cax = fig.add_axes([(i+0.1)*axwidth/figwidth,
                            cbardepth/figdepth/2,
                            0.8*axwidth/figwidth,
                            cbardepth/figdepth/2])
        plt.figtext((i+0.5)*axwidth/figwidth, 0.99, labels[i], fontsize=30,
                    ha="center", va="top")
        plotData(ax, cax, nodes[i], values[i])
    fig.savefig(OUTPUTDIR.joinpath("step_%04d.png" % step))
    plt.show()
    plt.close(fig)

def init_parser():
    # Create the parser
    parser = argparse.ArgumentParser(description='2D postprocessing',
                                     allow_abbrev=True)

    # Add the arguments
    parser.add_argument('folder',
                        metavar='folder',
                        # default=BASEDIR,
                        type=pathlib.Path,
                        # nargs='?',
                        help='the path to result folder')

    # Execute the parse_args() method
    return parser.parse_args()


if __name__ == "__main__":
    args = init_parser()

    BASEDIR = args.folder
    OUTPUTDIR = BASEDIR.joinpath('postprocessing')

    if not BASEDIR.exists():
        print('The path specified does not exist:', BASEDIR)
        sys.exit()

    nodes = loadNodes("nodes.inp")
    aux = loadNodes("auxNodes.inp")
    vertices = loadNodes("vertices.inp")
    OUTPUTDIR.mkdir(parents=True, exist_ok=True)
    xylim = np.zeros([2, 2])
    xylim[0, :] = [min(vertices[:, 0]), max(vertices[:, 0])]
    xylim[1, :] = [min(vertices[:, 1]), max(vertices[:, 1])]
    dist = xylim[:, 1]-xylim[:, 0]
    xylim[:, 0] = xylim[:, 0] - 0.05*dist
    xylim[:, 1] = xylim[:, 1] + 0.05*dist

    for step in range(1, 6):
        translations = loadVariables("translations", step)[0:len(nodes), :]
        pressure = loadVariables("pressure", step)[len(nodes)+len(aux): len(nodes)+len(aux)+len(vertices)]
        masterPlot(step, [nodes, nodes, vertices], [translations[:, 0]*1000.,
                   translations[:, 1]*1000., pressure],
                   [r"$u_x$ [mm]", r"$u_y$ [mm]", "pressure [Pa]"], xylim)
