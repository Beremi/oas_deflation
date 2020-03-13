import numpy as np
import matplotlib.pyplot as plt
import sys, os, math

def printAvaiableData(file):
    f = open(file)
    header = f.readline()
    print('Available data:')
    for i in range(len(header.split('\t'))):
        print('%d: %s' %(i,header.split('\t')[i]))
    return header

if __name__ == '__main__':


    if len(sys.argv) > 2:
        file = sys.argv[1]
        print('Loading %s' %file)
        header = printAvaiableData(file)

        plotIdxs = []
        for i in range (2, (len(sys.argv))-1, 2):
            ix = int(sys.argv[i])
            iy = int(sys.argv[i+1])
            print ('graph #%d: x:%d / y:%d' %((i/2), ix, iy) )
            plotIdxs.append(np.array([ix,iy]))

    elif len(sys.argv) == 2:
        file = sys.argv[1]
        print('Loading %s' %file)
        printAvaiableData(file)
        print('Now call me again and add indices x0 y0 x1 y1 ... xi yi.')
        print('For reversing values, use negative indices.')
        sys.exit()

    else:
        print('Give me FILENAME as first argument and call me again!')
        sys.exit()




    data = np.loadtxt(file)

    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2, left=0.2)
    titleX = ''
    titleY = ''

    for pl in plotIdxs:
        ix = pl[0]
        iy = pl[1]

        dx = np.copy(data[:,np.abs(ix)])
        if (ix<0): dx *= -1
        titleX += '%s, ' %header.split('\t')[np.abs(ix)]

        dy = np.copy(data[:,np.abs(iy)])
        if (iy<0): dy *= -1
        titleY += '%s, ' %header.split('\t')[np.abs(iy)]

        ax.plot(dx, dy, label= 'X:%s / Y:%s' %(header.split('\t')[np.abs(ix)], header.split('\t')[np.abs(iy)]))


    ax.set(xlabel=titleX, ylabel=titleY)
    ax.grid()
    leg = ax.legend()
    fig.savefig("%s-%s.png" %(titleX, titleY))
    plt.show()




























     #
