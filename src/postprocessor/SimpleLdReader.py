import numpy as np
import matplotlib.pyplot as plt
import sys, os, math

def printAvaiableData(file):
    f = open(file)
    header = f.readline()
    print('Available data:')
    for i in range(len(header.split('\t'))):
        if header[0]=='#':
            print('%d: %s' %(i,header.split('\t')[i]))
        else:
            print('No header. Follow your instinct. \nIf there are two columns, use ids 0 and 1.')
    return header

if __name__ == '__main__':
    #x nr_files file0 nr_plots0 xi yi ... xi yi filei nr_plotsi xi yi ... xi yi
    if (len(sys.argv) > 2 ):
        if ((len(sys.argv)-2)!=int(sys.argv[1])):
            nr_files = int(sys.argv[1])
            file_names = []
            files_plot_ids = []
            files_plot_multipliers = []
            headers = []

            j=0
            for i in range(nr_files):
                print('Loading %s' %sys.argv[2+j])
                file_names.append(sys.argv[2+j])
                f = open(sys.argv[2+j])
                h = f.readline()
                if h[0] != '#':
                    h='%sX\t%sY' %(sys.argv[2+j],sys.argv[2+j])
                headers.append(h)
                j+=1

                nr_plots_i = int(sys.argv[2+j])
                j+=1

                plot_ids = []
                plot_multipliers = []
                for p in range(nr_plots_i):
                    ix = int(sys.argv[2+j])
                    j+=1
                    mx= float(sys.argv[2+j])
                    j+=1
                    iy = int(sys.argv[2+j])
                    j+=1
                    my= float(sys.argv[2+j])
                    j+=1
                    print ('file %d graph #%d: x:%d / y:%d' %(i, p, ix, iy) )
                    plot_ids.append(np.array([ix,iy]))
                    plot_multipliers.append(np.array([mx,my]))
                files_plot_ids.append(plot_ids)
                files_plot_multipliers.append(plot_multipliers)

        else:
            nr_files = int(sys.argv[1])
            j=0
            print()
            for i in range(nr_files):
                file = sys.argv[2+j]
                print('Loading %s' %file)
                j+=1
                printAvaiableData(file)

            print('\nNow call me again like this: nr_files file0 nr_plots0 x0 mx0 y0 my0 ... xi mxi yi myi ... filei nr_plotsi x0 mx0 y0 my0 ... xi mxi yi myi')
            print('For reversing values, use negative indices.')
            sys.exit()

    else:
        print('Give me nr of files and filenames as arguments and call me again!')
        sys.exit()


    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2, left=0.2)
    titleX = ''
    titleY = ''
    plottedData = []
    maxima = []
    headerLine = ''
    fileHeaders = []

    for file in file_names:
        data = np.loadtxt(file)

        for pl in files_plot_ids[file_names.index(file)]:
            pl_id = files_plot_ids[file_names.index(file)].index(pl)
            ix = pl[0]
            iy = pl[1]

            dx = np.copy(data[:,np.abs(ix)])
            if (ix<0): dx *= -1
            dx *= files_plot_multipliers[file_names.index(file)][pl_id][0]
            titleX += '%s, ' %headers[file_names.index(file)].split('\t')[np.abs(ix)]
            headerLine += '%s\t' %headers[file_names.index(file)].split('\t')[np.abs(ix)]

            dy = np.copy(data[:,np.abs(iy)])
            if (iy<0): dy *= -1
            dy *= files_plot_multipliers[file_names.index(file)][pl_id][1]

            dyabs = np.abs(dy)
            ymax = dyabs.argmax()
            extr = np.array([dx[ymax], dy[ymax]])

            titleY += '%s, ' %headers[file_names.index(file)].split('\t')[np.abs(iy)]
            headerLine += '%s\t' %headers[file_names.index(file)].split('\t')[np.abs(iy)]

            tangent = dy[2] / dx[2]

            fileHeaders.append(headers[file_names.index(file)].split('\t')[np.abs(ix)]+'\t'+headers[file_names.index(file)].split('\t')[np.abs(iy)])

            pd=[]
            pd.append(dx)
            pd.append(dy)
            plottedData.append(pd)
            #ax.plot(dx, dy, label= 'X:%s / Y:%s \next: x=%.4e; y=%.4e / tg: %.4e' %(headers[file_names.index(file)].split('\t')[np.abs(ix)], headers[file_names.index(file)].split('\t')[np.abs(iy)], extr[0], extr[1], tangent))
            ax.plot(dx, dy, label= 'X:%s / Y:%s' %(headers[file_names.index(file)].split('\t')[np.abs(ix)], headers[file_names.index(file)].split('\t')[np.abs(iy)]))
            ax.scatter(extr[0], extr[1])
            ax.text(extr[0], extr[1], '%.3e'%extr[1])

    #ax.set(xlabel=titleX, ylabel=titleY)
    ax.grid()
    leg = ax.legend()
    ax.legend(loc=(-0.1,1.05))
    plt.tight_layout()
    fig.savefig("%s-%s.png" %(titleX, titleY))


    for i in range (len(plottedData)):
        fn = 'plottedData#%d.txt'%i
        fl=open(fn ,'w')
        pd = np.asarray(plottedData[i]).transpose()
        print(headerLine)
        np.savetxt(fl,  pd, delimiter='\t',   header = fileHeaders[i])
        fl.close()
        print('%s saved' %fn)
    #"""
    plt.show()



    """backup
    #nr_files file0 nr_plots0 xi yi filei nr_plotsi xi yi
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
    plottedData = []
    maxima = []
    headerLine = ''

    for pl in plotIdxs:
        ix = pl[0]
        iy = pl[1]

        dx = np.copy(data[:,np.abs(ix)])
        if (ix<0): dx *= -1
        titleX += '%s, ' %header.split('\t')[np.abs(ix)]
        headerLine += '%s\t' %header.split('\t')[np.abs(ix)]

        dy = np.copy(data[:,np.abs(iy)])
        if (iy<0): dy *= -1
        dyabs = np.abs(dy)
        ymax = dyabs.argmax()
        extr = np.array([dx[ymax], dy[ymax]])

        titleY += '%s, ' %header.split('\t')[np.abs(iy)]
        headerLine += '%s\t' %header.split('\t')[np.abs(iy)]

        tangent = dy[10] / dx[10]

        plottedData.append(dx)
        plottedData.append(dy)
        ax.plot(dx, dy, label= 'X:%s / Y:%s \next: x=%.4e; y=%.4e / tg: %.4e' %(header.split('\t')[np.abs(ix)], header.split('\t')[np.abs(iy)], extr[0], extr[1], tangent))
        ax.scatter(extr[0], extr[1])
        ax.text(extr[0], extr[1], '%.3e'%extr[1])

    ax.set(xlabel=titleX, ylabel=titleY)
    ax.grid()
    leg = ax.legend()
    ax.legend(loc=(-0.1,1.05))
    plt.tight_layout()
    fig.savefig("%s-%s.png" %(titleX, titleY))

    fl=open('plottedData.txt','w')
    plottedData = np.asarray(plottedData).transpose()
    #print(plottedData)
    np.savetxt(fl,  plottedData, delimiter='\t',   header = headerLine)
    fl.close()
    print('plottedData.txt saved')

    plt.show()

    #"""


























     #
