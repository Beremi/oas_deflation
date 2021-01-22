import numpy as np
import matplotlib.pyplot as plt
import os
from shutil import copy2
from distutils.spawn import find_executable
if find_executable('latex'):
    use_tex = True
    print("latex installed")
    custom_preamble = {
    "text.usetex": True,
    "text.latex.preamble": [
    r"\usepackage{amsmath}", # for the align enivironment
    ],
    }
    plt.rcParams.update(custom_preamble)
    plt.rcParams.update({'font.size': 10})
    plt.rcParams.update({'axes.linewidth': 2})
    plt.rcParams.update({'text.usetex': True})
    # plt.rcParams.update({'font.family': 'serif'})
    # plt.rcParams.update({'font.serif': 'Times New Roman'})
else:
    use_tex = False
    print("latex not-installed")
    print("will not plot out material parameters into plt")


MTI = 0.0393700787  # mm to inch


class Function(object):
    """docstring for Function."""
    name = "common"
    period = 0.01
    time_to_max = 1.0

    def __init__(self, ):
        super(Function, self).__init__()



def plot_results(ax, stress, strain, ls='-', col='b', normal=3,
                 print_x_label=False, print_y_label=False, label=None):
    if label is None:
        ax.plot(strain*1e2, np.array(stress)*1e-6, ls=ls, color=col)
    else:
        ax.plot(strain*1e2, np.array(stress)*1e-6, ls=ls, color=col, label=label)
    if normal == 0:
        if print_x_label:
            ax.set_xlabel(r"normal strain $\varepsilon_N$ [\%]")
        if print_y_label:
            ax.set_ylabel(r"normal stress $\sigma_N$ [MPa]")
    elif normal == 1:
        if print_x_label:
            ax.set_xlabel(r"shear strain $\varepsilon_T$ [\%]")
        if print_y_label:
            ax.set_ylabel(r"shear stress $\sigma_T$ [MPa]")
    elif normal == 2:
        if print_x_label:
            ax.set_xlabel(r"global strain $\varepsilon$ [\%]")
        if print_y_label:
            ax.set_ylabel(r"global stress $\sigma_{\mathrm{glob}}$ [MPa]")

    # offset = ax.yaxis.get_major_formatter().get_offset()
    # print("offset = %s" % offset)
    # print(type(offset))
    # if offset != 1.0:
    #     orig_label = ax.yaxis.get_label_text()
    #     print(orig_label)
    #     ax.yaxis.set_label_text("%s %lg" % (orig_label, offset))
    #     ax.yaxis.offsetText.set_visible(False)



    return


def genMechElemFile(folName="input_files"):
    elemFile = open("%s/mechElems.inp" % folName, 'w+')

    elemFile.write("# ElemType	nodeAidx	nodeBidx	nrOfVertices	vrtxAIdx	vrtxBIdx	Material\n")
    elemFile.write("LTCBEAM	0	1	2	11	10	0")

    elemFile.close()
    return


def genMechBCFile(folName="input_files"):
    bcFile = open("%s/mechBC.inp" % folName, 'w+')

    bcFile.write("# nodeIdx	KinTrX	KinTrY	KinRotZ	StTrX	StTrY	StRotZ\n")
    bcFile.write("NodalBC 0	0	0	0	-1	-1	-1\n")
    bcFile.write("NodalBC 1	1	2	0	-1	-1	-1")

    bcFile.close()
    return


def generateSolverFile(num_steps, folName="input_files"):
    solFile = open("%s/solver.inp" % folName, 'w+')

    solFile.write("SteadyStateNonLinearSolver\n")
    solFile.write("time_step %lg\n" % (1. / num_steps) )
    solFile.write("total_time 1.0")

    solFile.close()
    return


def generateFnFile(tensileLoad, shearLoad, folName="input_files", function=None):
    sym = 0
    fnFile = open("%s/functions.inp" % folName, 'w+')

    fnFile.write("PWLFunction 1 0 0\n")
    if function.name == "LinearFn":
        fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        fnFile.write("PWLFunction 2 0 1 0 %lg" % shearLoad)
    elif function.name == 'ConstSawToothFn':
        if tensileLoad - 1e-9 < 0:
            fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        else:
            fnFile.write("ConstSawToothFn value %lg period %lg sym %d\n"
                          % (tensileLoad, function.period, sym))
        fnFile.write("ConstSawToothFn value %lg period %lg sym %d\n"
                      % (shearLoad, function.period,  sym))
    elif function.name == 'LinSawToothFn':
        if tensileLoad - 1e-9 < 0:
            fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        else:
            fnFile.write("LinSawToothFn value %lg period %lg  time %lg sym %d\n"
                          % (tensileLoad, function.period, function.time_to_max, sym))
        fnFile.write("LinSawToothFn value %lg period %lg time %lg sym %d\n"
                      % (shearLoad, function.period, function.time_to_max, sym))
    else:
        print("function of type \'%s\' not implemented" % function)
        exit(1)
    fnFile.close()
    return

def generateNodeFile(length, perp_length=1.0, folName="input_files"):
    nodeFile = open("%s/nodes.inp" % folName, 'w+')

    nodeFile.write("Particle %lg %lg 0\n" % (0.5*length, 0.5*perp_length))
    nodeFile.write("Particle %lg %lg 0\n" % (1.5*length, 0.5*perp_length))

    nodeFile.write("AuxNode	%lg %lg\n" % (0.0*length, 0.5*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (0.5*length, 0.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (0.5*length, 1.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (2.0*length, 0.5*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (1.5*length, 1.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (1.5*length, 0.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (0.0*length, 0.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (0.0*length, 1.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (1.0*length, 0.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (1.0*length, 1.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg\n" % (2.0*length, 1.0*perp_length))
    nodeFile.write("AuxNode	%lg %lg" % (2.0*length, 0.0*perp_length))

    nodeFile.close()
    return


def genExpFile(length, perp_length=1.0, folName="input_files"):
    expFile = open("%s/exporters.inp" % folName, 'w+')

    expFile.write("DisplacementGauge LD displX ux %lg %lg %lg %lg\n" % (
        0.5*length, 0.5*perp_length, 1.5*length, 0.5*perp_length))
    expFile.write("DisplacementGauge LD displY uy %lg %lg %lg %lg\n" % (
        0.5*length, 0.5*perp_length, 1.5*length, 0.5*perp_length))

    expFile.write("ForceGauge LD loadX fx 1 1\n")
    expFile.write("ForceGauge LD loadY fy 1 1")

    expFile.close()
    return

def generateMaster(folName="input_files", ini_files="input_files",
                   matFile="materials.inp"):
    master_dst = "%s/master.inp" % folName
    masterFile = open(master_dst, 'w+')

    masterFile.write("Dimension	2\n")
    masterFile.write("Solver solver.inp\n")
    masterFile.write("NodeFiles	1 nodes.inp\n")

    if not os.path.isfile(matFile):
        print("material file %s not found" % matFile)
        exit(1)

    masterFile.write("MatFiles 1 ../%s\n" % matFile)

    masterFile.write("ElemFiles	1 mechElems.inp\n")
    masterFile.write("BCFiles 1 mechBC.inp\n")
    masterFile.write("FunctionFiles	1 functions.inp\n")
    masterFile.write("ExporterFiles	1 exporters.inp")

    masterFile.close()
    return master_dst


def generateInputs(tensileLoad, shearLoad,
                   length, perp_length, num_steps, folName,
                   matFile="materials.inp", function=None):

    genMechElemFile(folName)
    genMechBCFile(folName)
    generateSolverFile(num_steps, folName)
    generateFnFile(tensileLoad, shearLoad, folName, function=function)
    generateNodeFile(length, perp_length, folName)
    genExpFile(length, perp_length, folName)
    master = generateMaster(folName, matFile=matFile)
    return master


def runCalculation(angle, masterFile="input_files/master.inp",
                   binary="DiscreteModel",
                   nohup=True):
    if not os.path.isfile(binary):
        print("there is no binary file in this directory %s" % os.path.abspath("."))
        exit(1)
    if nohup:
        nohup_file = "nohup.out"
        cmd = "nohup ./%s %s > %s" % (binary, masterFile, nohup_file)

        if os.path.isfile(nohup_file):
            os.remove(nohup_file)
    else:
        cmd = "./%s %s" % (binary, masterFile)

    os.system(cmd)

    if nohup:
        print("calculation for angle %lg finished, for the course of the calculation, see %s"
              % (angle, nohup_file))

    return


def exportResults(folName="input_files"):
    LD_file = "%s/results/LD.out" % folName

    data = np.genfromtxt(LD_file, skip_header=1, usecols=(2, 3, 4, 5))

    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def run_and_export_single(angle_deg, final_displacement, num_steps, tension,
                          length, perp_length, folName,
                          axShear, axTensile, axCombined, nohup=False,
                          matFile="materials.inp", function=None):
    angle = angle_deg * np.pi / 180.
    cc = np.cos(angle)
    ss = np.sin(angle)
    shearLoad = final_displacement * cc
    tensileLoad = final_displacement * ss * (tension and 1 or -1)


    masterFile = generateInputs(tensileLoad, shearLoad,
                               length, perp_length, num_steps, folName,
                               matFile, function,
                               )

    runCalculation(angle_deg, masterFile, nohup=nohup)

    ###########################################################################

    displN, displT, loadN, loadT = exportResults(folName)


    stressN = np.append(0., loadN / (perp_length))
    stressT = np.append(0., loadT / (perp_length)) * (tension and 1 or -1)
    strainN = np.append(0., displN / length)
    strainT = np.append(0., displT / length) * (tension and 1 or -1)

    strains = np.sqrt( np.power(strainN, 2.) + np.power(strainT, 2.) ) * \
        (tension and 1 or -1)

    plot_results(axShear, stressT, strainT, col='b', normal=1,
                 # print_x_label=True
                 )

    plot_results(axTensile, stressN, strainN, col='r', normal=0,
                 # print_x_label=True
                 )

    # stressGlobal = np.sqrt( np.power(stressN, 2.) +
    #                        # alpha *
    #                        np.power(stressT, 2.) ) * (tension and 1 or -1)
    # plot_results(axCombined, stressGlobal, strains, col='g', normal=5,
    #              label=r"$\sigma_{\mathrm{glob}}$")
    project_N = ss * np.array(stressN)
    project_T = cc * np.array(stressT)
    stressProjected = project_N + project_T
    plot_results(axCombined, project_N, strains, col='r', normal=5,
                 label=r"$\sigma_{N,\mathrm{projected}}$")
    plot_results(axCombined, project_T, strains, col='b', normal=5,
                 label=r"$\sigma_{T,\mathrm{projected}}$")
    plot_results(axCombined, stressProjected, strains, col='k', normal=2,
                 # print_x_label=True,
                 label=r"$\sigma_{\mathrm{projected}}$")

    return




if __name__ == '__main__':
    ###########################################################################
    length = 30e-3  # element length
    perp_length = length  # facet size (in 2D it is just length)
    final_displacement = [1e-5, 1e-5]
    num_steps = 1000.
    # list of angles in degrees to plot curves for
    angle_deg_all = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., ]
    tension_both = [True, False]
    ###########################################################################
    # material file must be in the same directory as this python script
    ###########################################################################
    matFile = "material.inp"

    LinearFn = Function()
    LinearFn.name = "LinearFn"

    num_cycles = 100
    period = 1.0 / float(num_cycles)

    ConstSawToothFn = Function()
    ConstSawToothFn.name = "ConstSawToothFn"
    ConstSawToothFn.period = period

    LinSawToothFn = Function()
    LinSawToothFn.name = "LinSawToothFn"
    LinSawToothFn.period = period
    LinSawToothFn.time_to_max = 1.0

    ##########################################################################
    ### HERE CHOSE THE FUNCTION WHICH TO USE FOR LOADING #####################
    # function = LinearFn
    function = ConstSawToothFn
    # function = LinSawToothFn

    ##########################################################################
    ### HERE CHOSE VALUE OF FINAL DISPLACEMNT ################################
    if function.name == "LinearFn":
        final_displacement = [1e-3, 1e-3]
    elif function.name == "ConstSawToothFn":
        final_displacement = [5e-5, 5e-5]
    elif function.name == "LinSawToothFn":
        final_displacement = [1e-4, 1e-4]
    else:
        print("no such function defined %s" % function.name)

    matFileName = matFile.split(".")[0]
    folName = matFileName + "_calculation"
    if not os.path.isdir(folName):
        os.makedirs(folName)

    fig, axs = plt.subplots(len(angle_deg_all), 6,
                            figsize=(300. * MTI,
                                     50. * MTI * len(angle_deg_all)))


    ftsz = 16
    ftsz1 = 12
    for i, angle_deg in enumerate(angle_deg_all):
        for j, tension in enumerate(tension_both):
            if (len(angle_deg_all)==1):
                axShear = axs[(tension and 4 or 1)]
                axTensile = axs[(tension and 5 or 0)]
                axCombined = axs[(tension and 3 or 2)]
            else:
                axShear = axs[i, (tension and 4 or 1)]
                axTensile = axs[i, (tension and 5 or 0)]
                axCombined = axs[i, (tension and 3 or 2)]
            ###################################################################
            try:
                run_and_export_single(angle_deg, final_displacement[j],
                                      num_steps, tension,
                                      length, perp_length, folName,
                                      axShear, axTensile, axCombined,
                                      nohup=True, matFile=matFile,
                                      function=function)
            except Exception as e:
                try:
                    # this will probably not work on windows, therefore try catch is here
                    cpnohup_file = "nohup_angle_%lg.out" % angle_deg
                    os.system("cp nohup.out %s" % (cpnohup_file))
                    print("something went wrong during %s calculation under angle %lg, for details see %s"
                          % ((tension and "tensile" or "compressive"), angle_deg,
                             cpnohup_file))
                except Exception as ee:
                    pass
                print("Raised exception: ====================================")
                print(e)
                print("======================================================")

            ###################################################################
            if i == len(angle_deg_all)-1:
                axShear.set_xlabel(r"shear strain $\varepsilon_T$ [\%]",
                                   fontsize=ftsz1)
                axTensile.set_xlabel(r"normal strain $\varepsilon_N$ [\%]",
                                     fontsize=ftsz1)
                axCombined.set_xlabel(r"total strain $\varepsilon_T$ [\%]",
                                      fontsize=ftsz1)
            if i == 0:
                axCombined.set_title(tension and "tension" or "compression",
                                     fontsize=ftsz)
            if tension:
                axTensile.set_ylabel("%lg degrees" % angle_deg, fontsize=ftsz1)
                axTensile.yaxis.set_label_position("right")
            if i == len(angle_deg_all)//2 and not tension:
                axTensile.set_ylabel(
                    r"corresponding stress $\sigma$ [MPa]", fontsize=ftsz1,
                    labelpad=5)

    left_skip = 0.05
    fig.subplots_adjust(bottom=0.25/len(angle_deg_all),
                        top=1-(0.16/len(angle_deg_all)),
                        left=left_skip, right=0.98,
                        hspace=0.30, wspace=0.30)


    fig.savefig("%s_%s.pdf" % (matFileName, function.name))
    # plt.show()
    plt.close()
    try:
        os.system("rm nohup.out")
    except Exception as e:
        pass

    print("DONE")
