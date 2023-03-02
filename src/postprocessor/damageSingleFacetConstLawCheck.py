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
    symmetric = False
    lower = None

    def __init__(self, ):
        super(Function, self).__init__()



def plot_results(ax, strain, stress, ls='-', col='c',
                 x_label=None, y_label=None, label=None):
    if label is None:
        ax.plot(strain, np.array(stress), ls=ls, color=col)
    else:
        ax.plot(strain, np.array(stress), ls=ls, color=col, label=label)

    if x_label is not None:
        ax.set_xlabel(x_label)
    if y_label is not None:
        ax.set_ylabel(y_label)

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


def genMechBCFile(folName="input_files", force=False):
    bcFile = open("%s/mechBC.inp" % folName, 'w+')

    bcFile.write("# nodeIdx	KinTrX	KinTrY	KinRotZ	StTrX	StTrY	StRotZ\n")
    if force:
        bcFile.write("NodalBC 0	0	0	0	-1	-1	-1\n")
        bcFile.write("NodalBC 1	-1	-1	0	1	2	-1")
    else:
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
    sym = function.symmetric
    fnFile = open("%s/functions.inp" % folName, 'w+')

    fnFile.write("PWLFunction 1 0 0\n")
    if function.name == "LinearFn":
        fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        fnFile.write("PWLFunction 2 0 1 0 %lg" % shearLoad)
    elif function.name == 'ConstSawToothFn':
        if tensileLoad - 1e-9 < 0:
            fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        else:
            fnFile.write("ConstSawToothFn value %lg period %lg sym %d"
                          % (tensileLoad, function.period, int(sym)))
            if function.lower is not None:
                fnFile.write(" lower %lg" % (tensileLoad * function.lower))
            fnFile.write("\n")
        fnFile.write("ConstSawToothFn value %lg period %lg sym %d"
                      % (shearLoad, function.period,  int(sym)))
        if function.lower is not None:
            fnFile.write(" lower %lg" % (shearLoad * function.lower))
        fnFile.write("\n")
    elif function.name == 'LinSawToothFn':
        if tensileLoad - 1e-9 < 0:
            fnFile.write("PWLFunction 2 0 1 0 %lg\n" % tensileLoad)
        else:
            fnFile.write("LinSawToothFn value %lg period %lg  time %lg sym %d"
                          % (tensileLoad, function.period, function.time_to_max, int(sym)))
            if function.lower is not None:
                fnFile.write(" lower %lg" % (tensileLoad * function.lower))
            fnFile.write("\n")
        fnFile.write("LinSawToothFn value %lg period %lg time %lg sym %d"
                      % (shearLoad, function.period, function.time_to_max, int(sym)))
        if function.lower is not None:
            fnFile.write(" lower %lg" % (shearLoad * function.lower))
        fnFile.write("\n")
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
    expFile.write("ForceGauge LD loadY fy 1 1\n")
    expFile.write("ValueGauge LD damageN damageN\n")
    expFile.write("ValueGauge LD damageT damageT\n")
    expFile.write("ValueGauge LD nonElaStrainN strainPLN\n")
    expFile.write("ValueGauge LD nonElaStrainT strainPLT")

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
                   matFile="materials.inp", function=None, force=False):

    genMechElemFile(folName)
    genMechBCFile(folName, force)
    generateSolverFile(num_steps, folName)
    generateFnFile(tensileLoad, shearLoad, folName, function=function)
    generateNodeFile(length, perp_length, folName)
    genExpFile(length, perp_length, folName)
    master = generateMaster(folName, matFile=matFile)
    return master


def runCalculation(angle, masterFile="input_files/master.inp",
                   binary="OAS",
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

    try:
        data = np.genfromtxt(LD_file, skip_header=1)
    except Exception as e:
        print(e, end=' ...')
        print("but proceeding")
        data = np.genfromtxt(LD_file, skip_header=1, skip_footer=1)

    ##     time        displN      displT      loadN       loadT       damageN damageT
    return data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7], data[:, 8], data[:, 9]


def run_and_export_single(angle_deg, final_load, num_steps, tension, force,
                          length, perp_length, folName,
                          axDispl, axDamage, axSPi, axCombined, nohup=False,
                          matFile="materials.inp", function=None):
    angle = angle_deg * np.pi / 180.
    cc = np.cos(angle)
    ss = np.sin(angle)
    shearLoad = final_load * cc
    tensileLoad = final_load * ss * (tension and 1 or -1)


    masterFile = generateInputs(tensileLoad, shearLoad,
                               length, perp_length, num_steps, folName,
                               matFile, function, force,
                               )

    runCalculation(angle_deg, masterFile, nohup=nohup)

    ###########################################################################
    time_c, displN, displT, loadN, loadT, damageN, damageT, sPiN, sPiT = exportResults(folName)

    time_c = np.append(0., time_c)
    stressN = np.append(0., loadN / (perp_length)) * 1e-6
    stressT = np.append(0., loadT / (perp_length)) * (tension and 1 or -1) * 1e-6
    strainN = np.append(0., displN / length) * 1e2
    strainT = np.append(0., displT / length) * (tension and 1 or -1) * 1e2
    damageN = np.append(0., damageN)
    damageT = np.append(0., damageT)
    sPiN = np.append(0., sPiN)
    sPiT = np.append(0., sPiT)

    strains = np.sqrt( np.power(strainN, 2.) + np.power(strainT, 2.) ) * \
        (tension and 1 or -1)

    project_N = ss * np.array(stressN)
    project_T = cc * np.array(stressT)
    stressProjected = project_N + project_T
    plot_results(axCombined, strains, project_N, col='r',
                 label=r"$\sigma_{N,\mathrm{projected}}$")
    plot_results(axCombined, strains, project_T, col='c',
                 label=r"$\sigma_{T,\mathrm{projected}}$")
    plot_results(axCombined, strains, stressProjected, col='k',
                 # print_x_label=True,
                 label=r"$\sigma_{\mathrm{projected}}$")
    ###########################################################################
    plot_results(axDispl, time_c, strainN, col='r',
                 label=r"$\varepsilon_{N,\mathrm{projected}}$")
    plot_results(axDispl, time_c, strainT, col='c',
                 label=r"$\varepsilon_{T,\mathrm{projected}}$")
    plot_results(axDispl, time_c, strains, col='k',
                 # print_x_label=True,
                 label=r"$\varepsilon_{\mathrm{total}}$")
    ###########################################################################
    plot_results(axDamage, time_c, damageN, col='r',
                 label=r"$\omega_{N}$")
    plot_results(axDamage, time_c, damageT, col='c',
                 label=r"$\omega_{T}$")
    ###########################################################################
    plot_results(axSPi, time_c, sPiN, col='r',
                 label=r"$\varepsilon_{N}^{\pi}$")
    plot_results(axSPi, time_c, sPiT, col='c',
                 label=r"$\varepsilon_{T}^{\pi}$")
    return




if __name__ == '__main__':
    ###########################################################################
    length = 30e-3  # element length
    perp_length = length  # facet size (in 2D it is just length)
    num_steps = 1000.
    # list of angles in degrees to plot curves for
    angle_deg_all = [0., ]
    angle_deg_all = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90., ]
    tension = True  # for symmetric cyclic function does not matter, only for loading from zero to some value
    ## when force load applied, do not forget to set corresponding value of final_load
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # DO NOT USE YET !!!!
    force = False  # if True loading by force applied
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ###########################################################################
    # material file must be in the same directory as this python script
    ###########################################################################
    matFile = "material.inp"

    LinearFn = Function()
    LinearFn.name = "LinearFn"

    ###########################################################################
    ### HERE SPECIFY NUMBER OF CYCLES
    num_cycles = 10
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
    function = LinearFn
    function = LinSawToothFn
    function = ConstSawToothFn
    function.symmetric = False
    function.lower = 0.2

    ##########################################################################
    ### HERE CHOSE VALUE OF FINAL DISPLACEMNT ################################
    if function.name == "LinearFn":
        final_load = 1e-3  # if compression, put here the absolute value and set the tension var. to False (line 312)
    elif function.name == "ConstSawToothFn":
        final_load = 5e-5
    elif function.name == "LinSawToothFn":
        final_load = 1e-4
    else:
        print("no such function defined %s" % function.name)

    if force:  # set value of BC for force load (DO NOT USE YET !!!)
        final_load = 2.5e6

    matFileName = matFile.split(".")[0]
    folName = matFileName + "_calculation"
    if not os.path.isdir(folName):
        os.makedirs(folName)

    fig, axs = plt.subplots(len(angle_deg_all), 4,
                            figsize=(250. * MTI,
                                     50. * MTI * len(angle_deg_all)))


    ftsz = 16
    ftsz1 = 12
    for i, angle_deg in enumerate(angle_deg_all):

        if (len(angle_deg_all)==1):
            axCombined = axs[0]
            axDispl = axs[1]
            axDamage = axs[2]
            axSPi = axs[3]
        else:
            axCombined = axs[i, 0]
            axDispl = axs[i, 1]
            axDamage = axs[i, 2]
            axSPi = axs[i, 3]
        ###################################################################
        try:
            run_and_export_single(angle_deg, final_load,
                                  num_steps, tension, force,
                                  length, perp_length, folName,
                                  axDispl, axDamage, axSPi, axCombined,
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
                print(ee)
                pass
            print("Raised exception: ====================================")
            print(e)
            print("======================================================")

        ###################################################################
        if i == len(angle_deg_all)-1:
            axCombined.set_xlabel(r"total strain $\varepsilon$ [\%]",
                                  fontsize=ftsz1)
        if i == 0:
            axCombined.set_title(tension and "tension" or "compression",
                                 fontsize=ftsz)

        axSPi.set_ylabel("%lg degrees" % angle_deg, fontsize=ftsz1)
        axSPi.yaxis.set_label_position("right")
        if i == len(angle_deg_all)//2:
            try:
                axCombined.set_ylabel(
                    r"corresponding stress $\sigma$ [MPa]", fontsize=ftsz1,
                    labelpad=5)
            except Exception as e:
                print(e, end=' ...')
                print(" label setting not succesfull ... but proceed")

    left_skip = 0.05
    fig.subplots_adjust(bottom=0.25/len(angle_deg_all),
                        top=1-(0.16/len(angle_deg_all)),
                        left=left_skip, right=0.98,
                        hspace=0.30, wspace=0.30)


    fig.savefig("damage_%s_%s.pdf" % (matFileName, function.name))
    # plt.show()
    plt.close()
    try:
        os.system("rm nohup.out")
    except Exception as e:
        pass

    print("DONE")
