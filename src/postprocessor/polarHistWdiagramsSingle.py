import numpy as np
import matplotlib.pyplot as plt
import glob, os
import re
import anglesIO

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

folder = 'in1'
folder1 = 'in1'


N = 270
width = (2*np.pi) / N
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)

bounds = anglesIO.createBoundsList(6)

macroStep, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroEpsX, macroEpsY, macroGammaXY, macroPrincStress1, macroPrincStress2, macroPrincStrain1, macroPrincStrain2, macroPrincEnergy1, macroPrincEnergy2, macroEnergyX, macroEnergyY, macroEnergyXY, macroEnergy = anglesIO.loadMacroData (folder, 'PUCstrain_stress.out')
#availableBeamTimes = anglesIO.scanAvailableStepData(folder, macroTime)

availableBeamTimes = []

energyTotNA = []
energyTotTA = []

doneIdxA = []

idA = -666
idxA = -667
oldIdA = -66

displayStep = 1e-3
fileA = None
fileB = None
curdir = os.getcwd()

i = 0
for time in np.arange (0.0, np.amax(macroTime)*1.01, displayStep):
    idxA = int( macroStep[find_nearest(macroTime, time)] )

    #najit soubor  A s timhle cislem
    os.chdir(curdir)
    os.chdir(folder)
    oldFile = fileA
    for file in sorted(glob.glob("*.out")):
        if (file!='PUCstrain_stress.out'):
            str = [int(s) for s in re.findall(r'\d+', file)]
            idA = int(str[0])

            if (idA == idxA):
                fileA = file
                print('found A file: %s' %fileA )
                break
    if (oldFile == fileA):
        idxA= idA
        idA=oldIdA

    if (idA in doneIdxA):
        print('%d  already done, should skip...' %(idA))

    elif (idA == idxA ):
        print('idA %d idxA %d  ' %(idA, idxA))

        doneIdxA.append(idA)

        oldIdA = idA

        #input('press enter')
        stepTime = macroTime[idA]
        print(stepTime)
        availableBeamTimes.append(stepTime)


        os.chdir(curdir)
        os.chdir(folder)
        dx, dy, angle, dmgT, dmgN, slip, slipPi, strainN, strainT, energyTotNA, energyTotTA = anglesIO.loadBeamData(fileA, energyTotNA, energyTotTA)

        fig = anglesIO.assembleDataFigureSingle(angle,  N, idA-1,  bounds,  dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableBeamTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotNA, energyTotTA)


        print('Saving fig %s' % file)
        print()
        fig.savefig('out/'+file+'.png', bbox_inches='tight', dpi=100)
        i += 1
        plt.close()
        #"""
