import numpy as np
import matplotlib.pyplot as plt
import glob, os
import re
import anglesIO

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

folder = 'in'
folder1 = 'in1'


N = 2
width = (2*np.pi) / N
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)

bounds = anglesIO.createBoundsList(6)
Bbounds = anglesIO.createBoundsList(6)


macroStep, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroEpsX, macroEpsY, macroGammaXY, macroPrincStress1, macroPrincStress2, macroPrincStrain1, macroPrincStrain2, macroPrincEnergy1, macroPrincEnergy2, macroEnergyX, macroEnergyY, macroEnergyXY, macroEnergy = anglesIO.loadMacroData (folder, 'PUCstrain_stress.out')
#availableBeamTimes = anglesIO.scanAvailableStepData(folder, macroTime)

BmacroStep, BmacroTime, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStress1, BmacroPrincStress2, BmacroPrincStrain1, BmacroPrincStrain2, BmacroPrincEnergy1, BmacroPrincEnergy2, BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroEnergy = anglesIO.loadMacroData (folder1, 'PUCstrain_stress.out')
#BavailableBeamTimes = anglesIO.scanAvailableStepData(folder1, BmacroTime)


availableBeamTimes = []
BavailableBeamTimes = []

energyTotNA = []
energyTotTA = []
energyTotNB = []
energyTotTB = []

doneIdxA = []
doneIdxB = []

idA = -666
idxA = -667
oldIdA = -66
idB = -666
idxB = -667
oldIdB = -66

displayStep = 1e-3
fileA = None
fileB = None
curdir = os.getcwd()
m = np.maximum(np.amax(macroTime), np.amax(BmacroTime)) * 1.01
i = 0
for time in np.arange (0.0, m, displayStep):
    #print('time: %f' %time)
    #closest beam data index
    idxA = int( macroStep[find_nearest(macroTime, time)] )
    idxB = int( BmacroStep[find_nearest(BmacroTime, time)] )

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

    #najit soubor  B s timhle cislem
    os.chdir(curdir)
    os.chdir(folder1)
    oldFile = fileB
    for file in sorted(glob.glob("*.out")):
        if (file!='PUCstrain_stress.out'):
            str = [int(s) for s in re.findall(r'\d+', file)]
            idB = int(str[0])

            if (idB == idxB ):
                fileB = file
                print('found B file: %s' %fileB )
                break
    if (oldFile == fileB):
        idxB = idB
        idB=oldIdB

    if (idA in doneIdxA and idB in doneIdxB):
        print('%d and %d already done, should skip...' %(idA, idB))

    elif (idA == idxA or idB == idxB):
        print('idA %d idxA %d      idB %d idxB %d' %(idA, idxA, idB, idxB))

        doneIdxA.append(idA)
        doneIdxB.append(idB)

        oldIdA = idA
        oldIdB = idB

        #input('press enter')
        stepTime = macroTime[idA]
        print(stepTime)
        availableBeamTimes.append(stepTime)

        BstepTime = BmacroTime[idB]
        BavailableBeamTimes.append(BstepTime)

        os.chdir(curdir)
        os.chdir(folder)
        dx, dy, angle, dmgT, dmgN, slip, slipPi, strainN, strainT, energyTotNA, energyTotTA = anglesIO.loadBeamData(fileA, energyTotNA, energyTotTA)

        os.chdir(curdir)
        os.chdir(folder1)
        Bdx, Bdy, Bangle, BdmgT, BdmgN, Bslip, BslipPi, BstrainN, BstrainT, energyTotNB, energyTotTB = anglesIO.loadBeamData(fileB, energyTotNB, energyTotTB)

        print(availableBeamTimes)
        #fig = anglesIO.assembleDataFigure(angle, N, stepNr, bounds, dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2)
        #fig = anglesIO.assembleDualDataFigure(angle, Bangle, N, idxA-1, idxB-1, bounds, dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableBeamTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2,
        #Bbounds, BdmgT, BdmgN, Bslip, BslipPi, BstrainN, BstrainT, BmacroTime, BavailableBeamTimes, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1, BmacroPrincStress2 , BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1,BmacroPrincStrain2,BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1,BmacroPrincEnergy2, energyTotNA, energyTotTA, energyTotNB, energyTotTB)

        fig = anglesIO.assembleDualDataFigureSparse(angle, Bangle, N, idA-1, idB-1, bounds, dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableBeamTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2,
        Bbounds, BdmgT, BdmgN, Bslip, BslipPi, BstrainN, BstrainT, BmacroTime, BavailableBeamTimes, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1, BmacroPrincStress2 , BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1,BmacroPrincStrain2,BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1,BmacroPrincEnergy2, energyTotNA, energyTotTA, energyTotNB, energyTotTB)


        print('Saving fig %s' % file)
        print()
        fig.savefig('out/'+file+'.png', bbox_inches='tight', dpi=100)
        i += 1
        plt.close()
        #"""
