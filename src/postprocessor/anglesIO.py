import numpy as np
import matplotlib.pyplot as plt
import glob, os
import re


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def createBoundsList (nrGraphs):
    bounds = []
    for i in range (nrGraphs):
        bounds.append(np.zeros(2))

    return bounds

#14 x y normal_x normal_y damageT damageN strainTY strainTZ strainPLTY strainPLTZ strainN strainT energy_totalN energy_totalT
####x y normal_x normal_y damageT damageN strainTY strainTZ strainPLTY strainPLTZ strainN strainT
def loadBeamData(file, energyTotN, energyTotT):
    print('Loading %s' %file)
    data = np.loadtxt(file)
    dx = data[:,2]
    dy = data[:,3]
    angle = np.arctan2(dy, dx)
    dmgT = data[:,4]
    dmgN = data[:,5]
    slip = data[:,6]
    slipPi = data[:,8]
    strainN = data[:,10]
    strainT = data[:,11]

    enTotN = np.sum(data[:,12])
    enTotT = np.sum(data[:,13])

    energyTotN.append(enTotN)
    energyTotT.append(enTotT)

    return dx, dy, angle, dmgT, dmgN, slip, slipPi, strainN, strainT,  energyTotN, energyTotT


def scanAvailableStepData(folder, macroTime):
    #os.chdir(folder)
    print('scanning...')
    availableBeamTimes = np.empty(1)
    for file in sorted(glob.glob("*.out")):
        if (file!='PUCstrain_stress.out'):
            str = [int(s) for s in re.findall(r'\d+', file)]
            print(str)
            stepNr = int(str[0])-1
            stepTime = macroTime[stepNr]
            #print(stepTime)
            availableBeamTimes.append(stepTime)

    return availableBeamTimes

def loadMacroData (folder, file):
    #step	time	sigma_x	sigma_y	tau_xy	eps_x	eps_y	gamma_xy
    macroscopicData = np.loadtxt(folder+'/'+file)
    macroscopicData = np.delete( macroscopicData, [len(macroscopicData)-1], axis=0)
    macroStep       = macroscopicData[:,0]
    macroTime       = macroscopicData[:,1]
    macroSigmaX     = macroscopicData[:,2]
    macroSigmaY     = macroscopicData[:,3]
    macroTauXY      = macroscopicData[:,4]
    macroEpsX       = macroscopicData[:,5]
    macroEpsY       = macroscopicData[:,6]
    macroGammaXY    = macroscopicData[:,7]

    macroPrincStress1 = (macroSigmaX + macroSigmaY)/2 + 0.5*np.sqrt(  (macroSigmaX-macroSigmaY) **2 + 4*macroTauXY**2    )
    macroPrincStress2 = (macroSigmaX + macroSigmaY)/2 - 0.5*np.sqrt(  (macroSigmaX-macroSigmaY) **2 + 4*macroTauXY**2    )

    macroPrincStrain1 = (macroEpsX + macroEpsY)/2 + 0.5*np.sqrt(  (macroEpsX-macroEpsY) **2 + 4*macroGammaXY**2    )
    macroPrincStrain2 = (macroEpsX + macroEpsY)/2 - 0.5*np.sqrt(  (macroEpsX-macroEpsY) **2 + 4*macroGammaXY**2    )

    macroPrincEnergy1 = np.zeros(len(macroscopicData))
    macroPrincEnergy2 = np.zeros(len(macroscopicData))

    macroEnergyX = np.zeros(len(macroscopicData))
    macroEnergyY = np.zeros(len(macroscopicData))
    macroEnergyXY = np.zeros(len(macroscopicData))

    for i in range (0,len(macroEnergyX),1):
        macroEnergyX[i] = np.trapz( macroSigmaX[0:i], macroEpsX[0:i] )
        macroEnergyY[i] = np.trapz( macroSigmaY[0:i], macroEpsY[0:i])
        macroEnergyXY[i] = np.trapz( macroTauXY[0:i], macroGammaXY[0:i] )
        macroPrincEnergy1[i] = np.trapz( macroPrincStrain1[0:i], macroPrincStress1[0:i] )
        macroPrincEnergy2[i] = np.trapz( macroPrincStrain2[0:i], macroPrincStress2[0:i] )

    macroEnergy =  macroEnergyX + macroEnergyY + macroEnergyXY

    return macroStep, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroEpsX, macroEpsY, macroGammaXY, macroPrincStress1, macroPrincStress2, macroPrincStrain1, macroPrincStrain2, macroPrincEnergy1, macroPrincEnergy2,macroEnergyX, macroEnergyY, macroEnergyXY, macroEnergy






def assemblePolarHistDmg(fig, angle, dmg, N, title, bounds, pos1, pos2, pos3):
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N
    n, bins = np.histogram(angle, bins = N)
    n1, bins1 = np.histogram(angle, bins = N, weights = dmg)
    dmgfin = n1 / n
    print (np.amax(dmgfin))
    radii, tick = np.histogram(angle, bins = N)
    ax1 = fig.add_subplot(pos1,pos2,pos3, polar=True)
    bars = ax1.bar(theta, dmgfin, width=width)
    ax1.title.set_text(title)
    min = 0
    max = np.amax(dmgfin)
    bounds[0] = np.minimum(bounds[0], min)
    bounds[1] = np.maximum(bounds[1], max)
    plt.ylim((bounds[0],bounds[1]))
    """
    plt.ylim((0,0.1))
    ticks = []
    ticks.append(0)
    ticks.append(0.1)
    labels = []
    labels.append ('')
    labels.append ('0.1')
    plt.yticks(ticks, labels)
    """
    for r, bar in zip(dmgfin, bars):
        bar.set_facecolor(plt.cm.Reds( r/bounds[1] ))
        bar.set_alpha(1)



def assemblePolarHistDmgDouble(fig, angleA, dmgA, N, titleA, boundsA, angleB, dmgB, titleB, boundsB, pos1, pos2, pos3):
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N

    n, bins = np.histogram(angleA, bins = N)
    n1, bins1 = np.histogram(angleA, bins = N, weights = dmgA)
    dmgfinA = n1 / n
    #radiiA, tickA = np.histogram(angleA, bins = N)

    n, bins = np.histogram(angleB, bins = N)
    n1, bins1 = np.histogram(angleB, bins = N, weights = dmgB)
    dmgfinB = n1 / n
    #radiiB, tickB = np.histogram(angleB, bins = N)

    ax1 = fig.add_subplot(pos1,pos2,pos3, polar=True)

    barsA = ax1.bar(theta, dmgfinA, width=width)
    barsB = ax1.bar(theta + np.pi/4, dmgfinB, width=width, edgecolor='black', color =(0,0,0,0))
    ax1.title.set_text(titleA +' / ' +titleB)

    min = 0
    max = np.amax(dmgfinA)
    boundsA[0] = np.minimum(boundsA[0], min)
    boundsA[1] = np.maximum(boundsA[1], max)

    min = 0
    max = np.amax(dmgfinB)
    boundsB[0] = np.minimum(boundsB[0], min)
    boundsB[1] = np.maximum(boundsB[1], max)

    for r, bar in zip(dmgfinA, barsA):
        bar.set_facecolor(plt.cm.Reds(0.5))# r/boundsA[1] ))
        bar.set_alpha(0.9)

    """
    for r, bar in zip(dmgfinB, barsB):
        bar.set_facecolor(plt.cm.PiYG( r/boundsB[1] ))
        bar.set_alpha(0.5)
    """
    plt.ylim((  0 , np.maximum(boundsA[1], boundsB[1]) ))



def assemblePolarHistGen(fig, angle, data, N, title, bounds, pos1, pos2, pos3):
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N
    n, bins = np.histogram(angle, bins = N)
    n1, bins1 = np.histogram(angle, bins = N, weights = data )
    datafin = n1 / n
    print (np.amax(datafin))

    #radii2, tick2 = np.histogram(angle, bins = N )
    ax2 = fig.add_subplot(pos1, pos2, pos3, polar=True)
    bars2 = ax2.bar(theta, abs(datafin), width=width)
    ax2.title.set_text(title)
    min = np.amin(datafin)
    max = np.amax(datafin)
    bounds[0] = np.minimum(bounds[0], min)
    bounds[1] = np.maximum(bounds[1], max)
    m = np.maximum( np.abs(bounds[0]), np.abs(bounds[1])  )
    bounds[0] = -m
    bounds[1] = m

    for r, bar in zip(datafin, bars2):
        bar.set_facecolor(plt.cm.seismic(  (r-bounds[0])/(bounds[1]-bounds[0])   ))
        bar.set_alpha(0.8)
    bds = np.maximum( abs(bounds[0]), abs(bounds[1]) )
    bounds[0] = 0
    bounds[1] = bds

    plt.ylim((bounds[0],bounds[1]))



def assemblePolarHistDouble(fig, angleA, dataA, N, titleA, boundsA, angleB, dataB, titleB, boundsB, pos1, pos2, pos3):
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N

    n, bins = np.histogram(angleA, bins = N)
    n1, bins1 = np.histogram(angleA, bins = N, weights = dataA)
    datafinA = n1 / n
    #radiiA, tickA = np.histogram(angleA, bins = N)

    n2, bins2= np.histogram(angleB, bins = N)
    n22, bins22 = np.histogram(angleB, bins = N, weights = dataB)
    datafinB = n22 / n2
    #radiiB, tickB = np.histogram(angleB, bins = N)

    ax1 = fig.add_subplot(pos1,pos2,pos3, polar=True)

    barsA = ax1.bar(theta, np.abs(datafinA), width=width)
    barsB = ax1.bar(theta+ np.pi/4, np.abs(datafinB), width=width, edgecolor='black', color =(0,0,0,0))
    ax1.title.set_text(titleA + ' / ' + titleB)

    min = np.amin(datafinA)
    max = np.amax(datafinA)
    boundsA[0] = np.minimum(boundsA[0], min)
    boundsA[1] = np.maximum(boundsA[1], max)
    mA = np.maximum( np.abs(boundsA[0]), np.abs(boundsA[1])  )
    boundsA[0] = -mA
    boundsA[1] = mA

    for r, bar in zip(datafinA, barsA):
        bar.set_facecolor(plt.cm.seismic(  (r-boundsA[0])/(boundsA[1]-boundsA[0])   ))
        bar.set_alpha(0.6)
    bds = np.maximum( abs(boundsA[0]), abs(boundsA[1]) )
    #boundsA[0] = 0
    #boundsA[1] = bds


    min = np.amin(datafinB)
    max = np.amax(datafinB)
    boundsB[0] = np.minimum(boundsB[0], min)
    boundsB[1] = np.maximum(boundsB[1], max)
    mB = np.maximum( np.abs(boundsB[0]), np.abs(boundsB[1])  )
    boundsB[0] = -mB
    boundsB[1] = mB

    #for r, bar in zip(datafinB, barsB):
        #bar.set_facecolor(plt.cm.PiYG(  (r-boundsB[0])/(boundsB[1]-boundsB[0])   ))
        #bar.set_alpha(1)
    bds = np.maximum( abs(boundsB[0]), abs(boundsB[1]) )
    #boundsB[0] = 0
    #boundsB[1] = bds

    #plt.ylim((  np.minimum(boundsA[0], boundsB[0])    , np.maximum(boundsA[1], boundsB[1]) ))




def assemblePolarHistGen(fig, angle, data, N, title, bounds, pos1, pos2, pos3):
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    width = (2*np.pi) / N
    n, bins = np.histogram(angle, bins = N)
    n1, bins1 = np.histogram(angle, bins = N, weights = data )
    datafin = n1 / n
    #print (np.amax(datafin))

    radii2, tick2 = np.histogram(angle, bins = N )
    ax2 = fig.add_subplot(pos1, pos2, pos3, polar=True)
    bars2 = ax2.bar(theta, abs(datafin), width=width)
    ax2.title.set_text(title)
    min = np.amin(datafin)
    max = np.amax(datafin)
    bounds[0] = np.minimum(bounds[0], min)
    bounds[1] = np.maximum(bounds[1], max)
    m = np.maximum( np.abs(bounds[0]), np.abs(bounds[1])  )
    bounds[0] = -m
    bounds[1] = m

    for r, bar in zip(datafin, bars2):
        bar.set_facecolor(plt.cm.seismic(  (r-bounds[0])/(bounds[1]-bounds[0])   ))
        bar.set_alpha(0.8)
    bds = np.maximum( abs(bounds[0]), abs(bounds[1]) )
    bounds[0] = 0
    bounds[1] = bds

    plt.ylim((bounds[0],bounds[1]))


def plotXYMovingGraph(subplt, label, x, y, stepNr, color, linestyle=None):
    #print(x)
    #print(y)
    subplt.plot(x, y, c=color,linestyle=linestyle,  linewidth=2, label=label+' [%e, %e]' %(x[stepNr], y[stepNr]))
    subplt.scatter(x[stepNr], y[stepNr],c=color, s=30)

def plotXYGraph(subplt, label, x, y, stepsArray, stepNr, color, linestyle=None):
    #print('steps')
    #print(stepsArray)
    #print('y')
    #print(y)
    subplt.plot(stepsArray, y, c=color,linestyle=linestyle,  linewidth=2, label=label+' [%e, %e]' %(stepsArray[len(y)-1], y[len(y)-1]))
    subplt.scatter(stepsArray[len(y)-1], y[len(y)-1], c=color, s=30)
    subplt.set_xlim((0, np.amax(x)))

def assembleStressesGraph(title, fig, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2, stepNr, pos1, pos2, pos3):
    stresses = fig.add_subplot(pos1,pos2,pos3)
    #stresses = plt.subplot(267)
    plotXYMovingGraph(stresses, 'SigmaX', macroTime, macroSigmaX, stepNr, 'red')
    plotXYMovingGraph(stresses, 'SigmaY', macroTime, macroSigmaY, stepNr, 'blue')
    plotXYMovingGraph(stresses, 'TauXY', macroTime, macroTauXY, stepNr, 'green')
    plotXYMovingGraph(stresses, 'Sigma1', macroTime, macroPrincStress1, stepNr, 'black', linestyle = ':')
    plotXYMovingGraph(stresses, 'Sigma2', macroTime, macroPrincStress2, stepNr, 'black', linestyle = ':')
    stresses.grid()
    stresses.set_title(title+' Macro Stresses')
    stresses.legend(loc='best')
    stresses.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

def assembleStrainsGraph(title, fig, macroTime, macroEpsilonX, macroEpsilonY, macroGammaXY, macroPrincStrain1,macroPrincStrain2, stepNr, pos1, pos2, pos3):
    strains = fig.add_subplot(pos1,pos2,pos3)
    #stresses = plt.subplot(267)
    plotXYMovingGraph(strains, 'EpsilonX', macroTime, macroEpsilonX, stepNr, 'red')
    plotXYMovingGraph(strains, 'EpsilonY', macroTime, macroEpsilonY, stepNr, 'blue')
    plotXYMovingGraph(strains, 'GammaXY', macroTime, macroGammaXY, stepNr, 'green')
    plotXYMovingGraph(strains, 'Strain1', macroTime, macroPrincStrain1, stepNr, 'black', linestyle = ':')
    plotXYMovingGraph(strains, 'Strain2', macroTime, macroPrincStrain2, stepNr, 'black', linestyle = ':')
    strains.grid()
    strains.set_title(title+' Macro Strains')
    strains.legend(loc='best')

def assembleEnergyGraph (title,fig, macroTime, availableTimes, macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotN, energyTotT, stepNr, pos1, pos2, pos3):
    energies = fig.add_subplot(pos1,pos2,pos3)
    #stresses = plt.subplot(267)
    """
    plotXYMovingGraph(energies, 'EnergyX', macroTime, macroEnergyX, stepNr, 'red')
    plotXYMovingGraph(energies, 'EnergyY', macroTime, macroEnergyY, stepNr, 'blue')
    plotXYMovingGraph(energies, 'EnergyXY', macroTime, macroEnergyXY, stepNr, 'green')

    plotXYMovingGraph(energies, 'Energy1', macroTime, macroPrincEnergy1, stepNr, 'black', linestyle = ':')
    plotXYMovingGraph(energies, 'Energy2', macroTime, macroPrincEnergy2, stepNr, 'black', linestyle = '--')
    plotXYMovingGraph(energies, 'EnergyTot', macroTime, macroPrincEnergy1+macroPrincEnergy2, stepNr, 'black', linestyle = '-.')
    """
    #print((energyTotN))
    #print((energyTotT))
    #print((energyTotN+energyTotT))
    plotXYGraph(energies, 'BeamE_N', macroTime, energyTotN, availableTimes, stepNr, 'red', linestyle = '--')
    plotXYGraph(energies, 'BeamE_T', macroTime, energyTotT, availableTimes, stepNr, 'black', linestyle = '--')
    plotXYGraph(energies, 'BeamE_Total', macroTime, (np.asarray(energyTotN)+np.asarray(energyTotT)), availableTimes, stepNr, 'blue', linestyle = ':')
    energies.grid()
    energies.set_title(title+' Energy')
    energies.legend(loc='best')


def assembleDataFigure(angle, N, stepNr, bounds,  dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2):
    fig = plt.figure(figsize=(30,10))
    #dmg t hist
    assemblePolarHistDmg(fig, angle, dmgT, N, 'DMG T', bounds[0], 2, 6, 1)
    #dmg n hist
    assemblePolarHistDmg(fig, angle, dmgN, N, 'DMG N', bounds[1], 2, 6, 2)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slip, N, 'Slip', bounds[2], 2, 6, 3)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slipPi, N, 'SlipPi', bounds[3], 2, 6, 4)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainN, N, 'Strain N', bounds[4], 2, 6, 5)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainT, N, 'Strain T', bounds[5], 2, 6, 6)

    #stresses graph
    stresses = assembleStressesGraph(fig, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2, stepNr, 2, 3, 4)
    #strains graph
    strains = assembleStrainsGraph(fig, macroTime, macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2, stepNr, 2, 3, 5)
    #strains graph
    energies = assembleEnergyGraph (fig, macroTime, availableTimes, macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, stepNr, 2, 3, 6)

    return fig



def assembleDualDataFigure(angle, Bangle, N, stepNrA, stepNrB, bounds,  dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, Bbounds, BdmgT, BdmgN, Bslip, BslipPi, BstrainN, BstrainT, BmacroTime, BavailableTimes, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1,BmacroPrincStress2 , BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1,BmacroPrincStrain2,BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1,BmacroPrincEnergy2, energyTotNA, energyTotTA,energyTotNB, energyTotTB):
    fig = plt.figure(figsize=(30,25))
    #dmg t hist
    assemblePolarHistDmg(fig, angle, dmgT, N, 'DMG T', bounds[0], 5, 6, 1)
    #dmg n hist
    assemblePolarHistDmg(fig, angle, dmgN, N, 'DMG N', bounds[1], 5, 6, 2)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slip, N, 'Slip', bounds[2], 5, 6, 3)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slipPi, N, 'SlipPi', bounds[3], 5, 6, 4)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainN, N, 'Strain N', bounds[4], 5, 6, 5)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainT, N, 'Strain T', bounds[5], 5, 6, 6)

    #dmg t hist
    assemblePolarHistDmg(fig, Bangle, BdmgT, N, 'DMG T', Bbounds[0], 5, 6, 7)
    #dmg n hist
    assemblePolarHistDmg(fig, Bangle, BdmgN, N, 'DMG N', Bbounds[1], 5, 6, 8)
    #dmg slip hist
    assemblePolarHistGen(fig, Bangle, Bslip, N, 'Slip', Bbounds[2], 5, 6, 9)
    #dmg slip hist
    assemblePolarHistGen(fig, Bangle, BslipPi, N, 'SlipPi', Bbounds[3], 5, 6, 10)
    #dmg slip hist
    assemblePolarHistGen(fig, Bangle, BstrainN, N, 'Strain N', Bbounds[4], 5, 6, 11)
    #dmg slip hist
    assemblePolarHistGen(fig, Bangle, BstrainT, N, 'Strain T', Bbounds[5], 5, 6, 12)


    #dmg t double
    assemblePolarHistDmgDouble(fig, angle, dmgT, N, 'DMGT/a', bounds[0], Bangle, BdmgT, 'DMGT/b', Bbounds[0], 5, 6, 13)
    #dmg n double
    assemblePolarHistDmgDouble(fig, angle, dmgN, N, 'DMGN/a', bounds[1], Bangle, BdmgN, 'DMGN/b', Bbounds[1], 5, 6, 14)
    #dmg slip double
    assemblePolarHistDouble(fig, angle, slip, N, 'Slip/a', bounds[2], Bangle, Bslip, 'Slip/b', Bbounds[2], 5, 6, 15)
    #dmg slippi double
    assemblePolarHistDouble(fig, angle, slipPi, N, 'SlipPI/a', bounds[3], Bangle, BslipPi, 'SlipPI/b', Bbounds[3], 5, 6, 16)
    #dmg strain n  double
    assemblePolarHistDouble(fig, angle, strainN, N, 'Strain N/a', bounds[4], Bangle, BstrainN, 'Strain N/b', Bbounds[4], 5, 6, 17)
    #dmg strain t double
    assemblePolarHistDouble(fig, angle, strainT, N, 'Strain T/a', bounds[5], Bangle, BstrainT, 'Strain T/b', Bbounds[5], 5, 6, 18)


    #stresses graph
    stresses = assembleStressesGraph(fig, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2, stepNrA, 5, 3, 10)
    #strains graph
    strains = assembleStrainsGraph(fig, macroTime, macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2, stepNrA, 5, 3, 11)
    #strains graph
    energies = assembleEnergyGraph (fig, macroTime, availableTimes, macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotNA, energyTotTA, stepNrA, 5, 3, 12)

    #stresses graph
    stresses = assembleStressesGraph(fig, BmacroTime, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1, BmacroPrincStress2, stepNrB, 5, 3, 13)
    #strains graph
    strains = assembleStrainsGraph(fig, BmacroTime, BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1, BmacroPrincStrain2, stepNrB, 5, 3, 14)
    #strains graph
    energies = assembleEnergyGraph (fig, BmacroTime, BavailableTimes, BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1, BmacroPrincEnergy2, energyTotNB, energyTotTB, stepNrB, 5, 3, 15)

    return fig

def assembleDataFigureSingle(angle,  N, stepNrA,  bounds,  dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotNA, energyTotTA):
    fig = plt.figure(figsize=(30,10))
    #dmg t hist
    assemblePolarHistDmg(fig, angle, dmgT, N, 'DMG T', bounds[0], 2, 6, 1)
    #dmg n hist
    assemblePolarHistDmg(fig, angle, dmgN, N, 'DMG N', bounds[1], 2, 6, 2)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slip, N, 'Slip', bounds[2], 2, 6, 3)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, slipPi, N, 'SlipPi', bounds[3], 2, 6, 4)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainN, N, 'Strain N', bounds[4], 2, 6, 5)
    #dmg slip hist
    assemblePolarHistGen(fig, angle, strainT, N, 'Strain T', bounds[5], 2, 6, 6)

    #stresses graph
    stresses = assembleStressesGraph('Stresses',fig, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2, stepNrA, 2, 3, 4)
    #strains graph
    strains = assembleStrainsGraph('Strains',fig, macroTime, macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2, stepNrA, 2, 3, 5)
    #strains graph
    energies = assembleEnergyGraph ('Energy', fig, macroTime, availableTimes, macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotNA, energyTotTA, stepNrA, 2, 3, 6)

    return fig


def assembleDualDataFigureSparse(angle, Bangle, N, stepNrA, stepNrB, bounds,  dmgT, dmgN, slip, slipPi, strainN, strainT, macroTime, availableTimes, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2 , macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2,macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, Bbounds, BdmgT, BdmgN, Bslip, BslipPi, BstrainN, BstrainT, BmacroTime, BavailableTimes, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1,BmacroPrincStress2 , BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1,BmacroPrincStrain2,BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1,BmacroPrincEnergy2, energyTotNA, energyTotTA,energyTotNB, energyTotTB):
    fig = plt.figure(figsize=(30,13))
    #dmg t double
    assemblePolarHistDmgDouble(fig, angle, dmgT, N, 'Sheared DMGT', bounds[0], Bangle, BdmgT, 'Principal DMGT', Bbounds[0], 3, 6, 1)
    #dmg n double
    assemblePolarHistDmgDouble(fig, angle, dmgN, N, 'Sheared DMGN', bounds[1], Bangle, BdmgN, 'Principal DMGN', Bbounds[1], 3, 6, 2)
    #dmg slip double
    assemblePolarHistDouble(fig, angle, slip, N, 'Sheared Slip', bounds[2], Bangle, Bslip, 'Principal Slip', Bbounds[2], 3, 6, 3)
    #dmg slippi double
    assemblePolarHistDouble(fig, angle, slipPi, N, 'Sheared SlipPI', bounds[3], Bangle, BslipPi, 'Principal SlipPI', Bbounds[3], 3, 6, 4)
    #dmg strain n  double
    assemblePolarHistDouble(fig, angle, strainN, N, 'Sheared Strain N', bounds[4], Bangle, BstrainN, 'Principal Strain N', Bbounds[4], 3, 6, 5)
    #dmg strain t double
    assemblePolarHistDouble(fig, angle, strainT, N, 'Sheared Strain T', bounds[5], Bangle, BstrainT, 'Principal Strain T', Bbounds[5], 3, 6, 6  )


    #stresses graph
    stresses = assembleStressesGraph('Sheared', fig, macroTime, macroSigmaX, macroSigmaY, macroTauXY, macroPrincStress1,macroPrincStress2, stepNrA, 3, 3, 4)
    #strains graph
    strains = assembleStrainsGraph('Sheared', fig, macroTime, macroEpsX, macroEpsY, macroGammaXY, macroPrincStrain1,macroPrincStrain2, stepNrA, 3, 3, 5)
    #strains graph
    energies = assembleEnergyGraph ('Sheared', fig, macroTime, availableTimes, macroEnergyX, macroEnergyY, macroEnergyXY, macroPrincEnergy1,macroPrincEnergy2, energyTotNA, energyTotTA, stepNrA, 3, 3, 6)

    #stresses graph
    stresses = assembleStressesGraph('Principal', fig, BmacroTime, BmacroSigmaX, BmacroSigmaY, BmacroTauXY, BmacroPrincStress1, BmacroPrincStress2, stepNrB, 3, 3, 7)
    #strains graph
    strains = assembleStrainsGraph('Principal', fig, BmacroTime, BmacroEpsX, BmacroEpsY, BmacroGammaXY, BmacroPrincStrain1, BmacroPrincStrain2, stepNrB, 3, 3, 8)
    #strains graph
    energies = assembleEnergyGraph ('Principal', fig, BmacroTime, BavailableTimes, BmacroEnergyX, BmacroEnergyY, BmacroEnergyXY, BmacroPrincEnergy1, BmacroPrincEnergy2, energyTotNB, energyTotTB, stepNrB, 3, 3, 9)

    return fig
