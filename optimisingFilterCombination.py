#!/usr/bin/env python

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import integrate
import attenuation as atten
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import sys,IP

""" Element density dictionary """
densities = {'Ag':10.5,'Al':2.2,'Au':19.3,'C':2.27,'Cu':8.96,'Fe':4.54,'Mn':7.21,'Mo':10.22,'Mylar':1.38,
    'Pb':11.35,'Polystyrene':1.06,'Sn':7.31,'Ta':16.65,'Teflon':2.2,'Ti':4.506,'W':19.25,'Zn':7.14}

""" Load in image plate data and calculate its properties; sensitivity and attenuation """
def loadIPdata(hybrid):
    if(hybrid):
        IPdata        = np.loadtxt("IPsensitivity_PSLperPhoton.txt")
        E_MeV         = IPdata[:,0]
        E_keV         = E_MeV*1e3
        IPsensitivity = IPdata[:,1]
    else:
        E_keV = np.arange(6,3001,1)
        E_MeV = E_keV * 1e-3
    surface_layer  = IP.layer(E_MeV,'Surface', 9.,  1.66,['Mylar'],[1.])
    phosphor_layer = IP.layer(E_MeV,'Phosphor',115.,3.31,['Bromine','Fluorine','Barium','Iodine'],[0.85,1,1,0.15])
    back_layer     = IP.layer(E_MeV,'Back',    190.,1.66,['Mylar'],[1.])
    ferrite_layer  = IP.layer(E_MeV,'Ferrite', 160.,2.77,['Zn','Mn','Fe','O','H','C'],[1.,2.,5.,40.,15.,10.])
    layers         = [surface_layer,phosphor_layer,back_layer,ferrite_layer]
    imagePlateAttenuation = 1.
    for layer in layers:
        imagePlateAttenuation *= layer.attenuation
    # for layer in layers:
    #     plt.plot(E_keV,np.ones(len(layer.attenuation))-layer.attenuation,label=layer.name)
    # plt.xlim([0,20])
    # plt.legend(loc="best")
    # plt.show()
    if not(hybrid):
        IPsensitivity = phosphor_layer.sensitivity
    return E_MeV, IPsensitivity, imagePlateAttenuation

""" Load in filter data and calculate the IP response curves (combination of filter atten and IP response) """
def loadFilterData(E_MeV,preFilter=False):
    filters = np.loadtxt('./responseFunctions/'+fileName+'/'+fileName+'_FILTERS.txt',dtype=str,delimiter=',')
    filterMaterials, filterThicknesses = filters[:,0], [float(i) for i in filters[:,1]]
    numChannels = len(filterMaterials)
    numFilters  = np.arange(1,numChannels+1,1)
    if(cannonType=='singleChannel'):
        responseFunctions = []
        totalAttenuation = 1.
        for filMat, filThick, filt in zip(filterMaterials, filterThicknesses, numFilters):
            filterAttenuation = atten.attenuation(E_MeV,filMat,densities[filMat],filThick)
            totalAttenuation *= filterAttenuation
            if(preFilter and filt==1):continue
            IPresponse = totalAttenuation*IPsensitivity
            totalAttenuation *= imagePlateAttenuation
            # totalAttenuation = 1.
            responseFunctions.append(IPresponse)
    elif(cannonType=='multiChannel'):
        responseFunctions = []
        for filMat, filThick in zip(filterMaterials, filterThicknesses):
            filterAttenuation = atten.attenuation(E_MeV,filMat,densities[filMat],filThick)
            IPresponse = filterAttenuation*IPsensitivity
            responseFunctions.append(IPresponse)
    else:
        print 'Cannon design not specified!'
    if(cannonType=='singleChannel' and preFilter):
        ch = np.arange(1,numChannels,1)
    else:
        ch = np.arange(1,numChannels+1,1)
    return filterMaterials, filterThicknesses, responseFunctions, ch

""" Plotting for the IP response curves """
def plotResponseFunctions(E_MeV,responseFunctions,filterMaterials,filterThicknesses,figNum,fileName,figName,log=True,save=False):
    mpl.rcParams['font.size']=24
    plt.figure(figNum)
    colours = cm.rainbow(np.linspace(0, 1, len(filterMaterials)))
    for resFunc,filMat,filThick,colour in zip(responseFunctions,filterMaterials,filterThicknesses,colours):
        if(normed):
            plt.plot(E_MeV*1e3,resFunc*(1./np.max(resFunc[E_MeV<0.5])),label='%s %.1f'%(filMat,filThick),color=colour)
        else:
            plt.plot(E_MeV*1e3,resFunc,label='%s %.1f'%(filMat,filThick),color=colour)
    plt.legend(loc='best',prop={"size":12},ncol=2)
    plt.xlabel('Photon Energy (keV)')
    plt.ylabel('IP Response (PSL/photon)')
    # plt.xlim([0,250])
    plt.xlim([1,500])
    # plt.ylim([1e-5*max(responseFunctions[0]),5e0*max(responseFunctions[0])])
    plt.tight_layout()
    if(log):
        plt.yscale('log')
        # plt.ylim([1e-10*max(responseFunctions[0]),5e0*max(responseFunctions[0])])
    # plt.tight_layout()
    if(save):
        plt.savefig('./responseFunctions/'+fileName+'/'+figName)
    return figNum+1

""" Saving the response functions """
def saveResFuncs(E_MeV,ch,fileName,txtName,save=False):
    resFuncOutput = np.zeros((len(E_MeV),len(ch)+1))
    for resFunc,chan in zip(responseFunctions,ch):
        resFuncOutput[:,chan] = resFunc
    resFuncOutput[:,0] = E_MeV
    if(save):
        np.savetxt('./responseFunctions/'+fileName+'/'+txtName,resFuncOutput,delimiter=',',fmt='%.5e')
    return

""" Function to generate the photon number spectrum dN/dE, N.B. E and T must be in same unit """
def funcPhotonNumSpec(E,A,T):
    photonNumSpec = A*(1./E)*np.exp(-E/T)
    photonNumSpec[E>10.*T] = 0.
    return photonNumSpec
# def funcPhotonNumSpec(E_keV,A,T_keV):
#     """ Zeff is the effective nuclear charge """
#     Zeff = 3.5
#     photonSource = A*np.exp(1.-E_keV/T_keV)*(5e11/(4.*np.pi))*(Zeff/79.)
#     return photonSource

""" Calculate the IP responses """
def imagePlateResponses(E_MeV,A,T_MeV,responseFunctions,ch):
    IPsignals = []
    for resFunc in responseFunctions:
        photonNumSpec = funcPhotonNumSpec(E_MeV,A,T_MeV)
        measuredSpectrum = resFunc*photonNumSpec
        IPsignal = integrate.trapz(measuredSpectrum, E_MeV)
        IPsignals.append(IPsignal)
    if(plotting):
        mpl.rcParams['font.size']=24
        plt.figure(figNum)
        plt.plot(ch,IPsignals,marker='+',linestyle=' ',markersize=12, color='red',zorder=14)
        plt.yscale('log')
        plt.xlabel('Channel Number')
        plt.ylabel('Channel Signal (mPSL)')
        plt.xticks((ch))
        # plt.show()
    return np.array(IPsignals)

""" Find the closest value in an array """
def find_closest(data, v):
	return (np.abs(np.array(data)-float(v))).argmin()

""" ############# """
""" Gen Res Funcs """
""" ############# """
if(1):
    hybrid   = False
    normed   = False
    plotting = True

    fileName = 'BMXS'
    figName  = fileName+'.png'
    if(hybrid):
        txtName = fileName+'_hybrid.txt'
    else:
        txtName = fileName+'.txt'
    saving = False
    cannonType = 'singleChannel'
    # cannonType = 'multiChannel'
    """ Load in image plate data and calculate its properties; sensitivity and attenuation """
    E_MeV, IPsensitivity, imagePlateAttenuation = loadIPdata(hybrid) #low res produces manageable file sizes
    """ Load in filter data and calculate the IP response curves (combination of filter atten and IP response """
    filterMaterials, filterThicknesses, responseFunctions, ch = loadFilterData(E_MeV,preFilter=True)
    """ Plot the response curves for the different filter thicknesses """
    figNum = 1
    figNum = plotResponseFunctions(E_MeV,responseFunctions,filterMaterials,filterThicknesses,figNum,fileName,figName,log=False,save=saving)
    """ Saving the response functions """
    # saveResFuncs(E_MeV,ch,fileName,txtName,save=saving)
    """ Show the graphs """
    plt.show()

if(1):
    """ Set B/s params """
    A_base = 100.0
    T_base = 40.0
    T_MeV  = T_base*1e-3
    E_keV  = E_MeV*1e3

    """ Input deck for scanning parameters """
    tempArray = np.linspace(10,70,100)
    numTemps  = round(len(tempArray))

    """ Construct coarse laser energy array """
    fracArray = np.linspace(1.0,200.0,100)
    numFracs = round(len(fracArray))

    fittedTemps = []
    fittedEngys = []
    tempErrors  = []
    engyErrors  = []
    np.random.seed(0)
    for i in range(250):
        print i
        """ Make fake signals and add random noise """
        IPsigsRaw = imagePlateResponses(E_MeV,A_base,T_MeV,responseFunctions,ch)
        # IPsignals = np.array([IPsig for IPsig in IPsigsRaw])
        IPsignals = np.array([IPsig+np.random.normal(0.0, 0.1*IPsig)+np.abs(np.random.normal(0.0,1e-4)) for IPsig in IPsigsRaw])
        stdDev    = 0.1*IPsignals+1e-4

        chiSqBest   = 1e99
        chiSqMatrix = np.zeros((int(numFracs),int(numTemps)))
        for T_keV,t in zip(tempArray,range(int(numTemps))):
            photonSource     = funcPhotonNumSpec(E_keV,1.0,T_keV)
            signals          = photonSource*responseFunctions
            measurementArray = integrate.trapz(signals,E_keV)
            for frac,N in zip(fracArray,range(int(numFracs))):
                compSignalAbs    = measurementArray*frac
                chiSq            = sum(((compSignalAbs-IPsignals)/stdDev)**2.)
                chiSqMatrix[N,t] = chiSq
        chiSqMin   = chiSqMatrix.min()
        minChSqInd = np.unravel_index(np.argmin(chiSqMatrix,axis=None),chiSqMatrix.shape)
        bestFrac   = fracArray[minChSqInd[0]]
        bestTemp   = tempArray[minChSqInd[1]]

        """ Best signals based on the fits """
        bestSignals = integrate.trapz(funcPhotonNumSpec(E_keV,1.0,bestTemp)*responseFunctions,E_keV)*bestFrac

        if(plotting):
            mpl.rcParams['font.size']=24
            plt.figure(figNum,figsize=(9,6.5))
            figNum += 1
            plt.plot(ch,IPsigsRaw,  marker='x',linestyle=' ',markersize=12, color='red',  zorder=14, label="Raw Signals")
            plt.plot(ch,IPsignals,  marker='+',linestyle=' ',markersize=12, color='blue', zorder=14, label="Rand Signals")
            plt.plot(ch,bestSignals,marker='.',linestyle=' ',markersize=12, color='green',zorder=14, label="Best Signals")
            # plt.ylim(bottom=1e-4)
            plt.yscale('log')
            plt.xlabel('Channel Number')
            plt.ylabel('Channel Signal (mPSL)')
            plt.xticks((ch))
            plt.legend(loc="best")
            plt.tight_layout()
            plt.show()

        """ Finding the minimum chi squared value """
        numberParameters = 2.
        numberDataPoints = len(ch)
        chiSqreduction   = 1./(numberDataPoints-numberParameters)
        if(0): #This section normalises the chi squareds to the minimum (setting it to 1) which helps with the errors
            chiSqMatrix /= chiSqMatrix.min()
            chiSqMin     = 1.
        else:
            chiSqMatrix *= chiSqreduction
            chiSqMin     = chiSqMatrix.min()
        minChSqInd   = np.unravel_index(np.argmin(chiSqMatrix,axis=None),chiSqMatrix.shape)
        bestFrac     = fracArray[minChSqInd[0]]
        bestTemp     = tempArray[minChSqInd[1]]

        """ Find the best fitting T for each value of n and plots the chi squared """
        chiSqVariedn = []
        chiSqAmp     = []
        for i in range(len(chiSqMatrix[0,:])):
            chiSqCol = chiSqMatrix[:,i]
            chiSqVariedn.append(chiSqCol.min())
            chiSqAmp.append(fracArray[chiSqCol.argmin()])

        chiSqVariedT = []
        chiSqTemp = []
        for i in range(len(chiSqMatrix[:,0])):
            chiSqCol = chiSqMatrix[i,:]
            chiSqVariedT.append(chiSqCol.min())

        try:
            crazyFig = True
            """ Finding the errors """
            fracMinusErr = fracArray[find_closest(chiSqVariedT[:minChSqInd[0]],chiSqMin+2.3)]
            fracPlusErr  = fracArray[find_closest(chiSqVariedT[minChSqInd[0]:],chiSqMin+2.3)+minChSqInd[0]]
            tempMinusErr = tempArray[find_closest(chiSqVariedn[:minChSqInd[1]],chiSqMin+2.3)]
            tempPlusErr  = tempArray[find_closest(chiSqVariedn[minChSqInd[1]:],chiSqMin+2.3)+minChSqInd[1]]
            ampArray,ampMinusErr,ampFrac,ampPlusErr = fracArray,fracMinusErr,bestFrac,fracPlusErr
            tempErr = ((tempPlusErr-bestTemp)+(bestTemp-tempMinusErr))/2.
            ampErr  = ((ampPlusErr-ampFrac)+(ampFrac-ampMinusErr))/2.
            ylabel = "E$_{hot}$ (J)"#'Total Hot Electron Energy [J]'
            # print 'Best temp. [keV] = %.1f, + %.1f, - %.1f'%(bestTemp,tempPlusErr-bestTemp,bestTemp-tempMinusErr)
            # print 'Best engergy [J] = %.2f, + %.2f, - %.2f'%(ampFrac,(ampPlusErr-ampFrac),(ampFrac-ampMinusErr))
            # print ' '
        except Exception as ExceptErr:
            crazyFig = False
            print ExceptErr
            print "Could not find minimum!"

        fittedTemps.append(bestTemp)
        fittedEngys.append(bestFrac)
        tempErrors.append(tempErr)
        engyErrors.append(ampErr)
        # continue

        """ Crazy figure for paper """
        if(crazyFig and plotting):
            mpl.rcParams['font.size']=18
            fig5    = plt.figure(figNum,figsize=(10,10))#,constrained_layout=True)
            widths  = [1, 4]
            heights = [4, 1]
            spec5   = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=widths,height_ratios=heights,figure=fig5)
            spec5.update(wspace=0.1, hspace=0.1)
            ax = fig5.add_subplot(spec5[0, 0])
            plt.plot(chiSqVariedT,ampArray,linewidth=3)
            plt.plot([chiSqMin+2.3,chiSqMin+2.3,chiSqMin+2.3],np.array([ampMinusErr,ampFrac,ampPlusErr]),linestyle='-',color='k',linewidth=3,marker='_',markersize=10,mew=2)
            plt.ylim([min(ampArray),max(ampArray)])
            # plt.xlim([0,5])
            plt.xscale('log')
            plt.ylabel(ylabel)
            plt.xlabel('$\widetilde{\chi}^{2}}$')
            # plt.yscale('log')
            # plt.axis('off')

            ax   = fig5.add_subplot(spec5[0, 1])
            X, Y = np.meshgrid(tempArray,ampArray)
            im   = plt.pcolormesh(X,Y,np.log(chiSqMatrix),cmap='viridis')
            bestFit, = plt.plot(bestTemp,ampFrac,linestyle="",marker="x",mew=16,ms=3,color="white",label=r"%i\pm%i, %i\pm%i"%(bestTemp,tempErr,ampFrac,ampErr))
            # plt.plot(tempArray,chiSqAmp, linewidth=2,linestyle='--',color='grey')
            # cbar = plt.colorbar()#(pad=0.)
            # cbar.set_label('log($\widetilde{\chi}^{2}}$)')
            contourSet1 = plt.contour(X,Y,chiSqMatrix,chiSqMin+1.0,colors='yellow')
            contourSet2 = plt.contour(X,Y,chiSqMatrix,chiSqMin+2.3,colors='orange')
            contourSet3 = plt.contour(X,Y,chiSqMatrix,chiSqMin+4.6,colors='red')
            labels = ["$\widetilde{\chi}^{2}_{min} + 1.0$ (25%)", "$\widetilde{\chi}^{2}_{min} + 2.3$ (68%)","$\widetilde{\chi}^{2}_{min} + 4.6$ (95%)", "T$_{hot}$ = %i$\pm$%i keV\nE$_{hot}$ = %i$\pm$%i J"%(bestTemp,np.ceil(tempErr),ampFrac,np.ceil(ampErr))]
            h1,l1 = contourSet1.legend_elements()
            h2,l2 = contourSet2.legend_elements()
            h3,l3 = contourSet3.legend_elements()
            # h4,l4 = plt.gca().legend_elements(bestFit)
            plt.legend([h1[0], h2[0], h3[0], bestFit], [labels[0], labels[1], labels[2], labels[3]], loc="upper right")
            # plt.legend(loc="upper right")
            plt.tick_params(axis='both',which='both',labelbottom=False,labeltop=False,labelleft=False,labelright=False)#,left=False,right=False,bottom=False,top=False)
            # plt.tick_params(axes='y',which='both',left=False,right=False,labelleft=False,labelright=False)
            # plt.xlabel('T$_{hot}$ [keV]')
            # plt.ylabel(ylabel)
            # plt.yscale('log')
            cb_ax = fig5.add_axes([.91,.2925,.02,0.5875])
            cbar  = fig5.colorbar(im,orientation='vertical',cax=cb_ax)
            cbar.set_label('log($\widetilde{\chi}^{2}}$)')

            ax = fig5.add_subplot(spec5[1, 0])
            plt.axis('off')

            ax = fig5.add_subplot(spec5[1, 1])
            # plt.axis('off')
            plt.plot(tempArray,chiSqVariedn,linewidth=3)
            plt.plot([tempMinusErr,bestTemp,tempPlusErr],[chiSqMin+2.3,chiSqMin+2.3,chiSqMin+2.3],linestyle='-',color='k',linewidth=3,marker='|',markersize=10,mew=2)
            # plt.ylim([0,5])
            plt.xlim([min(tempArray),max(tempArray)])
            plt.yscale('log')
            plt.xlabel('T$_{hot}$ (keV)')
            plt.ylabel('$\widetilde{\chi}^{2}}$',rotation=0)
            ax = plt.gca()
            ax.yaxis.set_label_coords(-0.1,0.3)
            figNum += 1
        if(plotting):
            plt.show()

    fittedTemps = np.array(fittedTemps)
    fittedEngys = np.array(fittedEngys)
    tempErrors  = np.array(tempErrors)
    engyErrors  = np.array(engyErrors)
    """ Fitted Temps """
    plt.figure(figNum)
    figNum += 1
    plt.title("Fitted Temps")
    plt.axvline(x=T_base,zorder=14,color="tab:grey",linestyle="--")
    plt.hist(fittedTemps,color="tab:blue")
    plt.xlabel("Temperature (keV)")
    plt.ylabel("Frequency")
    """ Fitted Energies """
    plt.figure(figNum)
    figNum += 1
    plt.title("Fitted Energies")
    plt.axvline(x=A_base,zorder=14,color="tab:grey",linestyle="--")
    plt.hist(fittedEngys,color="tab:red")
    plt.xlabel("Energy (J)")
    plt.ylabel("Frequency")
    """ Temp Errors """
    plt.figure(figNum)
    figNum += 1
    plt.title("Temp Errors")
    plt.hist(100.*tempErrors/fittedTemps,color="tab:green")
    plt.xlabel("Temperature Error (%)")
    plt.ylabel("Frequency")
    """ Energy Errors """
    plt.figure(figNum)
    figNum += 1
    plt.title("Energy Errors")
    plt.hist(100.*engyErrors/fittedEngys,color="tab:purple")
    plt.xlabel("Energy Errors (%)")
    plt.ylabel("Frequency")
    plt.show()
