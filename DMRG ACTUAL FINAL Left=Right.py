# -*- coding: utf-8 -*-
"""
Created on Tue May 31 23:03:06 2022

@author: SS
"""

import time
from DMRG_functions import *

# legend and data containers
lgd, energies, timeTaken, corrFuncValues, lst, entEntropy, sweepEnergy = initContainers()

for m in range(len(stateSize)):  # iterating over how many states to keep
    
    startTime = time.time()
    truncError = 0

    # initial Operators for 2 site blocks                                                           
    HL, HR, SzL, SzR, SuL, SuR, corrT, corrC = initialHalfChainOperators()
    #one cycle of DMRG
    HL,SzL,SuL,corrT,truncError = CycleOnce(0, HL, HR, SzL, SzR, SuL, SuR, corrT, corrC, m, energies, corrFuncValues, truncError, lst)
    # ______________________________________________________________________________

    for i in range(1, k):

        HL, HR, SzL, SzR, SuL, SuR, corrT, corrC = addNewSite(HL, SzL, SuL, corrT)
        
        #for debugging
        if i%10 == 0:
            print("reached:",i)
        
        HL,SzL,SuL,corrT,truncError = CycleOnce(i, HL, HR, SzL, SzR, SuL, SuR, corrT, corrC, m, energies, corrFuncValues, truncError, lst)
    # ______________________________________________________________________________
    
    getEntEntropy(lst, entEntropy, m) #get entanglement entropy
    #finiteCorrections(lst, sweepNum, sweepEnergy, m) #get finite energy corrections
    timeTaken[m] = time.time()-startTime
    print("Truncation Error for", stateSize[m], "states is: ", truncError)
    print("Energy converges to:", np.mean(energies[:,m]), "+-", np.std(energies[:,m]))
    #energies = np.append(energies,sweepEnergy,axis=0) #with finite algorithm energy corrections

printGraphs(lgd, energies, entEntropy, timeTaken, corrFuncValues)
np.savetxt('corr.csv', corrFuncValues, delimiter=",")
np.savetxt('entEntropy.csv', entEntropy, delimiter=",")
np.savetxt('energies.csv', energies, delimiter=",")
np.savetxt('sweepenergies.csv', sweepEnergy, delimiter=",")
