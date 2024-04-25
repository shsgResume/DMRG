# functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix

# constants
sigz = np.array([[0.5, 0], [0, -0.5]]).astype(np.single) #pauli spin z 
sigP = np.array([[0, 1], [0, 0]]).astype(np.single) #pauli spin up

# variables
J = 1.0
k = 300 #length of the half chain in DMRG
LancStates = 6 #number of states to keep in the Lanczos algorithm
sweepNum = 2 #number of sweeps in finite system DMRG
stateSize = [10] #array of the number of states to keep

def initContainers():
    # constants
    lgd = np.zeros(len(stateSize)).astype(np.single)
    for j in range(len(lgd)):
        lgd[j] = str(stateSize[j])
    energies = np.zeros((k, len(stateSize))).astype(np.single)
    timeTaken = np.zeros(len(stateSize)).astype(np.single)
    corrFuncValues = np.zeros((k, len(stateSize))).astype(np.single)
    sweepEnergy = np.zeros((sweepNum*(k-2), len(stateSize)))
    #for entanglement entropy
    lst = []
    entEntropy = np.zeros((k, len(stateSize))).astype(np.single)
    
    return lgd, energies, timeTaken,corrFuncValues, lst, entEntropy, sweepEnergy


def BasisNewNP(eigvectL, eigvalL, n=2):  # default keep 2 lowest states

    if n > len(eigvalL):
        n = len(eigvalL)

    # for left block
    temp = eigvalL
    tempv = eigvectL
    numiter = 0
    basisEigenvaluesL = np.zeros(n)
    basisEigenvectorsL = np.zeros((len(eigvectL), n))

    while numiter < n:
        currentIndex = np.argmax(temp)
        basisEigenvaluesL[numiter] = temp[currentIndex]
        basisEigenvectorsL[:, numiter] = tempv[:, currentIndex]
        # delete the highest value eigenvalue in this iteration
        temp = np.delete(temp, currentIndex, 0)
        # delete the associated eigenstate in this iteration
        tempv = np.delete(tempv, currentIndex, 1)
        numiter += 1

    return basisEigenvaluesL, basisEigenvectorsL


def ReducerNP(rho, n1, n2):

    #rho = rho.reshape([n1, -1], order="C")
    #ra = np.dot(rho, np.conjugate(np.transpose(rho)))
    #rb = ra
    
    rho_tensor=rho.reshape([n1, n2, n1, n2])
    ra = np.trace(rho_tensor, axis1=1, axis2=3) #rho_a
    rb = np.trace(rho_tensor, axis1=0, axis2=2) #rho_b

    return ra, rb


def entanglementEntropy(ev):

    return -np.inner(ev, np.log(ev)) 
    #log2 commonly used in information theory and just provides a different  scale factor


def truncateOperator(Op1, Op2, Op3, Op4, basis):

    O1 = np.conjugate(basis.T) @ (Op1 @ basis)
    O2 = np.conjugate(basis.T) @ (Op2 @ basis)
    O3 = np.conjugate(basis.T) @ (Op3 @ basis)
    O4 = np.conjugate(basis.T) @ (Op4 @ basis)

    return O1, O2, O3, O4


def superBlock(HL, HR, SzL, SzR, SuL, SuR):

    #HL, HR, SzL, SzR, SuL, SuR = csr_matrix(HL),csr_matrix(HR),csr_matrix(SzL),csr_matrix(SzR),csr_matrix(SuL),csr_matrix(SuR)
    H = np.kron(HL, np.identity(len(HR))) + np.kron(np.identity(len(HL)), HR)
    H = H + np.kron(SzL, SzR) + (1/2)*(np.kron(SuL, np.conjugate(SuR.T)))
    H = H + (1/2)*(np.kron(np.conjugate(SuL.T), SuR))
    return csr_matrix(H)


def initialHalfChainOperators():
    # initial is for 2 sites
    HL = np.kron(sigz, sigz)+(1/2)*np.kron(sigP, np.conjugate(sigP.T))+(1/2)*np.kron(np.conjugate(sigP.T), sigP).astype(np.single)
    HR = HL
    SzL = np.kron(np.identity(2), sigz).astype(np.single)
    SzR = SzL
    SuL = np.kron(np.identity(2), sigP).astype(np.single)
    SuR = SuL

    # correlation operator, one as truncated operator and one for calculation in 2 site system
    # initially start with 4 site system
    corrTrunc = np.kron(sigz, np.identity(2))
    corrCalc = np.kron(sigz, sigz)

    return HL, HR, SzL, SzR, SuL, SuR, corrTrunc, corrCalc


def addNewSite(HL, SzL, SuL, corrTrunc):
    # just one new site so the second block is one site operators
    HL = np.kron(HL, np.identity(2)) + np.kron(SzL, sigz) + (1/2)*np.kron(SuL,np.conjugate(sigP.T)) + (1/2)*np.kron(np.conjugate(SuL.T), sigP)
    SzL = np.kron(np.identity(len(SzL)), sigz)
    SuL = np.kron(np.identity(len(SuL)), sigP)

    corrCalc = np.kron(corrTrunc, sigz)
    corrTrunc = np.kron(corrTrunc, np.identity(2))

    HR = HL
    SzR = SzL
    SuR = SuL

    return HL, HR, SzL, SzR, SuL, SuR, corrTrunc, corrCalc


def plotGr(values, lgd, xlabel, ylabel, title):
    plt.plot(values)
    plt.title(title)
    plt.legend(lgd, loc="lower right")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


def keepOperators(lst,H,Sz,Su):
    
    #keeping each operator at each step in order for the ent entropy function
    lst.append({"H":H,
                "Sz":Sz,
                "Su":Su})


def getEntEntropy(lst,storageArray,stateIndex):
    #call at the end
    #get the entanglement entropy as a function of the system size
    startIndex = stateIndex*k
    endIndex = startIndex + k
    lst = lst[startIndex:endIndex]
    L = len(lst)
    g=0
    
    while g<L/2:
        
        HR,SzR,SuR = lst[L-1-g]["H"], lst[L-1-g]["Sz"], lst[L-1-g]["Su"]
        HL,SzL,SuL = lst[g]["H"], lst[g]["Sz"], lst[g]["Su"]
        Hoverall = superBlock(HL,HR,SzL,SzR,SuL,SuR)
        
        eigval, eigvect = eigsh(Hoverall, LancStates)
        gs = eigvect[:, 0]
        Dmat = np.outer(np.conjugate(gs), eigvect[:, 0])
        
        #get the reduced density matrix
        ra, rb = ReducerNP(Dmat, len(lst[g]["H"]), len(lst[L-g-1]["H"]))
        eivaA, eiveA = eigsh(ra,stateSize[stateIndex])
        
        #ent entropy
        storageArray[g,stateIndex] = entanglementEntropy(eivaA)
        g+=2#to avoid even odd oscillation 


def CycleOnce(iter, HL, HR, SzL, SzR, SuL, SuR, corrT, corrC, stateIndex, energyArr, corrArr, truncVal, OpList):
    
    # new superblock Hamiltonian
    H = superBlock(HL, HR, SzL, SzR, SuL, SuR)

    # new densitymat for ground state
    eigval, eigvect = eigsh(H, LancStates)
    energyArr[iter, stateIndex] = np.amin(eigval) / ((iter+2)*2)
    gs = eigvect[:, 0]
    Dmat = np.outer(np.conjugate(gs), eigvect[:, 0])

    # initial correlation function values
    corrArr[iter, stateIndex] = np.inner(np.conjugate(gs), np.matmul(np.kron(corrC, np.identity(len(corrC))), gs))
    # np.kron with identity to give the same dimensions as the ground state; combination of left and right block

    # new reduced density mat
    ra, rb = ReducerNP(Dmat, len(HL), len(HR))

    # entanglement spectrum
    eivaA, eiveA = eigsh(ra,stateSize[stateIndex])

    # new basis
    # keep 4 states in this iteration, max is 4 since thats the dimensions of the reduced density mat
    basisvalL, basisvectL = BasisNewNP(eiveA, eivaA, stateSize[stateIndex])
    truncVal = truncVal+(1-np.sum(basisvalL))

    # truncation and new operators
    SzL, SuL, HL, corrT = truncateOperator(SzL, SuL, HL, corrT, basisvectL)
    keepOperators(OpList, HL, SzL, SuL) #Storing operators in memory
    
    return HL,SzL,SuL,corrT,truncVal


def addNewSiteAlt(direction,H,Sz,Su):
    #for implementation of finite algorithm if needed
    if direction == 'r':#for right block add new site
        H = np.kron(np.identity(2),H)+np.kron(sigz,Sz)+(1/2)*np.kron(Su,np.conjugate(sigP.T)) + (1/2)*np.kron(np.conjugate(Su.T), sigP)
        Sz = np.kron(sigz,np.identity(len(Sz)))
        Su = np.kron(sigP,np.identity(len(Su))) 
    elif direction == 'l':
        H = np.kron(H, np.identity(2)) + np.kron(Sz, sigz) + (1/2)*np.kron(Su,np.conjugate(sigP.T)) + (1/2)*np.kron(np.conjugate(Su.T), sigP)
        Sz = np.kron(np.identity(len(Sz)), sigz)
        Su = np.kron(np.identity(len(Su)), sigP)
    else:
        print("Error in adding new site. Not a valid direction")
    
    return H, Sz, Su


def truncateOperatorAlt(Op1, Op2, Op3, basis):

    O1 = np.conjugate(basis.T) @ (Op1 @ basis)
    O2 = np.conjugate(basis.T) @ (Op2 @ basis)
    O3 = np.conjugate(basis.T) @ (Op3 @ basis)

    return O1, O2, O3


def finiteCorrections(lst, sweepNum, storageArray, stateIndex):
    startIndex = stateIndex*k
    endIndex = startIndex + k
    lst = lst[startIndex:endIndex]
    
    L = len(lst)
    sweep = 0
        
    while sweep<sweepNum:
        
        HL,SzL,SuL = lst[0]["H"], lst[0]["Sz"], lst[0]["Su"] #initial left block in finite algorithm
        #sweep left to right only
        for g in range(0,L-2):
            HR,SzR,SuR = lst[L-3-g]["H"], lst[L-3-g]["Sz"], lst[L-3-g]["Su"]
            HR,SzR,SuR = addNewSiteAlt('r', HR, SzR, SuR)
            
            HL,SzL,SuL = addNewSiteAlt('l', HL, SzL, SuL)
            Hoverall = superBlock(HL,HR,SzL,SzR,SuL,SuR)

            eigval, eigvect = eigsh(Hoverall, LancStates)
            storageArray[sweep*(L-2) + g, stateIndex] = np.amin(eigval) / k
            
            gs = eigvect[:, 0]
            Dmat = np.outer(np.conjugate(gs), eigvect[:, 0])
            ra, rb = ReducerNP(Dmat, len(HL), len(HR))
            eivaA, eiveA = eigsh(ra,stateSize[stateIndex])

            basisvalL, basisvectL = BasisNewNP(eiveA, eivaA, stateSize[stateIndex])
            SzL, SuL, HL = truncateOperatorAlt(SzL, SuL, HL, basisvectL)
            lst[g] = {"H":HL, "Sz":SzL, "Su":SuL} #overwrite the operators kept
        sweep += 1


def printGraphs(lgd, energies, entEntropy, timeTaken, corrFuncValues):
    
    #to rearrange the values in entanglement entropy array. this helps plotting for multiple basis states kept
    for i in range(len(lgd)):
        j=0
        m=0
        while j<(k): 
            entEntropy[m, i] = entEntropy[j,i]
            m+=1
            j+=2

    endIndex = int(k/4)
    entEntropy = entEntropy[0:endIndex,:]
    
    plotGr(energies, lgd, "r'", "Energy of Ground State per site","Energy of the ground state per site with increasing number of sites")
    plotGr(entEntropy, lgd, "r''", "Entropy of Left Block","Entanglement Entropy of the Left Block")
    plotGr(timeTaken, lgd, " ", "Time taken in seconds","Time taken for the iteration")
    plotGr(corrFuncValues, lgd, "r", "Value of correlation Function","Correlation function of the first site with the rest")

