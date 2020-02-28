"""
NNC: Near Neighbor Clusterization
"""

from array import array
import numpy as np
import ROOT
import root_numpy
ROOT.PyConfig.IgnoreCommandLineOptions = True

def NNCLabel(C,infoCloseT,X,rescale):
    matrix = np.zeros((rescale,rescale))
    hImageBinClu           = matrix.copy()
    aux=0

    for j in range(0,len(infoCloseT)):#range(1,np.shape(clu)[0]):
        aux  = aux+1
        xbox = C[infoCloseT[j]][:,2].astype(int)
        ybox = C[infoCloseT[j]][:,3].astype(int)

        hImageBinClu[xbox,ybox] = aux
        
    labels_ = hImageBinClu[X[:,1],X[:,0]]
    
    tag_= 2*np.ones((np.shape(labels_)[0]),dtype=int)
    return labels_,tag_

def NNClustering(points, thC):
# poi: vettore dei punti X,Y nel piano
    
    import numpy as np
    C = np.zeros((len(points), 4))
    for i in range(0,  len(points)):
        C[i]=[i, 0, points[i,1], points[i,0]]
    
    NofC  = 0
    NeC   = 0
    nloop = 0
    while True:
        i = 1
        nordered = 0
        while (i < len(C)):
            j = 0
            while (j < i):
                sBreak = False
                if abs(C[j,2]-C[i,2])<thC and abs(C[j,3]-C[i,3])<thC:   # close point i to j
                    NofCj = C[j,0]
                    NeCj  = (C[C[:,0]==NofCj][:,1]).max()
                    NofCi = C[i,0]
                    NeCi  = (C[C[:,0]==NofCi][:,1]).max()
                    
                    if NofCi != NofCj:
                        if NeCi == 0:
                            C[i,0]= NofCj
                            C[i,1]= NeCj+1
                        else:
                            if NofCi>NofCj:
                                Ci = np.where(C[:,0]==NofCi)
                                for iCi in range(0, len(Ci)):
                                    C[Ci[iCi], 0] = NofCj
                                    C[Ci[iCi], 1] = NeCj+iCi+1
                            else: 
                                Cj = np.where(C[:,0]==NofCj)[0]
                                for iCj in range(0, len(Cj)):
                                    C[Cj[iCj], 0] = NofCi
                                    C[Cj[iCj], 1] = NeCi+iCj+1
                        sBreak = True
                        nordered += 1
                        break
                j += 1
            i +=1
        if nordered == 0:
            break
        nloop += 1

    sorted_C = np.array(sorted(C, key=lambda x:x[0]))
    return sorted_C

def ClusteringInfo(C):
# ruturn NNClustering clusetr Info array, number of not clusterd, and size of lagre cluster 
    import numpy as np
    NC0  = 0
    NCL  = 0
    maxc = 0
    imax = 0
    info = []
    for i in range(0, len(C)):
        if C[i][1]>0:
            if C[i,1]==1:
                csize = np.where(C[:,0]==C[i,0])[0]
                if len(csize) > maxc:
                    maxc = len(csize)
                    imax = NCL
                info.append(csize)
                NCL +=1
            else:
                NC0 +=1 
    return maxc, imax, info 

class NNC:
    
    def __init__(self, rescale):
        self.rescale = rescale
        self.thC = 2
        
    def fit(self, X):
        
        C = NNClustering(X, self.thC)
        dCloseT, iCloseT, infoCloseT = ClusteringInfo(C)
        clust = NNCLabel(C, infoCloseT, X, self.rescale)
        
        self.labels_, self.tag_ = clust
        
        return self