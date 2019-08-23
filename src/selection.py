""" This module incorporates several methods to seamlessly implement selection
    in seqbreed software
    Main arguments
        - criterion: Random, Phenotype, BLUP, GBLUP, Sstep, ...
        - generation: discrete(nsize) / continuous
        - intensity: male, female
        - mating: assortative / random
    The main output is an extended pedigree with offspring of selected parents

    ok: gwas option, fdr option
    TODO:
    - optimize reading /writing snp data
    - Add qtn positions in gwas plot
    - add labels pca,
    - add xvalidation = cross_val_predict(lr, boston.data, y, cv=10)
    - ULL: check pedigree order in complex settings (BLUP / GSSTEP)
"""

'''
    Copyright (C) 2019  Miguel Perez-Enciso

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import numpy as np
import pandas as pd
import scipy as sp
import gzip
import re
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.decomposition import PCA
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import fdrcorrection as fdr

from genome import *

#Rstats = importr('stats')

# missing code, ULL: it must be float
MISSING = -999999.

#########
class Pca:
#########
    """ plots PCA of selected inds and markers """
    def __init__(self, X):
        self.X = X
        self.p = []

    def fit(self):
        ''' performs PCA on pop and chip snps
            scikit assumes nsamples x nvariables
        '''
        pca = PCA(n_components=2)
        self.p = pca.fit(self.X.T).fit_transform(self.X.T)

    def plot(self, plotFile=None, labels=None):
        plt.title('PCA decomposition')
        if labels is None:
            plt.scatter(self.p[:, 0], self.p[:, 1])
        else:
            for i in range(max(labels)+1):
                plt.scatter(self.p[labels==i, 0], self.p[labels==i, 1], label='Gen '+str(i))
            plt.legend(loc='best')
        if plotFile is None:
            plt.show()
        else:
            plt.savefig(plotFile)
        plt.close('all')

##########
class Gwas:
##########
    """ Implements GWAS
        - X: the genotypes used to perform the gwas
        - b(float array): regression coefficients
        - pvalue(float array)
        - qvalue(float array)
    """
    def __init__(self, X, chip):
        self.X = X
        self.chip = chip
        self.b = []
        self.se = []
        self.pvalue = []
        self.qvalue = []

    def fit(self, inds=None, y=None, vcfFile=None, trait=0):
        """ carry out GWAS
            y can be directly provided or obtained from pop object
            returns [ b, se, pvalue, fdr]
            can do for any trait=0:ntrait, or as provided in y vector
            - X(nind x nmkr array): contains genotypes
            - y(nind array): contains phenotypes
            - trait(int): trait index to be used, overriden if y provided [0]
        """
        # get phenotype if y not provided
        if y is None:
            nind = len(inds)
            y = np.array(list(inds[i].y[trait] for i in range(nind)))
        else:
            nind = len(y)

        '''
        # genotypes from vcfFile (deprecated)
        if X is None and vcfFile is None:
                sys.exit('genotypes must be in vcfFile or passed as X matrix')

        elif vcfFile is not None:
            f_open = gzip.open if vcfFile.endswith('.gz') else open
            sep = '\||/'
            with f_open(vcfFile) as f:
                for line in f:
                    try:
                        line = line.decode()
                    except AttributeError:
                        pass
                    if line.startswith('#'): continue
                    line = line.rstrip().split('\t')
                    genotypes = line[9:]
                    gt = []
                    for g in genotypes:
                        ig = g.split(':')[0]
                        gt.append(re.split(sep, ig))
                    # convert into 1D list and to int
                    gt = np.array(list((map(int, sum(gt, [])))))
                    # trick to join by individual genotypes, and then sum
                    gt = np.split(gt, nind)
                    gt = np.asarray(list(map(sum, gt)))
                    b, intercept, r_value, p_value, std_err = stats.linregress(gt, y)
                    out = np.append(out, [b, std_err, p_value])
        '''

        # output
        out = np.array([])
        # genotypes provided in X
        for i in range(self.X.shape[0]):
              # trick to sum genotypes over haplotypes
              x = np.split(self.X[i, :], nind)
              x = np.asarray(list(map(sum, x)))
              b, intercept, r_value, p_value, std_err = stats.linregress(x, y)
              out = np.append(out, [b, std_err, p_value])

        out = out.reshape(len(out) // 3, 3)
        self.b, self.se, self.pvalue = out[:,0], out[:,1], out[:,2]
        # FDR obtained from statsmodels package
        self.fdr = fdr(self.pvalue)[1]

    def plot(self, plotFile=None, fdr=None):
        """ GWAS plot, one color per chr
        - plotFile(str): file name where plot is printed, can be png , pdf ... [None]
          prints to STDOUT if undefined
        - fdr(bool): plots fdr instead of -log10 p-values
        - if chip provided, each chr by diff color
        """

        nchr = len(self.chip.nsnp)
        col = np.tile(['r','b'], nchr)[:nchr] # r/b color per chr
        col = np.repeat(col, self.chip.nsnp)
        title =  'GWAS ' + self.chip.name


        # p-values are printed
        if fdr is None:
            plt.scatter(list(range(len(self.pvalue))), -np.log10(self.pvalue), c=col, marker='o')
            plt.ylabel('-log10 P-value')
        # FDR values are printed
        else:
            plt.scatter(list(range(len(self.fdr))), -np.log10(self.fdr), c=col, marker='o')
            plt.ylabel('-log10 FDR')

        plt.title(title)
        plt.xlabel('SNP')
        if plotFile is not None:
            plt.savefig(plotFile)
        else:
            plt.show()
        plt.close('all')

    def print(self, genome, gwasFile=None):
        """ prints to file or STDOUT
        - genome(Genome object): [req]
        - chip (chip object): [req]
        - gwastFile(str): file name where gwas results are printed [None]
          prints to STDOUT if undefined
        """
        f = sys.stdout if gwasFile == None else open(qtnFile, 'w')
        f.write('#CHR POS COEFF SE PVALUE FDR' + '\n')
        isnp = -1
        for ichr in range(genome.nchr):
            cname = str(genome.chrs[ichr].name)
            for i in range(len(self.chip.ipos[ichr])):
                isnp += 1
                bpos = str(genome.chrs[ichr].pos[self.chip.ipos[ichr][i]])
                line = cname + ' ' + bpos + ' ' + str(self.b[i]) + ' ' + str(self.se[i]) \
                       + ' ' + str(self.pvalue[i]) + ' ' + str(self.fdr[i])
                f.write(line+'\n')
        pass



##########
def doGRM(X, grmFile=None, nh=2):
##########
    """ Computes GRM
        - X contains genotypes
        - grmFile(str): file name containing GRM [req]
          compressed if file name ends with 'gz'
        - nh: ploidy
        - maf(float): minimum maf for a snp to be considered [1e-6]
    """
    p = X.mean(axis=1)/nh
    c = sum(nh * p * (1.-p))
    s = StandardScaler(with_std=False)
    X = s.fit_transform(X.T).T
    G = np.matmul(X.T, X) / c
    return G


##############
def doAInverse(ped):
##############
    """
    Returns A-inverse using Henderson's rules, w/o considering inbreeding
    ped is a nind x 3 np.array with id, father and mother ids, 0 for unknown parent
    """
    w = np.array([1., -0.5, -0.5])
    res = np.array([2., 4./3., 1., 0.]) # why dim 4????
    nind = ped.shape[0]

    # a dummy 0 row and 0 column is created to facilitate addressing posns
    AI = np.zeros(shape=((nind+1),nind+1))

    c = np.zeros(nind, dtype=int)
    c[ped[:,1]==0] += 1
    c[ped[:,2]==0] += 1

    # ULL with id indices
    for i in range(nind):
        for k1 in range(3):
            for k2 in range(3):
                AI[ped[i,k1],ped[i,k2]] += w[k1] * w[k2] * res[c[i]]

    # rm row 0 and col 0
    return AI[1:,1:]

###########
def dogblup(h2, y, grmFile=None, G=None, invert=True):
###########
    """ (G)BLUP evaluation, returns EBVs
        G can be passed as argument or printed in grmFile
        G is assumed to be the inverse if Invert=False
        - h2(float): h2 used [req]
        - y(array float): phenotypes, can contain missing values as coded by MISSING [req]
        - grmFile(str): file containing Cov matrix of breeding values or its inverse [None]
        - G(nind x nind float array): np.array with Cov matrix of breeding values or its inverse [None]
        - invert(bool): True if G should be inverted or False if already inverted [True]
    """
    # G is read from file
    if grmFile!=None:
        G = np.array(pd.read_csv(grmFile, header=None, comment='#', sep='\s+'))

    nind = len(y)
    if nind != G.shape[0]:
        sys.exit('GRM matrix size must be equal to no. inds')

    # builds MME
    if invert: 
        np.fill_diagonal(G, np.diag(G) * 1.05)
        lhs = np.linalg.inv(G)

    else:
        lhs = G

    y[y==MISSING] = 0
    x1 = np.repeat(1.,nind)
    x1[y==0] = 0
    x2 = np.append(x1, len(x1[x1!=0])) # last element is no. of data

    lhs = lhs * ((1.-h2)/h2)
    np.fill_diagonal(lhs, np.diag(lhs)+x1)
    lhs = np.vstack([lhs, x1])  # add row
    lhs = np.hstack((lhs, x2.reshape(nind+1,1))) # add col
    rhs = np.append(y, np.sum(y))

    # solves, all els but last are the ebvs
    ebv = np.linalg.solve(lhs, rhs)[:-1]
    return ebv


###########
def dosstep(h2, y, im, AInv, A, G=None, grmFile=None):
###########
    """ Single step evaluation
        - h2(float): h2 used [req]
        - y(array float): phenotypes, can contain missing values as coded by MISSING [req]
        - im(int array): set of genotyped individuals (indexed starting with 0) [req]
        - ainvFile(str): file name that contains A inverse [req]
        - aFile(str): file name that contains A [req]
        - grmFile(str): file name that contains GRM of genotyped individuals [req]
    """
    nind = len(y)
    if nind != A.shape[0]: sys.exit('NRM matrix size must be equal to # inds')

    ngind = len(im)
    if grmFile is not None:
       G = np.array(pd.read_csv(grmFile, header=None, comment='#', sep='\s+'))
    np.fill_diagonal(G, np.diag(G) * 1.05)
    if ngind != G.shape[0]: sys.exit('GRM matrix size must be equal to # genotyped inds')

    # builds H-1
    Hinv = AInv
    Ginv = np.linalg.inv(G)
    A22inv = np.linalg.inv(A[im,:][:,im])
    # ULLLL: this did not work!!!
    # Hinv[im,:][:,im] = Hinv[im,:][:,im] + Ginv - A22inv
    Z = Ginv - A22inv
    j=0
    for i in im:
        Hinv[i,im] += Z[j,:]
        j+=1

    ebv = dogblup(h2, y, G=Hinv, invert=False)
    return ebv

#########
def doEbv(pop, criterion='random', h2=None, X=None, grmFile=None, mkrIds=None, yIds=[], nh=2, itrait=0):
#########
    """ returns estimated breeding values (EBVs), assume selection is on first trait (itrait=0)
    - pop (Population object)
    - criterion(str): evaluation method: 'random', 'phenotype', 'blup', 'gblup', 'sstep' ['random']
    - h2(float): used h2, required in blup, gblup or sstep [None]
    - X: contains genotypes
    - grmFile(str): file containing Cov matrix of breeding values [None]
    - mkrIds(int array): set of genotyped individuals (indexed starting with 0) [None, req for sstep]
    - yIds(int array): integer array specifying individuals with data (indexed starting with 0) [all]
    - trait(int): trait index for which evaluation is performed [0]
    Either grmFile or G is required for gblup and sstep. In sstep, marker files should contain
    only information for genotyped individuals and in the same order.
    """
    criterion = criterion.lower()
    nind = len(pop.inds)

    # remove missing phenotypes ( none by default)
    if len(yIds)==0:
        yIds = np.arange(nind)
    y0 = np.array(list(pop.inds[i].y[itrait] for i in range(nind)))
    y = np.repeat(MISSING, nind)
    y[yIds] = y0[yIds]

    # check h2 is provided
    if criterion == 'blup' or criterion == 'gblup' or criterion=='sstep':
        if not h2: sys.exit('ULL: h2 must be specified for blup gblup or sstep')

    # computes A inverse & A required for blup or sstep
    if criterion == 'blup' or criterion=='sstep':
        ped = pop.ped()   #--> returns peigree
        if criterion == 'blup':
            AI = doAInverse(ped) #--> computes A-inv
        elif criterion == 'sstep':
            AI = doAInverse(ped)    # --> computes A-inv
            A =  np.linalg.inv(AI)  #--> computes A

    # computes GRM if undefined, required for glbup or sstep
    if (criterion == 'gblup' or criterion=='sstep'):
        G = doGRM(X, nh)

    # ---> Computes EBVs
    # pure drift
    if criterion == 'random':
        ebv = np.random.random(nind)
    # mass selection
    elif criterion == 'phenotype':
        ebv = y
    # BLUP
    elif criterion == 'blup':
        ebv = dogblup(h2, y, G=AI, invert=False)
    # GBLUP
    elif criterion == 'gblup':
        ebv = dogblup(h2, y, G=G, grmFile=None)
    # single step
    elif criterion == 'sstep':
        ebv = dosstep(h2, y, mkrIds, AInv=AI, A=A, G=G, grmFile=None)
    # unknown
    else:
        sys.exit('Unknown selection criterion ' + criterion)

    # assigns EBV to individuals
    for i in range(nind):
        pop.inds[i].ebv = ebv[i]


################
def ReturnNewPed(pop, nsel, famsize, mating='random', generation=0):
################
    """ Returns a new pedigree given ebv and selection restrictions
    The new pedigree contains the offspring of selected individuals
    - pop(Population object)
    - nsel(int 2 els. array): number of males and females selected
    - famsize(int): number of offspring per female
    - mating(str): 'random' (=='r') or 'assortative' (=='a') [random]
    - generation(int): only individuals from generation onwards are considered as potential parents [all]
    """
    nind = len(pop.inds)
    ebv = np.array(list(pop.inds[i].ebv for i in range(nind)))
    sex = np.array(list(pop.inds[i].sex for i in range(nind)))
    # ULLL: note id-1 is referring to i
    id = np.array(list(pop.inds[i].id for i in range(nind))) - 1
    t = pop.t[id]

    # merge in temporary array [ebv, sex, t, id, index]
    temp = np.stack((ebv, sex, t, id, list(range(nind)))).T

    # select individuals from avail generations (all by default)
    temp = temp[temp[:,2]>=generation,:]

    # ix is indices selected males, note -ebv for decreasing order
    males = temp[temp[:,1]==0,:]
    ixmales = np.argsort(-males[:,0])[:nsel[0]]

    females = temp[temp[:,1]==1,:]
    ixfemales = np.argsort(-females[:,0])[:nsel[1]]

    # if not assortative mating (=random)
    if not mating.lower().startswith('a'):
        ixfemales = np.random.permutation(ixfemales)

    ped = []
    i = 0
    ij = -1
    for i in range(len(ixmales)):
        # NOTE: we need to add +1 to obtain id from index
        papa = int(temp[ixmales[i],3]) + 1
        for j in range(len(ixfemales)//len(ixmales)):
            ij += 1
            mama = int(temp[ixfemales[ij],3]) + 1 # idem
            for k in range(famsize):
                nind += 1
                sex = np.random.randint(0,2)
                ped.append([nind, papa, mama, sex])

    # returns pedigree of next generation
    return np.array(ped, int)
