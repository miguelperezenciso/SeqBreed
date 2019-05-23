"""
  SeqBreed.py allows simulating quantitative traits of an arbitrary complexity
    in either diploid or polyploid organisms
    using complete sequence or any number of SNP arrays
  Selection is implemented by using functions in selection.py
  Authors: Miguel Perez-Enciso (miguel.perez@uab.es) & Laura M Zingaretti
"""

import numpy as np
import pandas as pd
from cython.parallel import prange
import re
import os
import io
import sys
import datetime
import gzip
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


############
def read_gen(genFile, ploidy, minMaf=0):
############
    """ This function parses a gen file, either plain text or gz compressed
         rm SNPs with maf < minMaf
        INPUT:
        - genFile(str): input vcf file name with (sequence) genotypes [req]
        - minMaf(float): minimum MAF for a SNP to be retained [0.]
        - ploidy
        OUTPUT
        - ploidy
        - nbase: no. of founders
        - gt: np array with genotypes
        - gf: dataframe with chr and position
        - f: allele frequency of alternative allele
    """
    if ploidy%2 != 0: sys.exit('STOP: ploidy must be even')

    Data = pd.read_csv(genFile, sep='\s+', header=None)

    # genotypes
    gt = Data[Data.columns[2:]]

    # rm genotypes with maf < minMaf
    f = gt.mean(axis=1)
    a, = np.where(np.logical_or(f < minMaf, 1. - f < minMaf))
    gt = np.array(gt.drop(gt.index[list(a)]))
    Data = Data.drop(Data.index[list(a)])
    f = np.delete(np.array(f), list(a))

    # SNP positions
    gf = Data[[0,1]]

    # f converted into maf
    f[np.where(f > 0.5)] = 1. - f[np.where(f > 0.5)]

    # n founders
    nbase = gt.shape[1]//ploidy
    print('ploidy:', ploidy)
    print('N base:', nbase)
    print('N snps:', gt.shape[0])

    return (ploidy, nbase, gt, gf, f)


############
def read_vcf(vcfFile, minMaf=0):
############
    """ This function parses vcf file, either plain text or gz compressed
        remove multiallelic snps and indels
        rm SNPs with maf < minMaf
        INPUT:
        - vcfFile(str): input vcf file name with (sequence) genotypes [req]
        - minMaf(float): minimum MAF for a SNP to be retained [0.]
        OUTPUT
        - ploidy
        - nbase: no. of founders
        - gt: np array with genotypes
        - gf: dataframe with chr and position
        - f: allele frequency of alternative allele
    """

    if vcfFile.endswith('.gz'):
        with gzip.open(vcfFile) as f:
            lines = [l.decode() for l in f if not l.decode().startswith('##')]
    else:
        with open(vcfFile) as f:
            lines = [l for l in f if not l.startswith('##')]

    m = lines[1].split('\t')[9]
    ploidy = 1 + m.count('|') + m.count('/')
    if ploidy % 2 != 0: print('ULL: ploidy must be even')

    # read data using pandas
    Data = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

    # remove indels and multiallelic SNPs
    i, = np.where(np.logical_or(Data['ALT'].apply(len) > 1, Data['REF'].apply(len) > 1))
    Data = Data.drop(Data.index[list(i)])
    # SNP positions
    gf = Data[['#CHROM', 'POS']]
    # genotypes without indels
    gt = Data[Data.columns[9:]]
    # retrieve genotype
    replace = lambda x: pd.Series([i for i in (x.str.split(":").str.get(0))])
    gt = gt.apply(replace)

    # split alleles, delimited by either / or |
    z = []
    for i in range(gt.shape[0]):
        b = [int(j) for j in np.array([x.split("|") for x in [y.replace('/', '|') for y in gt.iloc[i].values]]).ravel()]
        z.append(b)
    gt = pd.DataFrame(z)

    # rm genotypes with maf < minMaf
    f = gt.mean(axis=1)
    a, = np.where(np.logical_or(f < minMaf, 1. - f < minMaf))
    gt = np.array(gt.drop(gt.index[list(a)]))
    Data = Data.drop(Data.index[list(a)])
    gf = Data[['#CHROM', 'POS']]
    f = np.delete(np.array(f), list(a))
    # f converted into maf
    f[np.where(f > 0.5)] = 1. - f[np.where(f > 0.5)]

    # n founders
    nbase = gt.shape[1]//ploidy
    print('ploidy:', ploidy)
    print('N base:', nbase)
    print('N snps:', gt.shape[0])

    return (ploidy, nbase, gt, gf, f)

#############
def print_vcf(vcfFile, inds, genome, gfounders, chip, minMaf=1e-6):
#############
    """ writes vcf file of selected genotypes in chip of inds subset
        vcf is compressed if name ends with 'gz'
        skips if maf < 1e-6 by default
        - vcfFile(str): out vcf file name [req]
        - inds(individuals list): list of class Individual to be printed [req]
        - gfounders(GFounder object): contains founder genotypes [req]
        - chip(Chip object): contains list of snps to be printed [req]
        - minMaf(float): minimum MAF for a SNP to be printed [1e-6]
        WARNING: sloooow
    """
    header  = '##fileformat=VCFv4.3\n'  \
            + '##source=SeqBreed.py\n' \
            + '##reference=chipfile:' + str(chip.name) + '\n' \
            + '##fileDate=' + str(datetime.datetime.now()) + '\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
    header += '\t'.join(str(ind.id) for ind in inds) + '\n'
    vline = '.\tA\tT\t.\t.\t.\t.\t'
    nh   = genome.ploidy
    nind = len(inds)
    ichunks = np.zeros(len(inds * nh),int)
    origins = np.zeros(len(inds * nh),int)  # origin of first chunks
    isnps   = np.zeros(len(inds * nh),int)  # origin of first chunks

    f_open = gzip.open if vcfFile.endswith('.gz') else open
    with f_open(vcfFile, 'wb') as f:
        f.write(header.encode())
        for ichr in range(genome.nchr):
            if len(chip.ipos[ichr])>0:
                cname = genome.chrs[ichr].name
                for ipos in chip.ipos[ichr]:
                    # initializes
                    if ipos == chip.ipos[ichr][0]:
                        for i in range(nind):
                            for h in range(nh):
                                ichunks[nh*i+h] = inds[i].gchunk[ichr][h].isnp.searchsorted(ipos)
                                origins[nh*i+h] = inds[i].gchunk[ichr][h].origin[ichunks[nh*i+h]]
                                isnps[nh*i+h]   = inds[i].gchunk[ichr][h].isnp[ichunks[nh*i+h]]
                    # updates
                    else:
                        # checkkk < or <= AND ALSO NOTE ISNP starts at 0
                        who = np.array(np.nonzero(isnps[:]<ipos),int).flatten()
                        # translate into lambda function
                        if(len(who)>0):
                            for iwho in who:
                                i = iwho//nh
                                h = iwho%nh
                                ichunks[iwho] = inds[i].gchunk[ichr][h].isnp.searchsorted(ipos)
                                origins[iwho] = inds[i].gchunk[ichr][h].origin[ichunks[iwho]]
                                isnps[iwho]   = inds[i].gchunk[ichr][h].isnp[ichunks[iwho]]

                    # vector with genotypes
                    ii = gfounders.cumpos[ichr] + ipos
                    g  = gfounders.g[ii,origins].reshape(nind,nh)
                    maf = g.mean()
                    # skips SNPs with maf < minMaf
                    if min(maf, 1.-maf) > minMaf:
                        vg = ''
                        for ind in g:
                            vg += '|'.join([str(int(elem)) for elem in ind]) + '\t'
                        vg = vg[:-1] # rm last tab
                        line = str(cname) + '\t' + str(genome.chrs[ichr].pos[ipos]) + '\t' + vline + vg + '\n'
                        f.write(line.encode())


########
def do_X(inds, genome, gfounders, chip, minMaf=1e-6):
########
    """ completely rewritten, parallelized with cython
        loop by indiv, chr & hap
        returns genotype matrix
        skips if maf < 0.01 by default
        - inds(individuals list): list of class Individual to be printed [req]
        - genome(GFeatures object) contains genome data [req]
        - gfounders(GFounder object): contains founder genotypes [req]
        - chip(Chip object): contains list of snps to be printed [req]
        - minMaf(float): minimum MAF for a SNP to be returned [1e-6]
    """
    nh = genome.ploidy
    nind = len(inds)
    X = np.array([])

    # filter snps by maf / chr
    ipos_flt = []
    cum_ipos_flt = []
    for ichr in range(genome.nchr):
        ipos = chip.ipos[ichr] # chip snp indices
        if len(ipos)>0:
            cum_ipos = genome.cumpos[ichr] + ipos  # indices in gfounders.g
            ii = np.where(gfounders.f[cum_ipos] > minMaf)  # filtered chip indices
            ipos_flt.append(ipos[ii])
            cum_ipos_flt.append(ipos[ii] + genome.cumpos[ichr])
        else:
            ipos_flt.append([])
            cum_ipos_flt.append([])


    for i in prange(len(inds)):
        # print(i)
        for ichr in prange(genome.nchr):
            if len(ipos_flt[ichr]) > 0:
                g = np.zeros(len(ipos_flt[ichr]))
                for h in prange(nh):
                    ichunks = inds[i].gchunk[ichr][h].isnp.searchsorted(ipos_flt[ichr])  # vector with chunk ids of every ipos
                    # ichunks[nh*i+h] = bisect.bisect_left(inds[i].gchunk[ichr][h].isnp, ipos_flt[ichr])
                    origins = inds[i].gchunk[ichr][h].origin[ichunks]  # origins of every snp for that ind
                    g += gfounders.g[cum_ipos_flt[ichr], origins]  # genotype value
                X = np.append(X, g)

    nsnp = len(X)//nind
    X = X.reshape(nind,nsnp).T
    return X


############
def k2ped(t):
############
    """ returns pedigree starting with founders t generations ago 
        can be saved and reused with smaller t for efficiency
        but need renumbering!
        used for generating random recombinant individs
    """
    n = 2**t
    n0 = 1 # first indiv of current generation
    id0 = np.zeros(n,int)
    id1 = np.zeros(n,int)
    for i in range(t):
        id0 = np.append(id0, np.arange(n0, n0+n//2))
        id1 = np.append(id1, np.arange(n0+n//2, n0+n))
        n0 += n
        n = n//2
    return np.append(id0,id1).reshape(len(id0),2,order='F')


#################
def random_loaded(n, pos, weights):
#################
    """ returns n points distributed along weights
        weights is a cumsum of xover probs
    """
    rxn  = np.sort(np.random.sample(n))
    ipos = weights.searchsorted(rxn)
    bpos = pos[ipos - 1] + (rxn - weights[ipos - 1]) * (pos[ipos] - pos[ipos - 1])
    # trick because first interval (ipos=0) needs to be treated separately
    bpos[ipos == 0] = rxn[ipos == 0] * pos[0]
    return np.array(bpos, int)

###########
def meiosis(chunks, chrfeature, sex, ploidy=2, autopolyploid=True):
###########
    """ allows for polyploidy """
    if ploidy>2:
        xchunk = []
        # sample a random list of pairs
        pairs = list(range(ploidy))
        # free recombination if autopolyploids
        if autopolyploid:
            pairs = np.random.permutation(pairs)
        for i in list(range(0,ploidy,2)):
            hapchunk = mpeiosis([chunks[pairs[i]], chunks[pairs[i+1]]], chrfeature, sex)
            xchunk.append(hapchunk)
    else:
        xchunk = mpeiosis(chunks, chrfeature, sex)
    return xchunk


############
def mpeiosis(chunks, chrfeature, sex):
############
    """ returns a recombinant chromosome
        the chr structure is a list of origin, and position
        TO BE OPTIMIZED
    """
    MXOVER = 3
    REVERSE=[1,0]
    xchunk = ChunkChromosome()

    nxover = min(np.random.poisson(chrfeature.M[sex]), MXOVER)
    if nxover == 0:
        hap = np.random.randint(0, 2)
        xchunk = chunks[hap]

    else:
        # allocates positions according to a weighted random simulator
        xpos = random_loaded(nxover, chrfeature.mappos, chrfeature.mapsum[:,sex])
        # merges info a temporary file
        xc = np.array([chunks[0].pos, chunks[0].origin, np.zeros(len(chunks[0].pos)) - 1, chunks[0].isnp], dtype=int).T
        xc = np.append(xc, np.array([chunks[1].pos, np.zeros(len(chunks[1].pos)) - 1, chunks[1].origin, chunks[1].isnp], dtype=int).T, axis=0)
        xc = np.append(xc, np.array([xpos, np.zeros(nxover) - 2, np.zeros(nxover) - 2, np.zeros(nxover)], dtype=int).T,
                       axis=0)
        xc = xc[np.argsort(xc[:, 0]), :]
        # push up last origin and rm last row
        xc[-2, 2] = xc[-1, 2]
        xc = xc[:-1, :]
        # fills in broken down segments with same origin as above
        n = xc.shape[0]
        for i in list(reversed(range(n - 1))):
            if (xc[i, 1] == -1): xc[i, 1]   = abs(xc[i + 1, 1])
            if (xc[i, 2] == -1): xc[i, 2]   = abs(xc[i + 1, 2])
            if (xc[i, 1] == -2): xc[i, 1:3] = -abs(xc[i + 1, 1:3])

        # first haplotype
        hap = np.random.randint(0, 2)
        # traverse through xover events
        pos = np.array(xc[0, 0])
        ori = np.array(abs(xc[0, 1 + hap]))
        for i in range(1, n):
            # new xover event
            if (xc[i - 1, 1] < 0): hap = REVERSE[hap]
            pos = np.append(pos, xc[i, 0])
            ori = np.append(ori, abs(xc[i, 1 + hap]))

        # marks redundant chunks (those with consecutive identical origins)
        for i in range(n - 1):
            if ori[i] == abs(ori[i + 1]): ori[i] = -1
        # pack up
        xchunk.pos    = pos[ori >= 0]
        xchunk.origin = ori[ori >= 0]
        xchunk.isnp   = xc[ori >= 0, 3]

        # updates no. of snps in that chunk
        # ULL: use np.where(xchunk.isnp==0)[0] does not work!!
        for i in range(len(xchunk.isnp)):
           if (xchunk.isnp[i]==0): xchunk.isnp[i] = len(chrfeature.pos[chrfeature.pos<=xchunk.pos[i]])

    return xchunk


#####################
class ChunkChromosome:
#####################
    """ This class contains a pair position and origin to efficiently store chromosomes """
    def __init__(self, bpos=[], origin=[], nsnp=[]):
        self.pos    = np.array(bpos,dtype=int)
        self.origin = np.array(origin,dtype=int)
        self.isnp   = np.array(nsnp,dtype=int)

    def random(self, origins, chrfeature):
        """ inits chunk in random pieces """
        map = (chrfeature.mapsum[:,0]+chrfeature.mapsum[:,1])/2
        if len(origins)>1:
            self.pos = random_loaded(len(origins)-1, chrfeature.mappos, map)
            self.pos = np.append(self.pos, chrfeature.length)
        else:
            self.pos = np.array([chrfeature.length], int)
        self.origin = origins
        self.isnp = np.zeros(len(self.pos),int)
        # ascertain nsnp per chunk
        for i in range(len(origins)):
            self.isnp[i] = len(chrfeature.pos[chrfeature.pos<=self.pos[i]])

################
class Chromosome:
################
    """ This class contains basic recombination features and methods
        map file should have fields:
          chr bp_position [cM2Mb_sexAverage OR cM2Mb_males cM2Mb_females ]
        by default, the unassigned regions have cM2Mb=1 rec rate
        - name(str): chr name [req]
        - nsnp(int): no. of snps [0]
        - pos(int array): snp base pair position [empty]
        - length(int): bp length [1e10]
        - cM2Mb(float): recombination rate cM/Mb  [1.]
        - xchr(bool): if X chr [False]
        - ychr(bool): if Y chr [False]
        - mtchr(bool): if mitochondrion [False]
    """
    def __init__(self, name, nsnp=0, pos=[], length=1e10, cM2Mb=1., xchr=False, ychr=False, mtchr=False):
        self.name = name
        self.type = 'Autos'
        self.length = length
        self.X = xchr
        self.Y = ychr
        self.mt = mtchr
        self.cM2Mb = cM2Mb
        M = length * 1.e-8 / cM2Mb
        self.M = [M, M]
        self.mapsum = np.array([1.,1.],).reshape(1,2)
        self.mappos = np.array([length],int)
        self.nsnp = nsnp
        self.pos = pos
        if (xchr + ychr + mtchr) > 1:
            sys.exit('ULL: a chromosome cannnot be both sex and mt')
        # no recombination in mitochondria or Y chr
        if mtchr or ychr:
            self.M = [0., 0.]
            self.type = 'MTchr'
            if ychr: self.type = ' Ychr'
        # no recombination in male sexchr X
        elif xchr:
            self.M = [0., M]
            self.type=' Xchr'

    def set_map(self, mapFile):
        """ read and assigns p of xover by region and sex (optionally) """
        if os.path.isfile(mapFile):
            df = pd.read_csv(mapFile, header=None, comment='#', sep='\s+')
            df = df[(df[0].apply(str) == self.name) &  (df[1] <= self.length)]
            if not df.empty and not self.mt and not self.Y:
                df = df.drop([0],axis=1) # rm chr names
                if len(df.columns)>3:
                    df=df[[1,2,3]]
                elif len(df.columns)<3:
                    df=df[[1,2,2]]
                # convert to np
                map = df.values
                # completes map if necessary np.append
                if max(map[:,0]) < self.length:
                    map = np.append(map,[[self.length,self.cM2Mb,self.cM2Mb]],axis=0)
                # total length in Morgans
                d = map[:,0].copy() #ULL: otherwise assigns by reference
                d[1:] -= d[:(-1)]
                # now map stores cumulative rec rate (divided by total M)
                map[:,1] = (d*map[:,1] / sum(d*map[:,1])).cumsum()
                map[:,2] = (d*map[:,2] / sum(d*map[:,2])).cumsum()
                self.mappos = np.array(map[:,0],int)
                self.mapsum = map[:,1:3]
                # total length in Morgans
                M1 = np.dot(d, map[:,1]) * 1e-8
                M2 = np.dot(d, map[:,2]) * 1e-8
                if self.X: M1=0.
                self.M = [M1,M2]


#############
class Genome:
#############
    """ This class contains genome features (snp pos, ploidy, rec rate, sex / mt chr)
        - snpFile(str): snp pos file [req]
        - mapFile(str): file with recombination rate, either sex average or by sex [None]
        - XChr(str): name of X chr [None]
        - YChr(str): name of Y chr [None]
        - MTChr(str): name of mt chr [None]
        - ploidy(int): ploidy level [2]
        For ploidy>2, autoplyploidy is assumed
     """
    def __init__(self, snpFile, mapFile=None, XChr=None, YChr=None, MTChr=None, ploidy=2, autopolyploid=True):
        self.nchr  = 0
        self.ploidy = ploidy
        self.chrs = []
        self.cumpos = []
        self.dictionar = () # links chr name with order
        self.autopolyploid = autopolyploid

        if os.path.isfile(snpFile):
            df = pd.read_csv(snpFile, usecols=[0,1], header=None, comment='#', sep='\s+')
            df[0] = df[0].apply(str)
            chr_names = df[0].unique()
            self.nchr = len(chr_names)
            # dictionar return chr names order
            self.dictionar = dict(zip(chr_names, range(self.nchr)))
            # inits chrs
            for name in chr_names:
                pos = np.array(df[df[0]==name][1]).astype(int)
                nsnp = len(pos)
                length = pos[-1]
                X = False
                if name == str(XChr): X = True
                Y = False
                if name == str(YChr): Y = True
                mt = False
                if name == str(MTChr): mt = True
                if ploidy > 2 and (X or Y or mt):
                    sys.exit('Ploidy > 2 incompatible with sex or mt chrs')
                chr = Chromosome(name=name, pos=pos, nsnp=nsnp, length=length, xchr=X, ychr=Y, mtchr=mt)
                if mapFile != None:
                    chr.set_map(mapFile)
                self.chrs.append(chr)
            # cumulative no of snps / chr
            nsnp = list(self.chrs[i].nsnp for i in range(self.nchr))
            self.cumpos = np.cumsum(nsnp) - nsnp
            if len(self.chrs) != self.nchr:
                sys.exit('ULL: length of chrs must be nchr!!')

        else:
            sys.exit('ULL: file with snp positions missing')

    def print(self):
        """ prints basic genome features """
        print('GENOME FEATURES (ploidy = '+str(self.ploidy)+')')
        print('Chr  Type   bplength    Nsnps length_m length_f')
        length=0
        nsnp=0
        M=np.array([0.,0.])
        for ichr in range(self.nchr):
            c=self.chrs[ichr]
            print('{}  {}  {} {} {}'.format(c.name, c.type, c.length, c.nsnp, c.M))
            length += c.length
            nsnp   += c.nsnp
            M      += np.array(c.M)
        print('{}  {}  {} {}'.format('All      ', length, nsnp, M))
        print()

###############
class GFounder:
###############
    """ This class simply stores genotypes of founders arranged in nbase x ploidy x snp np.array by chr
        reads genotypes from vcfFile and infers ploidy and nbase, computes allele frequencies
        - snpFile(str): file name with SNP positions, obtained from read_vcf [req]
        - nbase(int): no. of individuals in gfounderFile
        - ploidy(int):
        - g(int): array with snp genotypes
    """
    def __init__(self, vcfFile, snpFile, minMaf=0, ploidy=2):
        self.nbase  = 0
        self.ploidy = ploidy
        self.nsnp = 0
        self.g = np.array([],np.uint8)
        self.f = np.array([])

        # automatically detects vcf file, otherwise gen format assumed
        if vcfFile.endswith(tuple(['vcf.gz', 'vcf'])):
            # reads vcf file and assigns ploidy, nbase, genotypes and snp positions
            self.ploidy, self.nbase, self.g, gf, self.f = read_vcf(vcfFile, minMaf)
        else:
            # reads gen file and assigns ploidy, nbase, genotypes and snp positions
            self.ploidy, self.nbase, self.g, gf, self.f = read_gen(vcfFile, ploidy, minMaf)

        self.nsnp = self.g.shape[0]
        # writes vcf positions into file, stored in gf
        gf.to_csv(snpFile, sep=' ', header=False, index=False)

##########
class Chip:
##########
    """ Contains chip genotyping information
        whole sequence is also a Chip class
        - chipFile(str): name of file containing chip snp positions [req]
        - genome(Genome object)
        - gfounders(GFounder object)
        - nsnp(int): no. of snps in chip to be chosen
        - unif(bool): whether snps uniformly or ranmdoly distributed
        - minMaf(float): min maf for a snp to be chosen
        - name(str): chip name
    """
    def __init__(self, genome, gfounders=None, chipFile=None, nsnp=0, unif=False, minMaf=0, name='chip_'):
        self.name = name
        self.nsnp = [] # list of n snsp for each chromosome
        self.ipos = [] # list of ipos for each chromosome

        # a random set of SNPs in the chip are simulated
        if chipFile is None:
            if nsnp == 0: sys.exit('Either chipFile or nsnp must be specified when defining a chip')
            # candidate posns
            jpos = np.arange(gfounders.nsnp)[gfounders.f > minMaf]
            nsnp = min(nsnp, len(jpos))  # just in case

            if unif:
                self.name += 'unif'
                stride = len(jpos) // nsnp
                jpos = jpos[::stride]   # https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html#basic-slicing-and-indexing
            else:
                self.name += 'random'
                jpos = np.sort(np.random.permutation(jpos)[:nsnp])

            # now get index within chr, cumpos[0]=0
            cumpos = genome.cumpos
            ichr = []
            for ic in range(genome.nchr-1):
                # nsnp in that chr
                nsnp = len(jpos[jpos<=cumpos[ic+1]]) - len(ichr)
                ichr = np.append(ichr, np.repeat(ic,nsnp))
            # fillout last chr
            ichr = np.append(ichr, np.repeat(ic+1, len(jpos)-len(ichr)))

            for ic in range(genome.nchr):
                ipos = np.array([])
                nsnp = len(ichr[ichr==ic])
                if nsnp > 0:
                    ipos = jpos[ichr==ic] - cumpos[ic]
                self.ipos.append(ipos.astype(int))
                self.nsnp.append(int(nsnp))

        # reads snp file and assigns snp order per chromosome
        elif os.path.isfile(chipFile):     # check if file exists
            df = pd.read_csv(chipFile, usecols=[0,1], header=None, comment='#', sep='\s+')
            for chr in genome.chrs:
                nsnp = 0
                ipos = np.array([])
                # check if any chip snp in that chr
                pos = np.array(df[(df[0].apply(str) == chr.name)][1])
                if pos.shape[0] > 0:
                    # now check order in dna sequence
                    ipos = chr.pos.searchsorted(pos)
                    # remove elements not in seqpos (index = last_1)
                    ipos = ipos[ipos<len(chr.pos)]
                    ipos = ipos[pos==chr.pos[ipos]]
                    nsnp = len(ipos)
                # append even if nsnp=0
                self.ipos.append(ipos.astype(int))
                self.nsnp.append(int(nsnp))
            if len(self.ipos) !=  genome.nchr:
                sys.exit('ULL: length of chip.ipos must be nchr!!')
        else:
            sys.exit('ULL: chip file does not exist:',chipFile)


    def print(self, genome):
        """ prints generic Chip features """
        print('CHIP FEATURES ' + self.name)
        print('Chr\tNsnps\tFirst posns ... Last')
        for ichr in range(len(self.nsnp)):
            if self.nsnp[ichr]>0:
                print(genome.chrs[ichr].name + '\t' \
                      + str(self.nsnp[ichr]) + 2*'\t' \
                      + ' '.join(map(str, genome.chrs[ichr].pos[self.ipos[ichr][:5]]))  \
                      + ' ... ' \
                      + str(genome.chrs[ichr].pos[self.ipos[ichr][-1]]) )
        print('All ' + str(sum(self.nsnp)))
        print('\n')
        # b= lambda a,index: np.fromiter((a[x,y] for _,(x,y) in enumerate(zip(index,range(0,len(index))))), int)

##########
class QTNs:
##########
    """ check if chr and positions exist in sequence
        chr ids are treated as str
        values are arranged in lists by chr
        error if qtn in Y chr
        - h2(float array): array containing trait h2, its length defines no. of of traits [req]
        - genome(Genome object)
        - qtnFile(str) name of file containing qtn positions and optionally add & dom effects [None]
          if None, qtn positions and effects are sampled (works only for one trait)
        - se(float array): contains environemntal sd, computed from observed genotypes with get_se function [None]
        - delta(float): dominance in olyploids
        - nqtn(int): no. of qtns, inferred from qtnFile, required if qtnFile undefined [0]
        Internal variables:
        - f(float array): qtn frequencies in founders
        - ipos(list of int arrays): for each chr, qtn order in sequence snps
        - va(float array nqtn x ntrait): add variance per qtn per trait
        - vd(float array nqtn x ntrait): dom variance per qtn per trait
        - ploidy
    """
    def __init__(self, h2, genome, qtnFile=None, se=None, nqtn=0, name=None):
        self.name = name
        self.h2 = h2
        self.ipos = []
        self.add = []
        self.dom = []
        self.nqtn = nqtn
        self.se = se
        self.va = []
        self.vd = []
        self.f = []
        self.ntrait = len(h2)
        self.ploidy = genome.ploidy
        self.delta = 1  # for polyploids

        if(len(h2)<1):
            sys.exit('ULL: h2 must be a vector of length>=1')

        if qtnFile and os.path.isfile(qtnFile):
            df = pd.read_csv(qtnFile, header=None, comment='#', sep='\s+')
            # infers no. traits from qtnFile columns; if none, effects are simulated
            ntrait_obs = (len(df.columns)-2) // 2
            if ntrait_obs == 0:
                print('ULL: no qtn effects defined in qtnFile and to be simulated from gamma')
            elif self.ntrait != ntrait_obs:
                print('ULL: defined and observed ntrait differ, minimum taken')
                self.ntrait = min(self.ntrait, ntrait_obs)
            for ichr in range(genome.nchr):
                dchr = df[df[0].apply(str) == genome.chrs[ichr].name]
                ipos = np.array([],int)
                add = []
                dom = []
                if not dchr.empty:
                    if genome.chrs[ichr].Y:
                        sys.exit('ULL: No QTN allowed in Y chr')
                    seqpos = genome.chrs[ichr].pos
                    temp = dchr.values
                    nqtn = temp.shape[0]
                    # iterates over qtn positions for that chr
                    # check positions are in sequence, rm otherwise
                    for iqtn in range(nqtn):
                        # ipos of qtl
                        isnp = seqpos.searchsorted(temp[iqtn,1])
                        # if qtn present in sequence
                        if isnp<len(seqpos) and temp[iqtn,1] == seqpos[isnp]:
                            ipos = np.append(ipos,isnp)
                            if ntrait_obs > 0:
                                # ULL with indices 2:3 = 2; 2:4 = 2,3
                                # THIS ASSUMES add, dom for each trait
                                add = np.append(add, temp[iqtn,2::2])
                                dom = np.append(dom, temp[iqtn,3::2])
                            else:
                                a = np.random.gamma(shape=0.2,scale=5) * np.random.permutation([-1,1])[0]
                                add = np.append(add,a)
                                dom = np.append(dom,0)

                    # add and dom reshaped to allow for ntrait > 1
                    n = len(ipos)   # no. qtn
                    m = len(add)//n # no. traits
                    add = add.reshape(n,m)
                    dom = dom.reshape(n,m)
                self.ipos.append(ipos)
                self.add.append(add)
                self.dom.append(dom)
                self.nqtn += len(ipos)
            # ntrait set back to 1
            #if self.ntrait==0: self.ntrait=1
            if len(self.ipos) !=  genome.nchr:
                print('ULL: length of qtn.ipos must be nchr!!')

        else:
            if not qtnFile: qtnFile=''
            print('QTN file ' + qtnFile + ' does not exist or not specified!\nQTN effects and positions randomly sampled')
            if self.ntrait > 1: sys.exit('ULL: random init only allowed with ntrait=1')
            if self.nqtn == 0: sys.exit('STOP: nqtn must be specified')
            add = np.random.gamma(shape=0.2, scale=5, size=nqtn)
            # randomly assign signs to effects
            add = add * np.random.permutation(np.repeat([1,-1],nqtn))[:nqtn]
            dom = np.repeat(0,nqtn)
            # reshape for coherence
            add = add.reshape(nqtn, 1)
            dom = dom.reshape(nqtn, 1)
            # array with snps per chr
            nsnp = np.array([genome.chrs[i].nsnp for i in range(genome.nchr)])
            # discard snps in Y chr
            for i in range(genome.nchr):
                if genome.chrs[i].Y: nsnp[i]=0
            nsnp_tot = sum(nsnp)
            cum_pos = nsnp.cumsum()
            # chooses positions at random
            ipos = np.sort(np.random.permutation(range(nsnp_tot))[:nqtn])
            ichr = []
            for iqtn in range(nqtn):
                ic = len(cum_pos[cum_pos <= ipos[iqtn]])
                ichr = np.append(ichr,ic)
                if ic>0: ipos[iqtn] -= cum_pos[ic-1]

            for ic in range(genome.nchr):
                self.ipos.append(ipos[ichr==ic])
                self.add.append(add[ichr == ic,:])
                self.dom.append(dom[ichr == ic,:])

            if len(self.ipos) != genome.nchr:
                print('ULL: length of qtn.ipos must be nchr!!')

    def set_se(self, inds):
        """ sets environmental var given h2 and inds subset ,
            all nbase individuals are used by default
        """
        self.se = np.zeros(self.ntrait)
        for itrait in range(self.ntrait):
            va = np.array([ind.g_add[itrait] for ind in inds]).var()
            self.se[itrait] = np.sqrt(va * (1./self.h2[itrait] - 1.))

    def get_var(self, genome, gfounders):
        """ computes observed va & vd per locus """
        va = np.array([])
        vd = np.array([])
        for ichr in range(len(self.ipos)):
            for i in range(len(self.ipos[ichr])):
                ipos = genome.cumpos[ichr] + self.ipos[ichr][i]
                # f is frequency of allele 1
                f = np.mean(gfounders.g[ipos,:])
                self.f = np.append(self.f, f)
                for itrait in range(self.ntrait):
                    a = self.add[ichr][i,itrait]
                    d = self.dom[ichr][i,itrait]
                    alpha = a + d*(1.-2.*f)
                    va = np.append(va , 2*f*(1.-f) * alpha**2)
                    vd = np.append(vd, (2*f*(1-f)*d)**2)
        self.va = va.reshape(self.nqtn, self.ntrait)
        self.vd = vd.reshape(self.nqtn, self.ntrait)

    def print(self, genome, qtnFile=None):
        """ print QTN allele effects to file so they can be reused if printFull=None
            qtnFile = name of file with qtn add & qtn effects
            printFull=True specifies freqs and
        """
        f = sys.stdout if qtnFile == None else open(qtnFile, 'w')
        f.write('#h2 by trait ' + ' '.join( str(self.h2[t])  for t in range(self.ntrait)) + '\n')
        if self.se is None:
           f.write('#se undefined\n')
        else:
           f.write('#se by trait ' + ' '.join(str(self.se[t]) for t in range(self.ntrait)) + '\n')
           f.write('#CHR POS FRQ ' + ' '.join('ADD_' + str(t) + ' DOM_' + str(t) + ' VA_'+str(t) + ' VD_'+str(t) \
                                               for t in range(self.ntrait)) + '\n')
        iq = -1
        for ichr in range(genome.nchr):
            cname = str(genome.chrs[ichr].name)
            for i in range(len(self.ipos[ichr])):
                iq += 1
                bpos = str(genome.chrs[ichr].pos[self.ipos[ichr][i]])
                add  = self.add[ichr][i,:]
                dom  = self.dom[ichr][i,:]
                va = self.va[iq,:]
                vd = self.vd[iq,:]
                frq = "{0:5.3f}".format(self.f[iq])
                line = cname + ' ' + bpos + ' ' + str(frq) + ' '
                line += ' '.join(str(add[t]) + ' ' + str(dom[t]) + ' ' + str(va[t]) + ' ' + str(vd[t]) \
                                     for t in range(self.ntrait))
                f.write(line+'\n')
        if qtnFile != None: f.close()
        print()

    def plot(self, plotFile=None):
        """ plot qtl variance effects
            for all traits by default (itrait=i)
            - freq vs Va / Vd
            - cum freq vs Va / Vd
            - density Va, Vd
        """
        maf = np.array([min(f,1-f) for f in self.f])
        for trait in range(self.ntrait):
            va = self.va[:,trait]
            vd = self.vd[:,trait]
            plt.subplot(3,2,1)
            plt.plot(maf, va, 'bo')
            plt.title('Add variance trait ' + str(trait))
            plt.xlabel('QTN MAF')
            plt.ylabel('Variance')

            plt.subplot(3, 2, 2)
            plt.plot(maf, vd, 'ro')
            plt.title('Dom variance trait ' + str(trait))
            plt.xlabel('QTN MAF')
            plt.ylabel('Variance')

            # histograms
            plt.subplot(3, 2, 3)
            plt.hist(va)
            plt.title('Add variance trait ' + str(trait))
            plt.xlabel('Var')
            plt.ylabel('Density')

            plt.subplot(3, 2, 4)
            plt.hist(vd, facecolor='r')
            plt.title('Dom variance trait ' + str(trait))
            plt.xlabel('Var')
            plt.ylabel('Density')

            # frequency vs cumulative variance
            index = np.argsort(maf)
            plt.subplot(3,2,5)
            plt.plot(maf[index], np.cumsum(va[index]), 'bo')
            plt.title('Cumulative Add variance trait ' + str(trait))
            plt.xlabel('QTN MAF')
            plt.ylabel('Cum Variance')

            plt.subplot(3,2,6)
            plt.plot(maf[index], np.cumsum(vd[index]), 'ro')
            plt.title('Cumulative Dom variance trait ' + str(trait))
            plt.xlabel('QTN MAF')
            plt.ylabel('Cum Variance')

            if plotFile is not None:
                plotfile = plotFile
                if self.ntrait > 1:
                    plotfile = plotFile[:-3] + str(trait) + '.' + plotFile[-3:]
                plt.savefig(plotfile)
            plt.show()



################
class Individual:
################
    """ This class contains individual information
        and associated methods
        - get_gvalue
        - get_y
        ULL: note that y, add, dom can be vectors (ntrait>1)
    """

    def __init__(self, id, ids=0, idd=0, sex=0, y=None, add=None, dom=None, ebv=None):
        self.id      = id
        self.id_sire = ids
        self.id_dam  = idd
        self.sex     = sex
        self.y       = y
        self.qtn     = []   # contain genotypes of each causal locus along chrs
        self.g_add   = add  # additive value for each trait
        self.g_dom   = dom
        self.ebv     = 'nan'   # estimated breeding value
        self.gchunk  = []

    def set_qtn(self, chip, genome, gfounders):
        """ generate genotype values, chip should be QTNs class """
        self.qtn = self.get_genotype(chip,genome,gfounders,hap=False)

    # should be moved to QTNs class to isolate map genotype phenotype
    def set_gvalues(self, qtns):
        """ generates g_add and g_dom values for each trait """
        if len(self.qtn)==0:
            sys.exit('ULL: qtn list empty; Need to call ind.set_qtn(qtns,genome,gbase)')
        self.g_add = np.zeros(qtns.ntrait)
        self.g_dom = np.zeros(qtns.ntrait)
        nh = qtns.ploidy
        for ichr in range(len(self.qtn)):
            if len(self.qtn[ichr])>0:
                ga = self.qtn[ichr] - nh//2
                gd = self.qtn[ichr]
                gd[gd==nh] = 0
                if nh>2:
                    gd[gd>=qtns.delta] = nh # relevant for polyploids only
                    gd = gd/2 # this is to scale dom effect in polyploids such that heteroz = homoz
                for itrait in range(qtns.ntrait):
                    self.g_add[itrait] += np.dot(ga, qtns.add[ichr][:,itrait])
                    self.g_dom[itrait] += np.dot(gd, qtns.dom[ichr][:,itrait])

    # should be moved to QTNs class to isolate map genotype phenotype
    def set_phenotypes(self, qtns):
        """ generates phenotypes """
        # check if genetic values defined
        if len(self.qtn)==0:
            sys.exit('Undefined error variance; Need to call qtns.set_se(inds)')
        if len(self.g_add) == 0: self.set_gvalues(qtns)
        self.y = np.zeros(qtns.ntrait)
        for itrait in range(qtns.ntrait):
            self.y[itrait] = self.g_add[itrait] + self.g_dom[itrait] + np.random.normal(0,qtns.se[itrait])

    # CHECKKKKKKKKKKKKKKKKKKKK
    def get_genotype(self, chip, genome, gfounders, hap='True'):
        """ generate actual genotypes from snp list (defined chip.ipos and chip.ichr)
            chip can be any class with field .ipos
            returns a list of nchr length, each being an array m x nh with allele codes
            or sum of haplotypes if hap='False'
            if hap='False' the sum along m is returned
        """
        g = []
        nh = gfounders.ploidy
        for ichr in range(len(genome.cumpos)):
            gh = np.array([], dtype=int)
            for h in range(nh):
                nchunk = len(self.gchunk[ichr][h].pos)
                # counter for ipos list (0 is first SNP)
                first_snp = -1
                for ichunk in range(nchunk):
                    last_isnp = self.gchunk[ichr][h].isnp[ichunk]  #### WRONG: should be ipos
                    # contain list of snps in that chunk
                    a = (chip.ipos[ichr]>first_snp) & (chip.ipos[ichr]<=last_isnp)
                    list_snp = chip.ipos[ichr][a] + genome.cumpos[ichr]
                    if sum(a) > 0:
                        iorigin = self.gchunk[ichr][h].origin[ichunk]
                        gh = np.append(gh,gfounders.g[list_snp,iorigin])
                    first_snp = last_isnp
            if len(gh)>0:
                gh = gh.reshape(nh,len(gh)//nh)
                # returns the sum of haplotype values per SNP
                if not hap:
                    gh = gh.sum(axis=0)
            g.append(gh)
        return g

    # generates chr chunks (should this be in meiosis or in def_init??)
    def initChunk(self, genome, origin, hap='All'):
        ''' Initializes individual chrs in chrChunk format
            origin should contain all origins, one per ploidy
        '''
        ploidy = len(origin)
        if ploidy != genome.ploidy: print('ERROR: ploidy and no. origins should match')
        for ichr in range(genome.nchr):
            hap = []
            for h in range(ploidy):
                chunk = ChunkChromosome([genome.chrs[ichr].length],[origin[h]],[genome.chrs[ichr].nsnp])
                hap.append(chunk)
            # ULL: no polyploids and sex/mt chrs
            if genome.chrs[ichr].mt:
                hap[0] = hap[1]
            elif genome.chrs[ichr].X and self.sex == 0:
                hap[0] = hap[1]
            # ULL: Y chr in females are generated but make no sense
            elif genome.chrs[ichr].Y:
                hap[1] = hap[0]
            self.gchunk.append(hap)

    def initRandomChunk(self, ibase, t, genome):
        """ Inits chr chunks assuming a poisson process t generations
            ibase should contain all possible nbase*ploidy origins
        """
        for ichr in range(genome.nchr):
            hap = []
            for h in range(genome.ploidy):
                # nxover is sum of Poissons, halved for X chr, none for Y or mt
                if genome.chrs[ichr].mt or genome.chrs[ichr].Y:
                    nxover = 0
                elif genome.chrs[ichr].X:
                    nxover = np.random.poisson(genome.chrs[ichr].M[h]*t*0.5)
                else:
                    nxover = np.random.poisson(genome.chrs[ichr].M[h]*t)
                origins = np.random.permutation(ibase)[:(nxover+1)]
                # dummy initialization
                chunk = ChunkChromosome([genome.chrs[ichr].length], origins, [genome.chrs[ichr].nsnp])
                chunk.random(origins, genome.chrs[ichr])
                hap.append(chunk)
            # ULL: no polyploids and sex/mt chrs
            if genome.chrs[ichr].mt:
                hap[0] = hap[1]
            elif genome.chrs[ichr].X and self.sex == 0:
                hap[0] = hap[1]
            # ULL: Y chr in females are generated but make no sense
            elif genome.chrs[ichr].Y:
                hap[1] = hap[0]
            self.gchunk.append(hap)

    def mate(self, parents, genome):
        """ generates offspring genome in chrChunk format
            ULL: what to do if one parent unknown
        """
        nh = genome.ploidy
        auto = genome.autopolyploid
        # for each chromosome
        for ichr in range(genome.nchr):
            # mt chr is female's origin w/o recombination
            if genome.chrs[ichr].mt:
                hap = [parents[1].gchunk[ichr][1], parents[1].gchunk[ichr][1]]
            # male offspring and sex chr is female's (sex=1) with recombination
            elif genome.chrs[ichr].X and self.sex==0:
                chunk = meiosis(parents[1].gchunk[ichr], genome.chrs[ichr], 1, nh)
                hap = [chunk, chunk]
            elif genome.chrs[ichr].Y:
                hap = [parents[0].gchunk[ichr][0], parents[0].gchunk[ichr][0]]
            # for X chr and female offs, no recombination in male map
            else:
                hap = []
                for parent in parents:
                    # OK check gchunk[ichr] contains all haplotypes
                    chunk = meiosis(parent.gchunk[ichr], genome.chrs[ichr], parent.sex, nh, auto)
                    hap.append(chunk)

            if(nh>2):
                temp = hap[0]
                for c in hap[1]: temp.append(c)
                hap = temp
            if len(hap) != nh:
                sys.exit('Length of hap should be ploidy level')
            self.gchunk.append(hap)

    def print(self, genome, hapFile=None):
        """ prints genotypes from ind in readable format
            writes hapFile format if present
        """
        f = sys.stdout if hapFile == None else open(hapFile, 'a')
        f.write('Ind Parents sex: ' + str(self.id) + ' ' + str(self.id_sire) + ' ' + str(self.id_dam) + '  ' + str(self.sex) + '\n')
        if self.y is not None:
            f.write('Y  g_add  g_dom: ' \
                + ' '.join(map(str, self.y)) + ' ' \
                + ' '.join(map(str, self.g_add)) + ' ' \
                + ' '.join(map(str, self.g_dom)) + '\n')
        for ichr in range(len(self.gchunk)):
            name = genome.chrs[ichr].name
            f.write(' Chr ' + name + '\n')
            for h in range(len(self.gchunk[ichr])):
                f.write('   Pos {}'.format(self.gchunk[ichr][h].pos) + '\n')
                f.write('   ori {}'.format(self.gchunk[ichr][h].origin) + '\n')
                f.write('   Isnp {}'.format(self.gchunk[ichr][h].isnp) + '\n')
        print('\n')


################
class Population:
################
    """ Contains population attributes eg pedigree selection method
        - genome(Genome object): [req]
        - pedFile(str): file name containing pedigree (not required if only founder individuals)
        - t(int array): generation of each individual, used to restrict selectionable individuals
        - mode(str): way to generate rec individuals, 'random' or 'pedigree' ['random']
        - generation
        - qtns(Qtn object)
        - gfounders(GFounder object)
    """
    def __init__(self, genome, pedFile=None, t=[], label=None, generation=None, qtns=None, gfounders=None):
        self.inds = []
        self.n = 0
        self.t = [] # generation per individual
        self.label = [] # arbitrary label, can be used to provide colors in say PCA
        self.nbase = gfounders.nbase # new variable can be useful ...
        self.nbase1 = gfounders.nbase # includes actual number w/o parents

        # reads pedigree and init individuals
        if pedFile is None or os.path.isfile(pedFile):
            # quick and dirty generate a pedigree with unrelated individuals if pedFile undefined
            if pedFile is None:
                n = self.nbase
                df = pd.DataFrame([list(range(1,n+1)), [0]*n, [0]*n]).T
            else:
                df = pd.read_csv(pedFile, header=None, comment='#', sep='\s+')
            self.nbase1 = sum(df[1]==0)
            nh = genome.ploidy
            # sex in pedigree, convert 2,1 --> 0,1
            if len(df.columns) > 3:
                sex = np.array(df[3],int)
                if max(sex)>1: sex = sex - 1

            i=0
            for id in df[0]:
                # female
                if id in df[2].values:
                    idsex = 1
                elif len(df.columns) > 3: # sex in pedigree
                    idsex = sex[i]
                else:
                    idsex = 0
                # base indiv if sire=0 and dam=0
                sire = df[1][i]
                dam = df[2][i]
                if sire==0 or dam==0:
                    if id<=self.nbase: # actual base individual
                        # ULL: origins start at 0
                        origin = np.array(list(range(2*i,2*i+(nh))))
                        self.addBaseInd(genome, origins=origin, id=id, sex=idsex, qtns=qtns, gfounders=gfounders)
                    else:  # randomly generated base individual
                        self.addRandomInd(genome, self.nbase, id=id, sex=idsex, t=5, mode='random', qtns=qtns,
                                     gfounders=gfounders)
                elif sire > 0 and dam > 0:
                    # ULL: ids start at 1 but ind.id list index start at 0
                    parents = [self.inds[sire-1],self.inds[dam-1]]
                    self.addInd(parents, genome, id=id, sex=idsex, qtns=qtns, gfounders=gfounders)
                else:
                    print(df[:][i])
                    sys.exit('STOP: either both parents known or both unknown')
                i += 1

            # assigns generation number, 0 for nbase, 1 for rest
            if generation is None:
                nbase1 = min(self.nbase1, self.n) # just in case not all base indivs in pedigree
                self.t = np.concatenate((np.repeat(0,nbase1), np.repeat(1,(self.n-nbase1))))
            elif len(generation) != self.n:
                sys.exit('generation must be same length as pedigree')
            else:
                self.t = generation
        else:
            sys.exit('ULL: no pedfile found')

        # determine Var_e if unknown and generates phenotypes
        if qtns:
            if qtns.se is None:
                qtns.set_se(self.inds[:self.nbase])
            for ind in self.inds:
                ind.set_phenotypes(qtns)

    def addBaseInd(self, genome, origins, id=None, sex=None, qtns=None, gfounders=None):
        """ add new individual as offspring of parents
            parents are Individual objects
        """
        if not id: id = len(self.inds)+1
        if not sex: sex = np.random.randint(0,2)

        ind = Individual(id, 0, 0, sex=sex)
        ind.initChunk(genome, origin=origins)

        if qtns:
            ind.set_qtn(qtns, genome, gfounders)  # get qtn genotypes
            ind.set_gvalues(qtns)         # get add and dom values
            if qtns.se: ind.set_phenotypes(qtns)
        self.inds.append(ind)
        self.n += 1
        self.t = np.append(self.t, 0)

    def addRandomInd(self, genome, nbase, id=None, sex=None, t=0, k=5, mode='random', qtns=None, gfounders=None):
        """ adds a new individual by recombining existing genomes
            take 2**k individuals and mate them randomly for k generations
            mode='random' random chunks of chrs are generated and origins randomly assigned
            mode='pedigree' generates random pedigree
        """
        mode = mode.lower()
        if id is None: id = len(self.inds)+1
        if sex is None: sex = np.random.randint(0,2)

        # range of potential parents
        ibase = list(range(nbase))
        if len(ibase) < 2**k:
            ibase = ibase * round((2**k / len(ibase) + 0.5))
        ibase = ibase[:(2**k + 1)]
        ibase = np.random.permutation(ibase)

        # mode for adding randomly generated individuals
        if mode == 'random':
            ind = Individual(id, 0, 0, sex=sex)
            ind.initRandomChunk(ibase, k, genome)

        # via a pedigree with 2**k ancestors
        elif mode == 'pedigree':
            nh = genome.ploidy
            ped = k2ped(k)  # generates pedigree of fake individual
            nfake = ped.shape[0]
            fake_pop = []
            # follow random pedigree
            for i0 in range(nfake):
                # id is always pop n + id
                ind = Individual(id, ids=0, idd=0, sex=sex)
                if ped[i0,1]==0:
                    i = ibase[i0] # random origin
                    ind.initChunk(genome, origin=np.array(list(range(2*i, 2*i + nh))))
                else:
                    ids = ped[i0,0]-1
                    fake_pop[ids].sex=0
                    idd = ped[i0,1]-1
                    fake_pop[idd].sex=1
                    parents = [fake_pop[ids], fake_pop[idd]]
                    ind.mate(parents, genome)
                fake_pop.append(ind)

        if qtns:
            ind.set_qtn(qtns, genome, gfounders)  # get qtn genotypes
            ind.set_gvalues(qtns)         # get add and dom values
            if qtns.se is not None: ind.set_phenotypes(qtns)
        # last individual is finally appended
        self.inds.append(ind)
        self.n += 1
        self.t = np.append(self.t, t)

    def addInd(self, parents, genome, gfounders=None, qtns=None, id=None, sex=None, t=0):
        """ add new individual as offspring of parents
            parents is a list of two Individual objects
        """
        if id is None:
            id = len(self.inds)+1
        if sex is None:
            sex = np.random.randint(0,2)
        ind = Individual(id, parents[0].id, parents[1].id, sex=sex)
        ind.mate(parents, genome)
        if qtns is not None:
            ind.set_qtn(qtns, genome, gfounders)  # get qtn genotypes
            ind.set_gvalues(qtns)         # get add and dom values
            if qtns.se is not None:
                ind.set_phenotypes(qtns)
        self.inds.append(ind)
        self.n += 1
        self.t = np.append(self.t, t)

    def addPed(self, ped, genome, qtns, gfounders, t=None):
        """ ped should be a np array with id, sire, dam, sex """
        n = ped.shape[0]
        # assumes ped makes up a new generation
        if t is None:
            t = max(self.t)+1
        for i in range(n):
            # ULL: index in inds[] is id in ped minus one
            id = ped[i,0]
            papa = ped[i,1]
            mama = ped[i,2]
            if papa==0 or mama==0:
                sys.exit('New individuals added in ped should have sire and dam known')
            if ped.shape[1]==4:
                sex = ped[i,3]
            else:
                sex = np.random.randint(0,2)
            parents = [self.inds[papa-1], self.inds[mama-1]]
            self.addInd(parents, genome, gfounders, qtns, id, sex, t)

    def ped(self):
        """ returns pedigree as an nx3 array """
        ped = np.array([],int)
        for i in range(self.n):
            ped = np.append(ped, [self.inds[i].id, self.inds[i].id_sire, self.inds[i].id_dam])
        ped = ped.reshape(self.n,3)
        return ped

    def summary(self, genome):
        """ prints general information """
        pass

    def ivalues(self):
        """ returns y,add,dom,ebv values
        return as pd frame? """
        n = self.n
        ntrait = len(self.inds[0].y)
        v = np.array([])
        for i in range(n):
            v = np.append(v, list(float(self.inds[i].y[j]) for j in range(ntrait) ))
            v = np.append(v, list(float(self.inds[i].g_add[j]) for j in range(ntrait)))
            v = np.append(v, list(float(self.inds[i].g_add[j]+self.inds[i].g_dom[j]) for j in range(ntrait)))
            v = np.append(v, self.inds[i].ebv)
        return v.reshape(n, 3*ntrait+1).astype('float')

    def print(self, popFile=None):
        ''' prints basic pop info '''
        f = sys.stdout if popFile is None else open(popFile, 'w')
        ntrait = len(self.inds[0].y)
        f.write('#ID FATHER MOTHER SEX GENERATION ' + ' '.join('Y_' + str(t) + ' ADD_' + str(t) + ' GEN_' + str(t) for t in range(ntrait)))
        f.write(' EBV \n')
        i=0
        for ind in self.inds:
            y = ind.y
            add = ind.g_add
            dom = ind.g_add + ind.g_dom
            f.write(str(ind.id) + ' ' + str(ind.id_sire) + ' ' + str(ind.id_dam) + ' ' + str(ind.sex) + ' ' + \
                    str(self.t[i]) + ' ' + \
                    ' '.join(str(y[t]) + ' ' + str(add[t]) + ' ' + str(dom[t]) for t in range(ntrait)) + ' ' +\
                    str(ind.ebv) + '\n')
            i+=1
        if popFile != None: f.close()

    def plot(self, y=None, ebv=None, itrait=0):
        """ basic boxplot """
        # boxplots phenotype by default, remove NA
        if ebv is None:
            y = np.array(list(self.inds[i].y[itrait] for i in range(self.n)))
            t = self.t
            Y = 'Phenotype ' + str(itrait)
        # ebvs can be missing
        else:
            y0 = np.array(list(float(self.inds[i].ebv) for i in range(self.n)))
            y = y0[np.logical_not(np.isnan(y0))]
            t = self.t[np.logical_not(np.isnan(y0))]
            Y = 'EBV'

        # create a data frame
        data = pd.DataFrame({'Generation': t, Y: y})

        ax = sns.boxplot(x='Generation', y=Y, data=data)
        plt.show()
        plt.close('all')
