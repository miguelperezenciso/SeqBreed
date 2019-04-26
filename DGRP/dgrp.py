"""
    Script illustrating genomic selection with thinned DGRP data

"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import copy

# SeqBreed modules
from SeqBreed import genome as gg
from SeqBreed.selection import selection as sel

############ main #################
# current dir
cdir = os.getcwd()

# working directory
wdir = cdir

# input file dir
ddir= cdir

# input files
chipfiles = [ ] # can include sequence as output from read_vcf
mapfile = ddir+'/dgrp.map' # note males do not recombine
vcffile = ddir+'/dgrp.vcf.gz'

# goto working directory
os.chdir(wdir)

# working files
seqfile = 'seq.pos'
gpickfile = 'gfounder.p'

# Reading vcf file can take a while for big files, but need to be done only once
# and gbase object can be saved eg in pickle
if not os.path.isfile(gpickfile):
    gbase = gg.GFounder(vcfFile=vcffile, snpFile=seqfile)
    pickle.dump(gbase, open(gpickfile, "wb"))
else:
    # read founder snps from pickle file
    gbase = pickle.load(open(gpickfile, "rb"))

# generates genome object with chr names, recombination map, etc
# requires snpFile generated in previous step, this file simply stores snp positions
gfeatures = gg.Genome(snpFile=seqfile, mapFile=mapfile, ploidy=gbase.ploidy)
# prints some basic info
gfeatures.print()

# 20 QTNs are randomly positioned along the genome
qtn = gg.QTNs(h2=[.5], genome=gfeatures, nqtn=30, qtnFile=None)
# get add and dom variance per locus, assuming equilibrium
qtn.get_var(gfeatures, gbase)

# print and plot main qtn features
qtn.print(gfeatures)
qtn.plot()

# chipfiles contains list of chips, sequence is an additional chip
chipfiles.append(seqfile)
# add chip information
chips = []
for file in chipfiles:
    chip = gg.Chip(chipFile=file, genome=gfeatures, name=file+'_chip')
    chips.append(chip)

# generate individual genomes through pedigree, it also computes se
# CHANGE TO INIT WITH NBASE > GBASE.NBASE, THEN INIT QTL , SET GVALUES AND Y FOR FOUNDER INDIVIDUALS
pop = gg.Population(gfeatures, pedFile=None, generation=None, qtns=qtn, gfounders=gbase)
print('No. of individuals is '+str(len(pop.inds)))

# SELECTION is implemented in cycles.

# four replicates, base populations is stored in pop0
nrep = 4

# Initializes and saves base populations
pop0 = copy.deepcopy(pop)

# no. of males and females selected, here 5 & 10
nsel = [5, 10] 
# no. of offspring per female, thus 100 individuals per generation
noffspring = 10

# no. of generations of selection
ngen = 3

# PHENOTYPIC SELECTION, discrete generations, random mating
print('\nMass selection, discrete generations, assortative mating')
for irep in range(nrep):
    print('Replicate ' + str(irep))
    # init population
    pop = copy.deepcopy(pop0)
    for t in range(ngen):
        print('  Generation, N: ' + str(t) + ' ' + str(len(pop.inds)))
        # step 1: do evaluation, adds ebv to pop.inds[:].ebv
        sel.doEbv(pop, criterion='phenotype')
        # step 2: pedigree with offspring of selected individuals, discrete generations
        newPed = sel.ReturnNewPed(pop, nsel, famsize=noffspring, mating='assortative',  generation=t)
        # step 3: generates new offspring (this function adds QTN genotypes, true bvs and y)
        pop.addPed(newPed, gfeatures, qtn, gbase)
print('No. of individuals is ' + str(len(pop.inds)))

# GBLUP USING SEQUENCE chip[0], SNPs MAF > 0.05, continuous generation
print('\n'+'GBLUP selection, continuous generations, random mating')
# init population
pop = copy.deepcopy(pop0)
for t in range(ngen):
    print('generation ' + str(t) + ' ' + str(len(pop.inds)))
    # step 0: generate marker data for evaluation
    X = gg.do_X(pop.inds, gfeatures, gbase, chips[0], minMaf=0.05)
    # step 1: predicts breeding values
    sel.doEbv(pop, criterion='gblup', X=X, h2=0.3, nh=gfeatures.ploidy)
    # step 2: pedigree with offspring of selected individuals
    newPed = sel.ReturnNewPed(pop, nsel, famsize=noffspring, mating='random',  generation=0)
    # step 3: generates new offspring (this function adds QTN genotypes, true bvs and y)
    pop.addPed(newPed, gfeatures, qtn, gbase)
pop.print()

# PCA plot with individuals colored by generation
X = gg.do_X(pop.inds, gfeatures, gbase, chips[0], minMaf=1e-2)
pca = sel.Pca(X)
pca.fit()
pca.plot(labels=pop.t)

