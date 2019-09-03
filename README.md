## [SeqBreed: A python tool to evaluate genomic prediction in complex scenarios](https://www.biorxiv.org/content/10.1101/748624v1)
#### [Miguel Perez-Enciso](miguel.perez@uab.es)
With contributions from [Laura M Zingaretti](m.lau.zingaretti@gmail.com) and [Lino C. Ramirez-Ayala](linocesar.ramirez@cragenomica.es)

---

SeqBreed.py is a python3 software to simulate populations under (genomic) selection. 
It inherits most of funcionality from SBVB ([Perez-Enciso et al., 2017](http://www.genetics.org/content/205/2/939.long), 
https://github.com/miguelperezenciso/sbvb1) and pSBVB ([Zingaretti et al. 2018](http://www.g3journal.org/content/9/2/327.long),
 https://github.com/lauzingaretti/pSBVB) softwares but code has been completely rewritten. It is now much more user friendly and adds numerous new features such as easy selection implementation (BLUP, GBLUP, SSTEP), plots (PCA, GWAS)... 

Its main target is to simulate genomic selection experiments but can be used as well to study the performance of GWAS or,
in general, study the dynamics of complex traits under numerous selection criteria: drift, mass selection, 
BLUP, GBLUP and single step GBLUP are currently implemented. 

It can accommodate any number of complex phenotypes controlled by an arbitrary number of loci 
(epistasis is not currently implemented though). Autopolyploid genomes can be simulated, as well as sex and mitochondrial chromosomes.

**NOTE:** SeqBreed is designed for use in short term breeding experiments only, as no new mutations are generated. We recommend to use real sequence data or high density SNP data as startup to realistically mimic variability and disequilibrium.  

### Requirements
See [requirements.txt](https://github.com/miguelperezenciso/SeqBreed/blob/master/requirements.txt).

### Last changes
* Dihaploid offspring can be generated now (https://github.com/miguelperezenciso/SeqBreed/blob/master/README.md#5-breeding-population)

### Installation
Clone the repository and run:

`sudo python3 -m pip install SeqBreed-XXX.whl`

Alternatively, you can install locally without sudo permits as:

`pip3 install --user SeqBreed-XXX.whl`

To access the module, include the following in the program: 

    from SeqBreed import genome as gg
    from SeqBreed.selection import selection as sel
    
Source code is provided in src, just in case. Files must be in working directory. To access them:

    import genome as gg
    import selection as sel

Program has been tested in mac and linux only, although it should run as any regular python script in windows.

As in any python project, it is recommended to use a separate environment to avoid conflicts between package versions. You can do that with **conda** as follows:

```
conda create -n seqbreed_environment
conda activate seqbreed_environment
   # do all instalation of required packages, eg,
   python3 -m pip install --user SeqBreed-XXX.whl
   # run python
   ....
# to finish
conda deactivate
```
### Quick startup
The basic phylosophy of **SeqBreed** is to have a file with SNP data from the founder population (in vcf or plink - like format), specify causal SNPs (QTNs) and their effects for every phenotype (either in a file or can be simulated by the program) and, optionally, a pedigree that is used for gene-dropping. In addition to sequence data, the user can specify a subset of SNPs (a chip) that can be used to implement genomic selection, do a PCA or a GWAS. Next, new individuals can be manually added to the extant population or a selection scheme can be automatically implemented. At any stage, data can be inspected, exported or plotted via a Principal Component Analysis (PCA). 

The minimum input for **SeqBreed** is the file with SNP genotypes and positions from the founder population. Addiitonal files mentioned add extra control but can be simulated (eg, QTN positions and effects or recombination rates). 

**WARNING:** All QTNs and all chip SNPs must be in the vcf/plink file, they are removed otherwise. 


### Examples
The best option to use the program is to follow the examples.
* [`main.py`](https://github.com/miguelperezenciso/SeqBreed/blob/master/main.py) and [`SeqBreed_tutorial.ipynb`](https://github.com/miguelperezenciso/SeqBreed/blob/master/SeqBreed_tutorial.ipynb) illustrates main functionalities.
* [`POTATO`](https://github.com/miguelperezenciso/SeqBreed/tree/master/POTATO) folder contains data from tetraploid potato and illustrates how to generate an F2, how to do a GWAS, 
PCA-corrected GWAS, how to simulate additional phenotypes...
* [`DGRP`](https://github.com/miguelperezenciso/SeqBreed/tree/master/DGRP) folder contains data from the Drosphila Genome Reference Panel (http://dgrp2.gnets.ncsu.edu/)
project and illustrates genomic selection, how to save and reuse big files with ```pickle```, etc. 


### Main classes
**SeqBreed** allows storing and accessing genomic and population information easily. 

* ```Population``` class allows accessing most items of interest. It has a series of methods allowing generating new individuals, including the different selection steps. 
* ```GFounder``` class allows storing genotypes from base population individuals, ie, those in the vcf / gen file.
* ```Genome``` contains all relevant genome features such as number of chromosomes, recombination rates, sequence snp positions...
* ```Chip``` contains a list of SNP positions.

Besides, the ```Selection``` module contains methods to perform GWAS, Genomic Prediction, etc.

### Input
**SeqBreed** minimally only requires a genotype file in vcf or 'gen' formats, either compressed or uncompressed.
For vcf format specifications, check https://en.wikipedia.org/wiki/Variant_Call_Format. The 'gen' format is 
simply 

    chr base_pair indiv1_allele1  indiv1_allele2  indiv2_allele1 ...

one row for each SNP, where alleles are coded as ```0, 1``` for diploids and ```0...(ploidy-1)``` for polyploids. 
For example, the following gen file contains information for 3 markers and 2 diploid individuals:

    chr1   1000 1 0 1 1
    chrX   100034 1 1 0 0
    chrX   200000 1 1 0 1

This format can be generated easily from plink .gen format file.

**NOTE:** **SeqBreed** recognizes vcf format if file name ends with 'vcf' or 'vcf.gz', everything else is
treated as 'gen' format. **SeqBreed** automatically recognizes whether the files are gz compressed.

**NOTE:** **SeqBreed** automatically recognizes ploidy from vcf files. Ploidy must be specified for gen files (see below).

**IMPORTANT: No missing genotypes are allowed and SeqBreed assumes genotypes to be phased.** This is simply because **SeqBreed** must know which are the parental haplotypes. Since this is unlikely, phasing and imputation with some standard algorithm will do. Errors in imputation and phasing should have virtually no effect on the conclusions to be drawn from **SeqBreed** simulations.

If no QTN file is provided, a predetermined number of randomly chosen QTNs can be generated by **SeqBreed**,
where additive qtn effects are sampled from a gamma. For more control, a qtn file is needed.  

Optional Files (see more details below and examples in the repository):
 
- **QTN file** specifying the genetic architecture: any number of loci and any number of traits are allowed. 
The format of this file is, per row, QTN position (chromosome and base pair) followed by additive and 
dominant effect for each trait. If a QTN does not affect the given trait, 0's must be employed.
- **File with (sex) specific recombination maps**: Sex and mitochondrial chromosomes can be specified. 
Auto polyploid genomes can also be specified. 
- **A starting pedigree**: If not provided, a new one is automatically generated where base individuals are unrelated. 

### Usage
A typical SeqBreed run consists of at least the following steps:

1-	Upload founder sequence genotypes and a `GFounder` object is created. A file with all SNP positions in sequence is generated.

2-	Initialize `Genome` class. Optionally, sex or mitochondrial chromosomes are specified as well as local recombination maps.

3-	Genetic architectures for every trait are specified via a `QTNs` object. Environmental variances are also inferred.

4-	A `Population` object is generated, optionally via gene-dropping along a predetermined pedigree.

Once `Population` is initialized, SeqBreed allows a number of operations to be performed, such as implementing selection, plotting or exporting phenotype and molecular data.

#### 1. Base population
`GFounder` class constructor simply needs vcf file, name of output file with snp positions and minimimum frequency (MAF) for a snp to be considered.

```gbase = gg.GFounder(vcfFile, snpFile, minMaf=0)```

It automatically determines no. of base individuals, ploidy (in vcf files), genotypes and allele frequencies.
For gen files, say 'test.gen', ploidy (say octoploid) must be specified:

```gbase = GFounder(vcfFile='test.gen', snpFile, minMaf=0, ploidy=8)```

In this case, 8 alleles by individual and SNP must be present.

**NOTE:** If unspecified, default ploidy is two. 

#### 2. Genome features
```Genome``` class allows specifying genome features

```gfeatures = gg.Genome(snpFile, mapFile=None, XChr=None, YChr=None, MTChr=None, ploidy=gbase.ploidy)```

internally, it holds info on snp positions for each chromosome, among other variables. Ploidy is automatically
inferred in previous step (```gbase = GFounder(.)```). The mapfile has the following format:

```chr base_pair rec_rate_males rec_rate_females```

where recombination rate is in cM / Mb, and refers to the genome region up to ```base_pair```. Rec rate is set to default
(1 cM/Mb) in missing regions from the map file. this value can be modified in the code (check ```Chromosome``` class 
constructor cM2Mb = 1).

#### 3. Specifying genetic architecture
QTN class allows specifying causative SNPs and their effects on an unlimited number of phenotypes. ```h2``` is an
array with dimension the number of traits that specify heritability for each trait.

```qtn = QTN(h2=[h2, ...], genome=gfeatures, qtnFile=None, se=None, nqtn=0, name=None)```

QTNS can be specified in three ways:

**OPTION 1:** A random number of ```nqtn``` additive qtls genome wide distributed, additive effects are sampled 
from a gamma, valid only for one phenotype.

```qtn = gg.QTNs(h2=[0.7], genome=gfeatures, nqtn=10)```

**OPTION 2:** qtnFile contains only QTN positions. add effects are sampled form a gamma, only one phenotype

```qtn = gg.QTNs(h2=[0.5], genome=gfeatures, qtnFile=qtnfile1)```

qtnFile format is:

```chr base_pair```

One line for each QTN

**OPTION 3:** qtn positions and effects are read from qtnFile, any number of phenotypes and additive /
dominance action can be specified.

```qtn = gg.QTNs(h2=[0.9, 0.3], genome=gfeatures, qtnFile=qtnfile2)```

qtnFile format is:

```chr base_pair add_eff_trait1 dom_eff_trait1 add_eff_trait2 dom_eff_trait2 ... ```

One line for each QTN.

#### 4. Genotyping Chips
The ```Chip``` class allows efficient manipulation and storing of a genotyping array. Chips can be defined either by
reading SNP positions from a file (```chipFile```) or can be generated from a random or uniform set of SNPs. In the 
latter case, the number of SNPs in the chip must be specified and a minimum MAF can also be set.

    chip = gg.Chip(gfeatures, chipFile=None, nsnp=0, unif=False, minMaf=0, name='chip_') 
    
where:

- ```gfeatures``` is a Genome object
- ```gbase``` is a GFounders object
- ```chipFile```(str) is the name of pos file containing chip snp posiitons
- ```nsnp```(int) no. of snps to be generated randomly or uniformly
- ```unif```(bool) if True, snps chosen are uniformly distributed; randomly otherwise
- ```minMaf```(float) minimum MAF
- ```name```(str) chip name

If chipFile is specified, it overrides the generation option. `Chip` objects can be generated from a file containing SNP positions (`chipfile`):

`chip = gg.Chip(genome=gfeatures, chipFile=chipfile) `

or randomly generated. The following generates a chip with N equally spaced SNPs with MAF>maf

`chip = gg.Chip(gfeatures, gbase, nsnp=N, unif=True, minMaf=maf)`

The following generates a chip with N randomly chosen SNPs without any restriction on MAF:

`chip = gg.Chip(gfeatures, gbase, nsnp=N, unif=False)`


#### 5. Breeding population
Most of information needed is in ```Population``` class, which is basically a collection of individuals plus
some little extra information.

```pop = Population(genome=gfeatures, pedFile=None, t=[], label=None, qtns=qtn, gfounders=gbase)```

pedFile contains pedigree information, ie,

```id id_father id_mother```

where ids **MUST** be integer consecutive numbers from 1 onwards; id_father and id_mother **MUST** be 0 for
all base individuals, ie, the number of individuals in the vcf file. If pedFile is not specified, 
a dummy pedigree for all base individuals is generated (```id 0 0``` for id=1:nbase).

Individuals can be accessed by 

```Pop.inds[i]```

where ```Pop.inds[i]``` contains information of individual ```id=i+1``` (python array indices start at 0). Phenotypes can
be accessed by 

```Pop.inds[i].y```

Pop and ind info can be printed with:

```
# prints generic info
pop.print()
# prints detailed indiv info
pop.inds[i].print(gfeatures)
# boxplots average phenotype values per generation for phenotype itrait (= 0:ntrait)
pop.plot(trait=itrait) 
```

It is usually difficult to find real sequence data to obtain a reasonably sized founder (base) population. 
An interesting feature of **SeqBreed** is the possibility of generating ‘dummy’ founder individuals 
by randomly combining recombinant founder haplotypes. The following function adds a randomly generated individual:

```pop.addRandomInd(gfeatures, gbase.nbase, k=5, mode='pedigree', qtns=qtn, gfounders=gbase)```

where ```mode``` can be 'pedigree' or 'random', and ```k``` specifies the number of recombination generations. If 
```mode``` is 'pedigree', a random pedigree consisting of 2^k founder individuals and k generations is simulated, and
genedropping is performed along this pedigree. The resulting individual is added to the ```pop``` object.
If ```mode``` is 'random', a recombinant chromosome with x ~ Poisson(0.5 k L), L being genetic lenth, recombinant breaks
is simulated, and each non-recombinant stretch is assigned a random founder haplotype.

**Adding extra individuals** manually is also possible. The following adds a male offspring (sex=0) from 24th and first individual:

    parents = [pop.inds[23], pop.inds[0] ]
    pop.addInd(parents, genome=gfeatures, gfounders=gbase, qtns=qtn, sex=0)

The following generates a **dihaploid** offspring of last individual in population:

    parents = [pop.inds[-1], pop.inds[-1]]
    pop.addInd(parents, dihap=True, genome=gfeatures, gfounders=gbase, qtns=qtn, sex=0)

The following adds a set of new individuals structured as a pedigree:

    pop.addPed(ped, genome=gfeatures, qtns=qtn, gfounders=gbase)

where `ped` is an n x 3 or n x 4 array containing `id, father, mother, sex` per row. The last column (`sex`) is optional. In this case, `sex` is randomly sampled.

#### 6. Retrieving genotype data 
**SeqBreed** internally keeps individual genomes only as a collection of recombination blocks and founder origins of
those blocks ([Perez-Enciso et al. 2000](https://gsejournal.biomedcentral.com/articles/10.1186/1297-9686-32-5-467)), therefore enormoulsy saving cpu time and memory in all operations on Individual objects. 
Actual genotypes must be retrieved though for GBLUP, GWAS or PCA analyses. This is achieved with ```do_X``` function

```X = do_X(inds, gfeatures, gbase, chip, minMaf=1e-6)```

where:

- ```X``` is an integer p x n np.array, p = n. markers, n = n inds.
- ```inds``` is an array of Individual objects (say ```pop.inds```)
- ```gfeatures``` is a ```Genome``` object
- ```gbase``` is a ```GFounders``` object
- ```chip``` is a ```Chip``` object
- ```minMaf``` is the minimum allele frequency for a SNP to be considered [1e-6]

**WARNING!!!!:** ```do_X``` has been parallelized with ```cython.parallel``` module but still can be expensive and memory
 demanding. Beware with large SNP datasets.

#### 7. Implementing selection
Selection proceeds in three steps (see code above):
* First, estimated breeding values (EBVs) are obtained for the current
population. If molecular data are needed (GBLUP, ssGBLUP), an X matrix containing genotypes must be 
generated with ```do_X``` function. 
* Second, the best individuals are mated and offspring are generated. 
Computationally, this is done by extending the current pedigree. 
* Third, new genomes and phenotypes are generated for the new offspring (but not their EBVs). 

Estimated Breeding Values are obtained with function

```sel.doEbv(pop, criterion='random', h2=None, X=None, mkrIds=None, yIds=[], nh=2, itrait=0)```

where:

- ```criterion``` (str) is one of: 'random', 'phenotype', 'blup', 'gblup' or 'sstep' (single step).
- ```h2``` is heritability used for genetic evaluation (needed for blup, gblup and sstep).
- ```X``` is an array with genotypes (needed for gblup and sstep), and can be obtained with ```do_X``` function
- ```mkrIds``` (integer numpy array) indices of genotyped individuals (indexed starting with 0) [required for sstep]
- ```yIds``` (numpy int array): integer array specifying individuals with data (indexed starting with 0) [all]
- ```trait```(int): trait index for which evaluation is performed [0]

**In sstep, marker files should contain only information for genotyped individuals and in the same order.**

Function ```doEBV``` assigns EBVs to ```pop.inds[:].ebv```. Users can replace the function by any of their
choice such that selection is based on the custom defined criterion. For instance, if vector ```ebvs``` contains the custom values,
then

    nind = len(pop.inds)
    for i in range(nind): pop.inds[i].ebv = ebvs[i]

assigns new EBVs to ```pop``` object. The next step is to generate offspring from selected parents. 
This is performed with function 

```
newPed = sel.ReturnNewPed(pop, nsel, famsize, mating='random', generation=0)
```

where:

- ```nsel``` (int numpy two element array) contains the number of males and females to be selected as parents.
- ```famsize```(int): number of offspring per female
- ```mating```(str): 'random' (=='r') or 'assortative' (=='a'). 'random' means males and females are mated at random; 'assortative' means that best males are mated to best females.
- ```generation```(int): only individuals from generation onwards are considered as potential parents [all = 0] ,
this can be used to specify discrete or continuous generations.

`ReturnNewPed` method returns offspring of selected parents (```newPed``` is an array with this pedigree), and `pop` object is extended with ```Population``` method `addPed` as follows:

    pop.addPed(newPed, gfeatures, qtn, gbase, t=None)
    
where:

- ```gfeatures``` is a ```Genome``` object
- ```qtn``` is a ```QTNs``` object
- ```gbase``` is a ```GFounders``` object
- ```t``` is an integer with generation, and stored in vector ```pop.t``` , by default, t is increased by one and all
new individuals are assigned ```t = max(current pop.t) + 1```.

This method actually generates the genomes and phenotypes of new offspring and adds it to `pop` object. Generation counter `t` is automatically increased by 1.

See accompanying script ```main.py``` for examples. 
  
#### 8. PCA and GWAS
Both PCA and GWAS plots can be obtained, GWAS simply performs a regression SNP by SNP. In either case, you must first 
generate the genotypic data X (This has been parallelized but can be expensive though).

```
X = gg.do_X(pop.inds, gfeatures, gbase, chip=chip1)
```

```
# Returns 2D PCA
pca = sel.Pca(X)
pca.fit()
pca.plot(labels=labels) # labels is an integer class vector of dimension # of inds to plot in colors
```

```
itrait = 0    # first trait is analyzed
gwas = sel.Gwas(X, chip1)  # WARNING: chip1 must match that used to get X
gwas.fit(inds=pop.inds, trait=itrait)

# alternatively, you can pass the phenotypes directly
y = np.array(list(pop.inds[i].y[itrait] for i in range(pop.n)))
gwas.fit(y=y, trait=itrait)

gwas.plot(fdr=True)        # plots FDR
gwas.plot()                # plots pvalue
gwas.print(gfeatures)      # prints gwas results
```
By default, `gwas` uses raw phenotypes. Example in [potato.py](https://github.com/miguelperezenciso/SeqBreed/blob/master/POTATO/potato.py) illustrates how to correct for PC components:

```
pca = sel.Pca(X)
pca.fit()
# get phenotype vector for given trait
itrait=0
y = np.array(list(ind.y[itrait] for ind in pop.inds))
# correct phenotypes using first two PCs
y = y - pca.p[:,0] - pca.p[:,1] 
gwas.fit(y=y, trait=itrait)
gwas.plot()
```

### Citation
Please cite this if you use or reuse the code:

M. Perez-Enciso, L.C. Ramirez-Ayala,  L.M. Zingaretti. 2019. SeqBreed: a python tool to evaluate genomic prediction 
in complex scenarios. [Submitted](https://www.biorxiv.org/content/10.1101/748624v1). 

### How to contribute
Please send comments, suggestions or report bugs to miguel.perez@uab.es. From 2020 on, I will be working at
INIA in Madrid (Spain), check at [INIA](www.inia.es) website or in the internet for my new email.

#### Disclaimer
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

