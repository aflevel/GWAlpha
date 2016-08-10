#! /usr/bin/env python
import numpy as np
import sys
import os

np.set_printoptions(linewidth=500)

if os.name=='nt':
	brk="\r\n"
else:
	brk="\n"

bincuts=(6.4,23.5,76.5,93.6)
lam=40

if "-help" in sys.argv or "--help" in sys.argv or "-h" in sys.argv:
	print("usage is:"+brk+" python GWAlpha_RDC.py [INPUT_PHENOTYPE] [INPUT_GENOTYPE] -option_name <value>..."+brk+" [INPUT_PHENOTYPE]: a *.csv file containing the trait values"+brk+" [INPUT_GENOTYPE]: a *.csv file containing the genotype values"+brk+" and the options are:")
	print("                    -nBins : number of bins the phenotypes were pooled into, this option assumes the bins were of even sizes, integer, default is 5")
	print("                    -speBins : custom specified cut-offs percentage used to define the bins, this option overwrites the -nBins option,"+brk+"                               set of n-1 continous to represent the n bins, default (in %) is 6.4 23.5 76.5 93.6")
	sys.exit()

if "-nBins" in sys.argv:
	nBins=int(sys.argv[sys.argv.index("-nBins")+1])
	bincuts=np.linspace(0,100,nBins+1)[1:nBins]
	bincuts=tuple(bincuts.reshape(1, -1)[0])

#Missing Data Threshold
Miss_Tol=.3

#get the SNP data
SNP_handle=open(sys.argv[2],'r')
Lines_geno=np.array(SNP_handle.readline().replace('\n','').split(',')[2:])
X = SNP_handle.readlines()
X = np.array([x.replace('\n','').split(',') for x in X])
SNP_names=X[:,:2]
X=X[:,2:]

#get the phenotype data
Pheno_handle=open(sys.argv[1],'r')
Y = Pheno_handle.readlines()
Y = np.array([x.replace('\n','').split(',') for x in Y])
Lines_pheno=Y[1:,0]
Y=Y[1:,1].astype("float")

#remove lines with missing phenotype
Y=Y[~np.isnan(Y)]
Lines_pheno=Lines_pheno[~np.isnan(Y)]
#sort the lines in Phenotype
Pheno_sort=np.argsort(Lines_pheno)
Y=Y[Pheno_sort]
#remove the lines with no genotypes
Lines_in_pheno=np.in1d(Lines_pheno,Lines_geno)
Y=Y[Lines_in_pheno]
Lines_pheno=Lines_pheno[Lines_in_pheno]

#remove lines with missing phenotype and too many missing genotypes
Lines_in_geno=np.in1d(Lines_geno,Lines_pheno)
X=X[:,Lines_in_geno]
miss_mask=np.array(2*np.sum(X=='nan',axis=1),dtype='float')/X.shape[1]
X=X[miss_mask<Miss_Tol]
SNP_names=SNP_names[miss_mask<Miss_Tol,:]

#sort the lines in SNP
SNP_sort=np.argsort(Lines_geno[Lines_in_geno])
X=np.transpose(X[:,SNP_sort])
X=X.astype('float')

#And removing monomorphic SNPs
SNPsum=np.nansum(X,axis=0)
nInd=np.sum(~np.isnan(X),axis=0)
mask_fix=np.all([SNPsum!=0,SNPsum!=2*nInd],axis=0)
X=X[:,mask_fix]
SNP_names=SNP_names[mask_fix,:]
nSNP=X.shape[1]

#Computing the coverage based on a Poisson distribution with default parameter lam=40
#Coverage=np.random.poisson(lam,(len(bincuts)+1)*nSNP).reshape(len(bincuts)+1,nSNP)
#OR even coverage of 200X
Coverage=np.repeat(200,(len(bincuts)+1)*nSNP).reshape(len(bincuts)+1,nSNP)

#Generate the data for the simulated dataset
cuts=np.insert(np.percentile(Y,bincuts),0,np.min(Y))
cuts=np.insert(cuts,len(cuts),np.max(Y))

SNP_sync=np.concatenate((SNP_names[:,0],SNP_names[:,1],np.repeat('N',nSNP)),axis=1)
SNP_sync=SNP_sync.reshape(3,nSNP)

for k in range(1,len(cuts)):
	bottom=cuts[k-1]
	top=cuts[k]
	X_bin=X[np.all([Y>=bottom,Y<=top],axis=0),:]
	freq_bin=np.nansum(X_bin,axis=0)/np.array(2*np.sum(~np.isnan(X_bin),axis=0),dtype='float')
	sync_bin=np.array(Coverage[k-1,:]*freq_bin,dtype='int')
	alt_bin=np.array(Coverage[k-1,:]*(1-freq_bin),dtype='int')
	SNP_sync=np.concatenate((SNP_sync,
			                 np.core.defchararray.add(
			                 np.core.defchararray.add(
			                 np.core.defchararray.add(np.array(sync_bin,dtype='str'),
			                                          np.repeat(':',nSNP)),
			                                          np.array(alt_bin,dtype='str')),
			                                          np.repeat(':0:0:0:0',nSNP)).reshape(1,nSNP)
			               ),axis=0)

SNP_sync=SNP_sync[:,1:]

filename=sys.argv[1].replace(".csv","").replace("_","") +"w"+str(len(bincuts)+1)+"bins.sync"
np.savetxt(filename, np.transpose(SNP_sync), delimiter="\t", fmt="%s")

PERC=str(np.array(bincuts)/100)
CUTS=str(cuts[1:(len(cuts)-1)])
while "  " in PERC:
	PERC=PERC.replace("  "," ")
PERC=PERC.replace(" ",",").replace("[,","[").replace(",]","]")
while "  " in CUTS:
	CUTS=CUTS.replace("  "," ")
CUTS=CUTS.replace(" ",",").replace("[,","[").replace(",]","]")

PHENO = open(sys.argv[1].replace(".csv","").replace("_","") +"w"+str(len(bincuts)+1)+'bins_pheno.py', 'w')
PHENO.write('#Name of the phenotype'+brk+'Pheno_name=\"' + sys.argv[1] +'\"'+brk)
PHENO.write('#Standard deviation'+brk+'sig='+str(np.std(Y))+brk)
PHENO.write('#Minimum value for the phenotype'+brk+'MIN='+str(np.min(Y))+brk)
PHENO.write('#Maximum value of the phenotype'+brk+'MAX='+str(np.max(Y))+brk)
PHENO.write('#Position of the percentile cut-off positions'+brk+'perc=' + PERC + brk)
PHENO.write('#Phenotype values corresponding to each percentile cut-off positions' + brk+ 'q=' + CUTS + brk)
PHENO.close()


