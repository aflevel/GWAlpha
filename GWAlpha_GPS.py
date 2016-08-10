#! /usr/bin/env python
import numpy as np
#from scipy import stats
import sys
import os
import string
import random
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
	return ''.join(random.choice(chars) for _ in range(size))

if os.name=='nt':
	brk="\r\n"
else:
	brk="\n"

nSim=1
h2=.5
nFX=5
nInd=150
nSNP=100
minFreq=.05
bincuts=(6.4,23.5,76.5,93.6)
lam=40

status="y"

if len(sys.argv)==1:
	print(""+brk+"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"+brk+"")
	print(""+brk+"Welcome to the GWAlpha Genotype to Phenotype Simulator!!!"+brk+" You are currently running the simulator with default parameters"+brk+" Please type:"+brk+" python GWAlpha_GPS.py -help"+brk+" to print the options"+brk+" or refer to the README file or the manual for more details"+brk)
	print(""+brk+"vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"+brk+"")
	print(brk)
	status = raw_input("Would you like to proceed with default parameters? y/n:")

if "y" not in status and "Y" not in status:
	sys.exit()

if "-help" in sys.argv or "--help" in sys.argv or "-h" in sys.argv:
	print("usage is:"+brk+" python GWAlpha_GPS.py -option_name <value>..."+brk+" and the options are:")
	print("                    -h2 : value for the heritability of the trait, continuous ranging from 0 to 1, default is 0.5")
	print("                    -nSNP : number of SNP to be simulated, integer, default is 100")
	print("                    -nInd : number of individuals in the population, integer, default is 150")
	print("                    -nFX : number of genetic effects (ie SNP affecting the phenotype so that their cumulated effect sums to 1), integer, default is 5")
	print("                    -lambdaCoverage : lambda parameter of the poisson distribution of the sequencing coverage, continuous positive, default is 40")
	print("                    -minFreq : minimun allele frequency for the SNP, continuous from 0 to 1, default is 0.05")
	print("                    -nBins : number of bins the phenotypes were pooled into, this option assumes the bins were of even sizes, integer, default is 5")
	print("                    -speBins : custom specified cut-offs percentage used to define the bins, this option overwrites the -nBins option,"+brk+"                               set of n-1 continous to represent the n bins, default (in %) is 6.4 23.5 76.5 93.6")
	print("                    -nSim : number of simulation to be carried, integer, default is 1")
	print("                    -glm : runs a general linear model between the phenotype and the SNP data, returning a *_glm.txt file with minor allele frequencies, p-values and ranking of associations")
	print("                    -fet : runs a Fisher exact test between the two extreme bins of the data, returning a *_fet.txt file with minor allele frequencies, p-values and ranking of test statistics")
	sys.exit()

if "-nSim" in sys.argv:
	nSim=int(sys.argv[sys.argv.index("-nSim")+1])

if "-h2" in sys.argv:
	h2=float(sys.argv[sys.argv.index("-h2")+1])

if "-nFX" in sys.argv:
	nFX=int(sys.argv[sys.argv.index("-nFX")+1])

if "-nSNP" in sys.argv:
	nSNP=int(sys.argv[sys.argv.index("-nSNP")+1])

if "-nInd" in sys.argv:
	nInd=int(sys.argv[sys.argv.index("-nInd")+1])

if "-lambdaCoverage" in sys.argv:
	lam=float(sys.argv[sys.argv.index("-lambdaCoverage")+1])

if "-minFreq" in sys.argv:
	minFreq=float(sys.argv[sys.argv.index("-minFreq")+1])

if "-nBins" in sys.argv:
	nBins=int(sys.argv[sys.argv.index("-nBins")+1])
	bincuts=np.linspace(0,100,nBins+1)[1:nBins]
	bincuts=tuple(bincuts.reshape(1, -1)[0])

if "-speBins" in sys.argv:
	speBins=sys.argv[sys.argv.index("-speBins")+1:]
	nextargv=[s for s in speBins if "-" in s]
	if (len(nextargv)==0):
		bincuts=np.array(speBins,dtype="float")
		bincuts=tuple(bincuts.reshape(1, -1)[0])
	else:
		endBins = speBins.index(nextargv[0])
		bincuts=np.array(speBins[:endBins],dtype="float")
		bincuts=tuple(bincuts.reshape(1, -1)[0])

if "-VerboseOFF" not in sys.argv:
	print("Running the GWAlpha Genotype to Phenotype Simulator x" + str(nSim) +" times with parameters:")
	print("h2=" + str(h2))
	print("nSNP=" + str(nSNP))
	print("nInd=" + str(nInd))
	print("nFX=" + str(nFX))
	print("minFreq=" + str(minFreq))
	print("lambdaCov=" + str(lam))
	print("BinCuts=" + str(bincuts).replace(""+brk+"",""))

if nFX>nSNP:
	sys.exit("Error: you are trying to model more genetic effects than SNP.")

if len(bincuts)+2>nInd:
	sys.exit("Error: you are trying to model more phenotype bins than individuals.")

if len(bincuts)==0:
	print("Warning: GWAlpha will not run with less than two bins")

for z in range(nSim):

	#Generating the allele frequency distributon
	(a,b)=(2,2)
	freq=np.random.beta(a,b,nSNP)

	#Generating the genotypes
	X=np.random.binomial(2,freq,(nInd,nSNP))

	#And removing monomorphic SNPs
	SNPsum=X.sum(axis=0)
	mask=np.all([SNPsum!=0,SNPsum!=2*nInd],axis=0)
	X=X[:,mask]
	nSNP=X.shape[1]

	#estimating the allele frequencies in the data
	freq_hat=np.array(SNPsum[mask],dtype="float")/(2*nInd)
	mask=np.ndarray.flatten(np.array(np.all([freq_hat>minFreq,freq_hat<(1-minFreq)],axis=0)).astype("float"))

	#Computing the coverage based on a Poisson distribution with default parameter lam=40
	Coverage=np.random.poisson(lam,(len(bincuts)+1)*nSNP).reshape(len(bincuts)+1,nSNP)
	#or assuming that coverage is even in all bins across loci
	#Coverage=np.random.poisson(lam,nSNP)

	#drawing genetic effects and assigning associated SNPs
	GFX=np.random.random_sample(nFX)
	GFX=GFX/sum(GFX)
	AssoSNP=np.array(np.random.choice(range(nSNP),nFX,replace=False,p=mask/np.sum(mask)),dtype="int")

	#partionning the variance taking into account the linkage among associated SNPs
	AssoX=X[:,AssoSNP]
	Rho=np.corrcoef(AssoX,rowvar=0)
	XtX=GFX.reshape(1,nFX)*GFX.reshape(nFX,1)
	Vg=np.sum(Rho*XtX)/4
	Ve=Vg*(1/h2-1)

	#Generating the phenotypes based on the variance components
	G=AssoX*GFX
	Xb=G.sum(1)
	e=np.random.normal(0,Ve**(0.5),nInd)
	Y=Xb+e

	#Generate the data for the simulated dataset
	cuts=np.insert(np.percentile(Y,bincuts),0,np.min(Y))
	cuts=np.insert(cuts,len(cuts),np.max(Y))
	SNP_sync=np.concatenate((np.repeat('Chr',nSNP),np.arange(1,nSNP+1),np.repeat('N',nSNP)),axis=1).reshape(3,nSNP)

	#Use a standard LS fit to the data
	if "-glm" in sys.argv:
		from scipy import stats
		p_SNP=np.array([])
		for i in range(nSNP):
			slope, intercept, r_value, p_value, std_err = stats.linregress(X[:,i],Y)
			p_SNP=np.append(p_SNP,np.round(p_value,16))
		freq_SNP=.5-np.abs(np.sum(X,axis=0)/(2.*nInd)-.5)
		Rk_SNP=1+np.argsort(p_SNP).argsort()
		SNP_glm=np.vstack((np.repeat('Chr',nSNP).astype("|S20"),np.arange(1,nSNP+1),np.repeat('N',nSNP),freq_SNP,p_SNP,Rk_SNP))

	#Use a Fisher exact test to the data
	if "-fet" in sys.argv:
		from scipy import stats
		low=cuts[1]
		X_low=X[Y<=low,:]
		freq_low=X_low.sum(axis=0)/float(2*X_low.shape[0])
		high=cuts[len(cuts)-2]
		X_high=X[Y>=high,:]
		freq_high=X_high.sum(axis=0)/float(2*X_high.shape[0])
		def FET(i):
			table=[[int(freq_low[i]*Coverage[0,i]),int((1-freq_low[i])*Coverage[0,i])],[int(freq_high[i]*Coverage[0,i]),int((1-freq_high[i])*Coverage[0,i])]]
			return(stats.fisher_exact(table,alternative='two-sided')[1])
		p_SNP=np.array([])
		for i in range(nSNP):
			p_value= FET(i)
			p_SNP=np.append(p_SNP,np.round(p_value,16))
		freq_SNP=.5-np.abs(np.sum(X,axis=0)/(2.*nInd)-.5)
		Rk_SNP=1+np.argsort(p_SNP).argsort()
		SNP_fet=np.vstack((np.repeat('Chr',nSNP).astype("|S20"),np.arange(1,nSNP+1),np.repeat('N',nSNP),freq_SNP,p_SNP,Rk_SNP))

	for k in range(1,len(cuts)):
		bottom=cuts[k-1]
		top=cuts[k]
		X_bin=X[np.all([Y>=bottom,Y<=top],axis=0),:]
		freq_bin=X_bin.sum(axis=0)/float(2*X_bin.shape[0])
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

	FileHandle=id_generator()
	filename=FileHandle +"_simul.sync"
	np.savetxt(filename, np.transpose(SNP_sync), delimiter="\t", fmt="%s")
	
	if '-glm' in sys.argv:
		filename=FileHandle +"_glm.txt"
		np.savetxt(filename, np.transpose(SNP_glm), header="Chr\tPosition\tAllele\tFreq\tPvalue\tRank", delimiter="\t", fmt="%s")
	
	if '-fet' in sys.argv:
		filename=FileHandle +"_fet.txt"
		np.savetxt(filename, np.transpose(SNP_fet), header="Chr\tPosition\tAllele\tFreq\tPvalue\tRank", delimiter="\t", fmt="%s")
	
	PERC=str(np.array(bincuts)/100)
	CUTS=str(cuts[1:(len(cuts)-1)])
	ASSOSNP=str(AssoSNP)
	FXSNP=str(GFX)
	while "  " in PERC:
		PERC=PERC.replace("  "," ")
	PERC=PERC.replace(" ",",").replace("[,","[").replace(",]","]")
	while "  " in CUTS:
		CUTS=CUTS.replace("  "," ")
	CUTS=CUTS.replace(" ",",").replace("[,","[").replace(",]","]")
	while "  " in ASSOSNP:
		ASSOSNP=ASSOSNP.replace("  "," ")
	ASSOSNP=ASSOSNP.replace(" ",",").replace("[,","[").replace(",]","]")
	while "  " in FXSNP:
		FXSNP=FXSNP.replace("  "," ")
	FXSNP=FXSNP.replace(" ",",").replace("[,","[").replace(",]","]")

	PHENO = open(FileHandle+'_simul_pheno.py', 'w')
	PHENO.write('#Name of the phenotype'+brk+'Pheno_name=\"'+FileHandle+'\"'+brk)
	PHENO.write('#Standard deviation'+brk+'sig='+str(np.std(Y))+brk)
	PHENO.write('#Minimum value for the phenotype'+brk+'MIN='+str(np.min(Y))+brk)
	PHENO.write('#Maximum value of the phenotype'+brk+'MAX='+str(np.max(Y))+brk)
	PHENO.write('#Position of the percentile cut-off positions'+brk+'perc=' + PERC + brk)
	PHENO.write('#Phenotype values corresponding to each percentile cut-off positions' + brk+ 'q=' + CUTS + brk)
	PHENO.write('#indices of the Associated SNPs, starting from 0 in the sync file'+ brk +'AssoSNP=' + ASSOSNP + brk)
	PHENO.write('#Effects of the Associated SNPs, summing to 1'+brk+'FXSNP='+ FXSNP +brk)
	PHENO.write('#Number of individuals simulated'+brk+'nInd='+ str(nInd) +brk)
	PHENO.close()


