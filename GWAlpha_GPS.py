#! /usr/bin/env python
import numpy as np
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
nRep=1
h2=.5
nFX=5
nInd=150
nSNP=1000
minFreq=.05
bincuts=(6.4,23.5,76.5,93.6)
lam=40
sig=40
first=True

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
	print("                    -nSim : number of independent simulations to be carried, integer, default is 1")
	print("                    -glm : runs a general linear model between the phenotype and the SNP data, returning a *_glm.txt file with minor allele frequencies, p-values and ranking of associations")
	print("                    -fet : runs a Fisher exact test between the two extreme bins of the data, returning a *_fet.txt file with minor allele frequencies, p-values and ranking of test statistics")
	print("                    -cmh : runs a Cochran-Mantel-Haenszel test between the two extreme bins of multiple replicates of the same data, returning a *_cmh.txt file with minor allele frequencies, p-values and ranking of test statistics")
	print("                    -nRep : number of replicates of a simulation involving the same causal SNPs with the same genetic effects used for the CMH test, integer, default is 1")
	sys.exit()

if "-nSim" in sys.argv:
	nSim=int(sys.argv[sys.argv.index("-nSim")+1])

if "-nRep" in sys.argv:
	nRep=int(sys.argv[sys.argv.index("-nRep")+1])

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

if "-cmh" in sys.argv and "-nRep" not in sys.argv:
	sys.exit("Error: if you want to run the CMH test, you must provide a number of replicates using the -nRep flag")

def GenoPhenoSimulator(Masker=True):
	global nSNP
	global nInd
	global minFreq
	global nFX
	global nRep
	global first
	if nRep>1 and not first:
		global GFX
		global AssoSNP
		global freq
	
	#Generating the allele frequency distributon
	(a,b)=(2,2)
	if nRep==1 or first:
		freq=np.random.beta(a,b,nSNP)

	#Generating the genotypes
	X=np.random.binomial(2,freq,(nInd,nSNP))

	#And removing monomorphic SNPs
	SNPsum=X.sum(axis=0)
	if Masker:
		mask=np.all([SNPsum!=0,SNPsum!=2*nInd],axis=0)
		X=X[:,mask]
		nSNP=X.shape[1]
		freq=freq[mask]

	#estimating the allele frequencies in the data
		freq_hat=np.array(SNPsum[mask],dtype="float")/(2*nInd)
		mask=np.ndarray.flatten(np.array(np.all([freq_hat>minFreq,freq_hat<(1-minFreq)],axis=0)).astype("float"))

	#Computing the coverage based on a Poisson distribution with default parameter lam=40
	#Coverage=np.random.poisson(lam,(len(bincuts)+1)*nSNP).reshape(len(bincuts)+1,nSNP)
	#or based on a lognormal distribution
	Coverage=np.random.lognormal(np.log(lam/np.sqrt((sig**2.)/lam**2.+1)),np.log((sig**2.)/lam**2.+1),(len(bincuts)+1)*nSNP).reshape(len(bincuts)+1,nSNP)
	#or assuming that coverage is even in all bins across loci
	#Coverage=np.random.poisson(lam,nSNP)

	#drawing genetic effects and assigning associated SNPs (if multiple replicates of the experiments are computed, GFX and AssoSNP are recycled from the global env)
	if nRep==1 or first:
		#GFX=np.random.random_sample(nFX)
		GFX=np.random.exponential(2,nFX)
		GFX=GFX/sum(GFX)
		AssoSNP=np.array(np.random.choice(range(nSNP),nFX,replace=False,p=mask/np.sum(mask)),dtype="int")
		first=False

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

	return((X, Y, Coverage, GFX, AssoSNP, freq))

for z in range(nSim):
	
	#Generate the data for the simulated dataset
	(X, Y , Coverage , GFX, AssoSNP, freq)= GenoPhenoSimulator()
	cuts=np.insert(np.percentile(Y,bincuts),0,np.min(Y))
	cuts=np.insert(cuts,len(cuts),np.max(Y))
	try:
		SNP_sync=np.concatenate((np.repeat('Chr',nSNP),np.arange(1,nSNP+1),np.repeat('N',nSNP)),axis=1).reshape(3,nSNP)
	except IndexError:
		SNP_sync=np.concatenate((np.repeat('Chr',nSNP),np.arange(1,nSNP+1),np.repeat('N',nSNP)),axis=0).reshape(3,nSNP)

	#Use a standard LS fit to the data
	if "-glm" in sys.argv:
		from scipy import stats
		beta_SNP=np.array([])
		p_SNP=np.array([])
		for i in range(nSNP):
			slope, intercept, r_value, p_value, std_err = stats.linregress(X[:,i],Y)
			beta_SNP=np.append(beta_SNP,np.round(slope,16))
			p_SNP=np.append(p_SNP,np.round(p_value,16))
		freq_SNP=.5-np.abs(np.sum(X,axis=0)/(2.*nInd)-.5)
		Rk_SNP=1+np.argsort(p_SNP).argsort()
		SNP_glm=np.vstack((np.repeat('Chr',nSNP).astype("|S20"),np.arange(1,nSNP+1),np.repeat('N',nSNP),freq_SNP,beta_SNP,p_SNP,Rk_SNP))

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

	if "-cmh" in sys.argv:
		from scipy import stats
		def calc_cmh_chisq(counts,cont_correct=.5):
			global nRep
			numer_sum = 0
			denom_sum = 0
	
			for i in range(nRep):
				a = counts[0,:,i]
				b = counts[1,:,i]
				c = counts[2,:,i]
				d = counts[3,:,i]
				n = a + b + c + d
				numer_sum += a - (a+b)*(a+c)/(1.0*n)
				denom_sum += (a+b)*(a+c)*(b+d)*(c+d)/(1.0*n*n*(n-1))
	
			chi_sq = ((abs(numer_sum) - cont_correct)*(abs(numer_sum) - cont_correct))/(1.0*denom_sum)
			return(chi_sq)

		low=cuts[1]
		X_low=X[Y<=low,:]
		freq1_low=X_low.sum(axis=0)/float(2*X_low.shape[0])
		freq2_low=1-freq1_low
		count1_low=freq1_low*Coverage[1,:]
		count2_low=freq2_low*Coverage[1,:]
		high=cuts[len(cuts)-2]
		X_high=X[Y>=high,:]
		freq1_high=X_high.sum(axis=0)/float(2*X_high.shape[0])
		freq2_high=1-freq1_high
		count1_high=freq1_high*Coverage[Coverage.shape[0]-1,:]
		count2_high=freq2_high*Coverage[Coverage.shape[0]-1,:]
		counts=np.vstack((count1_low,count2_low,count1_high,count2_high))

		for r in range(nRep-1):
			if X.shape[1]!=nSNP:
				nSNP=X.shape[1]
			(X, Y, Coverage, GFX, AssoSNP, freq)= GenoPhenoSimulator(Masker=False)
			cuts=np.insert(np.percentile(Y,bincuts),0,np.min(Y))
			cuts=np.insert(cuts,len(cuts),np.max(Y))
			low=cuts[1]
			X_low=X[Y<=low,:]
			freq1_low=X_low.sum(axis=0)/float(2*X_low.shape[0])
			freq2_low=1-freq1_low
			count1_low=freq1_low*Coverage[1,:]
			count2_low=freq2_low*Coverage[1,:]
			high=cuts[len(cuts)-2]
			X_high=X[Y>=high,:]
			freq1_high=X_high.sum(axis=0)/float(2*X_high.shape[0])
			freq2_high=1-freq1_high
			count1_high=freq1_high*Coverage[Coverage.shape[0]-1,:]
			count2_high=freq2_high*Coverage[Coverage.shape[0]-1,:]
			counts_int=np.vstack((count1_low,count2_low,count1_high,count2_high))
			counts=np.dstack((counts,counts_int))

		x2_SNP = calc_cmh_chisq(counts)
		p_SNP = stats.chisqprob(x2_SNP,1)
		freq_SNP=.5-np.abs(freq-.5)
		Rk_SNP=1+np.argsort(p_SNP).argsort()
		SNP_cmh=np.vstack((np.repeat('Chr',nSNP).astype("|S20"),np.arange(1,nSNP+1),np.repeat('N',nSNP),freq_SNP,x2_SNP,p_SNP,Rk_SNP))
	
	
	if "-outputName" in sys.argv:
		FileHandle=str(sys.argv[sys.argv.index("-outputName")+1])+"_"+id_generator()
	else:
		FileHandle=id_generator()
	
	filename=FileHandle +"_simul.sync"
	np.savetxt(filename, np.transpose(SNP_sync), delimiter="\t", fmt="%s")
	
	if '-glm' in sys.argv:
		filename=FileHandle +"_glm.txt"
		np.savetxt(filename, np.transpose(SNP_glm), header="Chr\tPosition\tAllele\tFreq\tbeta\tPvalue\tRank", delimiter="\t", fmt="%s")
	
	if '-fet' in sys.argv:
		filename=FileHandle +"_fet.txt"
		np.savetxt(filename, np.transpose(SNP_fet), header="Chr\tPosition\tAllele\tFreq\tPvalue\tRank", delimiter="\t", fmt="%s")
	
	if '-cmh' in sys.argv:
		filename=FileHandle +"_cmh.txt"
		np.savetxt(filename, np.transpose(SNP_cmh), header="Chr\tPosition\tAllele\tFreq\tChisq\tPvalue\tRank", delimiter="\t", fmt="%s")
	
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


