#! /usr/bin/env python
from numpy import *
seterr(all='ignore')
from scipy import stats
from scipy.optimize import minimize
import sys
import re

Pheno_File=sys.argv[1].split('.')[0]
execfile(re.sub(r'_tmp\d+', '', Pheno_File)+"_pheno.py")

Pheno_Dir=sys.argv[1].rsplit('/',1)[0]
Pheno_File=sys.argv[1].rsplit('/',1)[1]

MAF=.03
noP=False

if "-MAF" in sys.argv:
	MAF=float(sys.argv[sys.argv.index("-MAF")+1])

if "noP" in sys.argv:
	noP=True

if noP:
	def PEN(pA):
		return(1)
else:
	def PEN(pA):
		return(((pA*(1-pA))**.5))

if "ML" in sys.argv:
	method="ML"
elif "LS" in sys.argv:
	method="LS"
else: method=raw_input('please indicate which estimatiton method to use (either type ML for maximum likelihood or LS for least-squared): ML/LS: ')

def QuantiSeq(sig,MIN,MAX,perc,q,freqA):
	BIN=insert(perc,len(perc),1)-insert(perc,0,0)
	pA=sum(freqA*BIN)
	BINA=freqA*BIN/pA
	BINA=BINA/sum(BINA)
	BINB=(repeat(1,len(freqA))-freqA)*BIN/(1-pA)
	BINB=BINB/sum(BINB)
	percA=array(BINA[0])
	percB=array(BINB[0])
	for i in range(1,len(BINA)):
		percA=insert(percA,i,sum(BINA[:i+1]))
		percB=insert(percB,i,sum(BINB[:i+1]))
	
	q_prime=(q-repeat(MIN,len(q)))/(MAX-MIN)
	q_prime=insert(q_prime,0,0)
	percA1=delete(percA,len(percA)-1)
	percB1=delete(percB,len(percB)-1)
	percA0=insert(percA1,0,0)
	percB0=insert(percB1,0,0)
	def fn(x):
		if (method=="LS"):
			return(sum(square(percA-stats.beta.cdf(q_prime,x[0],x[1])))+sum(square(percB-stats.beta.cdf(q_prime,x[2],x[3]))))
		else :
			return(-sum(log(stats.beta.cdf(percA,x[0],x[1])-stats.beta.cdf(percA0,x[0],x[1])))-sum(log(stats.beta.cdf(percB,x[2],x[3])-stats.beta.cdf(percB0,x[2],x[3]))))
	try:
		x0 = array([1,1,1,1])
		sol=minimize(fn,x0, method='nelder-mead', options={'xtol': 1e-8, 'disp': False})['x']
		muA_hat=MIN+(MAX-MIN)*sol[0]/(sol[0]+sol[1])
		muB_hat=MIN+(MAX-MIN)*sol[2]/(sol[2]+sol[3])
		Alfa=PEN(pA)*abs(muB_hat-muA_hat)/(2*sig)
	except ValueError:
		Alfa=0
	return(Alfa)

def colsplit(x):
	return(x.split(":"))

perc_bins=[y-x for y, x in zip(perc+[1],[0]+perc)]

Freq_File=open(sys.argv[1],'rb').readlines()
SNP_call=['A','T','C','G','N','D']

GWAlpha_out=[]

for SNP in Freq_File:
	SNP=SNP.replace('\n','').split('\t')
	SNP_name=[SNP[0],SNP[1]]
	SNP_freq=SNP[3:]
	SNP_freq=map(colsplit,SNP_freq)
	SNP_freq=array(SNP_freq,dtype=float)
	freq=SNP_freq/sum(SNP_freq,axis=1)[:,None]+1e-6
	freq_max=amax(freq,axis=0)
	for i in [0,1,2,3,5]:
		if [freq_max>MAF][0][i]:
			if [freq_max!=max(freq_max)][0][i] or sum([freq_max>MAF][0])>2:
				Allele_freq=round(sum(freq[:,i]*perc_bins),3)
				if Allele_freq<MAF: continue
				GWAlpha=QuantiSeq(sig,MIN,MAX,perc,q,freq[:,i])
				GWAlpha_line=[SNP_name[0],SNP_name[1],SNP_call[i],GWAlpha,Allele_freq]
				GWAlpha_out.append(GWAlpha_line)
			elif sum([freq_max==max(freq_max)][0])>1:
				Allele_freq=round(sum(freq[:,i]*perc_bins),3)
				if Allele_freq<MAF: continue
				GWAlpha=QuantiSeq(sig,MIN,MAX,perc,q,freq[:,i])
				GWAlpha_line=[SNP_name[0],SNP_name[1],SNP_call[i],GWAlpha,Allele_freq]
				GWAlpha_out.append(GWAlpha_line)
				break

header="Chromosome,Position,Mutation,Alpha,MAF"
filename=Pheno_Dir+"/GWAlpha_"+Pheno_File+"_out.csv"
savetxt(filename, array(GWAlpha_out), delimiter=",",header=header,fmt="%s")

