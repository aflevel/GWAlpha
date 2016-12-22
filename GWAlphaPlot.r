#!/usr/bin/Rscript

args=commandArgs(TRUE)
options(warn=-1)

RAW=T

if ("pval" %in% args) {
	RAW=F
	if ("qqplot" %in% args) QQ=T
}

LL=function(x,y) sum(dnorm(y,x[1],x[2],log=T))

Pheno.name=unlist(strsplit(args[1],"/"))
Pheno.name=Pheno.name[length(Pheno.name)]
Pheno.name=gsub("GWAlpha_","",Pheno.name)
Pheno.name=gsub("_out.csv","",Pheno.name)

GWAlpha=read.csv(args[1],header=T)
GWAlpha=GWAlpha[!is.na(as.numeric(as.character(GWAlpha[,ncol(GWAlpha)]))),]

if ("2L" %in% GWAlpha[,1]) {
CHR.labs=c("5"="4","1"="2L","2"="2R","3"="3L","4"="3R","6"="X")
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
} else {
CHR.labs=unique(as.character(GWAlpha[,1]))
names(CHR.labs)=unique(as.character(GWAlpha[,1]))
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
}

OUT=GWAlpha[GWAlpha[,1] %in% CHR.labs,]
CHR=OUT[,1]

for (i in 1:length(CHR.labs)) CHR=gsub(CHR.labs[i],names(CHR.labs)[i],CHR)
CHR.labs=CHR.labs[names(CHR.labs) %in% unique(CHR)]
CHR.labs=CHR.labs[order(names(CHR.labs))]

Position=as.numeric(as.character(OUT[,2]))
CHR.start=vector()
CHR.end=vector()
for (i in unique(CHR)) {
	ID=which(CHR==i)
	CHR.start=c(CHR.start,ID[1])
	CHR.end=c(CHR.end,ID[length(ID)])
}
CHR.summary=cbind(Position[CHR.start],Position[CHR.end])
rownames(CHR.summary)=unique(CHR)

Position.cum=CHR.summary[,2]-CHR.summary[,1]
CHR.cum=0
for (i in 2:length(unique(CHR))) CHR.cum=c(CHR.cum,sum(Position.cum[1:(i-1)]))
names(CHR.cum)=unique(CHR)

Position.cum=vector()
for (i in unique(CHR)) Position.cum=c(Position.cum,Position[which(CHR==i)]+as.numeric(CHR.cum[which(names(CHR.cum)==i)])-1)

CHR.mid=CHR.cum+(c(CHR.cum[-1],max(CHR.cum)+Position[length(Position)])-CHR.cum)/2

CHR=as.numeric(as.factor(CHR))
GWAlpha=as.numeric(as.character(OUT$Alpha))

if (RAW) {
	png(file=gsub("_out.csv",".png",args[1]),width=1200,height=400)
	plot(0,pch='',xlim=c(0,max(Position.cum)),xlab="Chromosome",xaxt='n',ylab=expression(alpha^2),ylim=c(0,max(GWAlpha^2,na.rm=T)),bty="n",las=2,main=Pheno.name)
	axis(side=1,at=CHR.mid[!is.na(names(CHR.mid))],labels=CHR.labs,pos=0)
	points(GWAlpha[CHR %% 2 !=0]^2~Position.cum[CHR %% 2 !=0],pch=20, col="#104E8B")
	points(GWAlpha[CHR %% 2 ==0]^2~Position.cum[CHR %% 2 ==0],pch=20, col="#ADD8E6")
	#abline(v=CHR.cum[!is.na(names(CHR.cum))][-1],col=2,lty=2,ylim=c(0,max(GWAlpha^2,na.rm=T)))
	if (length(GWAlpha)>100)	abline(h=quantile(GWAlpha^2,1-(100/length(GWAlpha))),col=3)
	dev.off()
} else {
	EST=optim(par=c(mu=0,sd=sd(GWAlpha)),fn=LL,y=GWAlpha,control=list(fnscale=-1,reltol=1e-10))$par
	pval=ifelse(pnorm(GWAlpha,EST[1],EST[2]) > 0.5, 2*(1-pnorm(GWAlpha,EST[1],EST[2])), 2*pnorm(GWAlpha,EST[1],EST[2]))
	p_score=-log10(pval)
	png(file=gsub("_out.csv",".png",args[1]),width=1200,height=400)
	plot(0,pch='',xlim=c(0,max(Position.cum)),xlab="Chromosome",xaxt='n',ylab=expression(-log[10](italic(p))),ylim=c(0,max(p_score,na.rm=T)),bty="n",las=2,main=Pheno.name)
	axis(side=1,at=CHR.mid[!is.na(names(CHR.mid))],labels=CHR.labs,pos=0)
	points(p_score[CHR %% 2 !=0]~Position.cum[CHR %% 2 !=0],pch=20, col="#104E8B")
	points(p_score[CHR %% 2 ==0]~Position.cum[CHR %% 2 ==0],pch=20, col="#ADD8E6")
	#abline(v=CHR.cum[!is.na(names(CHR.cum))][-1],col=2,lty=2)
	if (length(GWAlpha)>100)	abline(h=quantile(p_score,1-(100/length(GWAlpha))),col=3)
	dev.off()
	if (QQ) {
		png(file=gsub("_out.csv","_qqplot.png",args[1]),width=1200,height=400)
		p_order=-log10(pval[order(pval,decreasing=T)])
		p_exp=-log10(seq(1,0,length=length(pval)))
		plot(p_order~p_exp,ylab="Observed p-value",xlab="Expected p-value")
		abline(0,1,col=2)
		dev.off()
	}
}


