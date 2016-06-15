#!usr/bin/R

args=commandArgs(TRUE)

Int.size=5000
SNPinKB=1*Int.size/1000

Pheno.name=args[1]
Pheno.name=gsub("GWAlpha_","",Pheno.name)
Pheno.name=gsub("_out.csv","",Pheno.name)

GWAlpha=read.csv(paste("GWAlpha_",Pheno.name,"_out.csv",sep=""),header=T)
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

GWAlpha.order=GWAlpha[order(GWAlpha)]
GWAlpha.QQ=predict(smooth.spline(GWAlpha.order,spar=.6),deriv=2)$y
Thresh=GWAlpha.order[which(GWAlpha.QQ==max(GWAlpha.QQ))]

png(file=paste("GWAlpha_",Pheno.name,".png",sep=""),width=1200,height=600)
plot(0,pch='',xlim=c(0,max(Position.cum)),xlab="Chromosome",xaxt='n',ylab="Alpha^2",ylim=c(0,max(GWAlpha^2,na.rm=T)), main=Pheno.name)
points(GWAlpha[CHR %% 2 !=0]^2~Position.cum[CHR %% 2 !=0],pch=20, col="#104E8B")
points(GWAlpha[CHR %% 2 ==0]^2~Position.cum[CHR %% 2 ==0],pch=20, col="#ADD8E6")
abline(v=CHR.cum[!is.na(names(CHR.cum))][-1],col=2,lty=2)
abline(h=quantile(GWAlpha^2,1-(0.05/length(GWAlpha))),col=3)
axis(side=1,at=CHR.mid[!is.na(names(CHR.mid))],labels=CHR.labs,las=2)
dev.off()

