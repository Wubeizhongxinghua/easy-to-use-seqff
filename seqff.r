##################################################################################################################################################
##################################################################################################################################################
##### SeqFF.R
##### Sung Kim, Phd
##### Sequenom, Inc.
#####
##### Revised by Li Mingyang
##################################################################################################################################################
##### Instructions
##### Command Line 
##### R CMD BATCH 
##### SeqFF.R
library(argparse)
library(Rsamtools)
parser <- ArgumentParser() #创建参数解析对象

parser$add_argument("-f","--file.name",help='Input bam file name') #添加参数
parser$add_argument("-o","--output.filename",help='output file name') #添加参数
args <- parser$parse_args()

ff.pred <- function(gc.norm.bc.61927, B, mu, parameter.1, parameter.2){
  gc.norm.bc.61927[is.na(gc.norm.bc.61927)] <- 0
  gc.norm.bc <- gc.norm.bc.61927[grepl('chr[0-9]', names(gc.norm.bc.61927))]
  gc.norm.bc.c <- gc.norm.bc - mu
  y.hat <- matrix(c(1, gc.norm.bc.c), nrow = 1) %*% B
  y.hat.rep <- sum(y.hat, na.rm = T) / sum(gc.norm.bc)
  ff.hat <- (y.hat.rep + parameter.1) * parameter.2
  return(ff.hat)
}

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##### COMMAND LINE ARGUMENTS

file.name = args$file.name
output.filename = args$output.filename

what <- c('rname','pos')
param <- ScanBamParam(what=what)

load("pd4615-sup-0008-file1.rdata")
bininfo = read.csv("pd4615-sup-0010-table2.csv")
colnames(bininfo)[1]="binName"
bininfo$binorder=c(1:61927)

##### READ IN DATA
print(file.name)
dat <- data.frame(scanBam(file.name,param=param))
colnames(dat)=c("refChr","begin")

print("Data read! Now processing...")

dat=dat[dat$refChr!="*" & dat$refChr!="chrM" ,]
binindex = ifelse((dat$begin%%50000)==0,floor(dat$begin/50000)-1,floor(dat$begin/50000))

fastbincount = table(paste(dat$refChr,binindex, sep="_"))
newtemp=as.data.frame(matrix(NA, ncol=2, nrow=length(fastbincount)))
newtemp[,1]=paste(names(fastbincount))
newtemp[,2]=as.numeric(paste(fastbincount))
colnames(newtemp)=c("binName","counts")


bininfo = merge(bininfo, newtemp, by="binName",all.x=T)
bininfo=bininfo[order(bininfo$binorder),]

##################################################################################################################################################
##### DATA PROCESSING
autosomebinsonly = bininfo$BinFilterFlag==1 & bininfo$CHR!="chrX" & bininfo$CHR!="chrY"
alluseablebins = bininfo$BinFilterFlag==1 
autoscaledtemp  <- bininfo$counts[autosomebinsonly]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
allscaledtemp  <- bininfo$counts[alluseablebins]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
# additive loess correction
mediancountpergc <- tapply(
autoscaledtemp,bininfo$GC[autosomebinsonly], function(x) median(x, na.rm=T))

print('All procession and correction done, now predicting...')

## prediction 
loess.fitted  <- predict( loess(mediancountpergc ~ as.numeric(names(mediancountpergc))), bininfo$GC[alluseablebins]) 
normalizedbincount <- allscaledtemp + ( median(autoscaledtemp, na.rm=T) - loess.fitted )  

bincounts=rep(0,61927)
names(bincounts) = bininfo$binName
bincounts[alluseablebins] <- (normalizedbincount/sum(normalizedbincount, na.rm=T)) * length(normalizedbincount)
bincounts[is.na(bincounts)] <- 0
wrsc=ff.pred(bincounts,B,mu,parameter.1,parameter.2)
enet = matrix(bincounts,nrow=1) %*% matrix(elnetbeta) + elnetintercept
ff=c((wrsc+enet)/2, enet, wrsc)
names(ff)=c("SeqFF","Enet","WRSC")
save(bincounts,file=paste0(output.filename,"_bincounts.RData"))

write.csv(ff, file=output.filename)


