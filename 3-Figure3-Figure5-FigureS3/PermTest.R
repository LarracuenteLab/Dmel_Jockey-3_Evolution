library(regioneR)

set.seed(12345) 

##########################################
########### G3/Jockey3 #####################
##########################################

########### Long reads
data=read.table("g2_LongReadWindow_counts.txt", header=T)
data=data[which(data$Ychrom=="other"),]

### Define regions
Genome=toGRanges(data)
Centromere=toGRanges(data[which(data$chromatin=="Centromere"),])
NoHeterochromatin=toGRanges(data[-which(data$chromatin=="Heterochromatin"),])
NoEuchromatin=toGRanges(data[-which(data$chromatin=="Euchromatin"),])

### Jockey3 insertions
data2=data[which(data$G2.count>0),] # remove 0 count
data2=data.frame(data2[rep(seq_len(dim(data2)[1]), data2$G2.count), 1:4, drop = FALSE], row.names=NULL) # one raw per count
Insertions=toGRanges(data2)

### Permutation test
permTest(A=Insertions, B=Centromere, universe=Genome, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoEuchromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoHeterochromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)

########### Short reads
data=read.table("g2_ShortReadWindow_counts.txt", header=T)
data=data[which(data$Ychrom=="other"),]

### Define regions
Genome=toGRanges(data)
Centromere=toGRanges(data[which(data$chromatin=="Centromere"),])
NoHeterochromatin=toGRanges(data[-which(data$chromatin=="Heterochromatin"),])
NoEuchromatin=toGRanges(data[-which(data$chromatin=="Euchromatin"),])

### Jockey3 insertions
data2=data[which(data$G2.count>0),] # remove 0 count
data2=data.frame(data2[rep(seq_len(dim(data2)[1]), data2$G2.count), 1:4, drop = FALSE], row.names=NULL) # one raw per count
Insertions=toGRanges(data2)

### Permutation test
permTest(A=Insertions, B=Centromere, universe=Genome, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoEuchromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoHeterochromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)


##########################################
########### R1 #####################
##########################################

########### Long reads
data=read.table("r1_shortReadWindow_counts.txt", header=T)
data=data[which(data$Ychrom=="other"),]

### Define regions
Genome=toGRanges(data)
Centromere=toGRanges(data[which(data$chromatin=="Centromere"),])
NoHeterochromatin=toGRanges(data[-which(data$chromatin=="Heterochromatin"),])
NoEuchromatin=toGRanges(data[-which(data$chromatin=="Euchromatin"),])

### Jockey3 insertions
data2=data[which(data$G2.count>0),] # remove 0 count
data2=data.frame(data2[rep(seq_len(dim(data2)[1]), data2$G2.count), 1:4, drop = FALSE], row.names=NULL) # one raw per count
Insertions=toGRanges(data2)

### Permutation test
permTest(A=Insertions, B=Centromere, universe=Genome, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoEuchromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)
permTest(A=Insertions, B=Centromere, universe=NoHeterochromatin, ntimes=10000, randomize.function=resampleRegions, evaluate.function=numOverlaps, force.parallel = FALSE)

