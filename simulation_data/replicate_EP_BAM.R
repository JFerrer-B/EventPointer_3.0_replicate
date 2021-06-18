# replicate_EP_BAM

## index of this step:

# 1. prepare BAM files (example for 1M is displayed)
# 2. detect and classify events (example for 1M is displayed)
# 3. compute EP_BAM statistic (example for 1M is displayed)

# 1. prepare BAM files (example for 1M is displayed) ####

## to store the results for 1M:

command <- "mkdir ./simulation_data/ep_old/1M"
system(command)

nameBAM <- dir("./simulation_data/BAM_rMATS/1M",pattern = "*Coord.out.bam$")
pathBAM <- "./simulation_data/BAM_rMATS/1M"
ref <- "GTF"
pathref <- "./simulation_data/annotation/refseq_hg19.formatted.gtf"
cores <- 12

SG_RNASeq_1M<-PrepareBam_EP(Samples=nameBAM,
                            SamplePath=pathBAM,
                            Ref_Transc=ref,
                            fileTransc=pathref,
                            cores=cores)

## it is recommended to save this data in a .RData file, since the calculation time is long.
# save(SG_RNASeq_1M,file="./suppa2_simulacion/ep_old/1M/SG_RNASeq_1M.RData")

# 2. detect and classify events (example for 1M is displayed) ####

AllEvents_RNASeq_1M<-EventDetection(SG_RNASeq_1M,cores=16,Path="./suppa2_simulacion/ep_old/1M/")

## it is recommended to save this data in a .RData file, since the calculation time is long.
# save(AllEvents_RNASeq_1M,file = "./suppa2_simulacion/ep_old/1M/AllEvents_RNASeq_1M.RData")

# 3. compute EP_BAM statistic (example for 1M is displayed) ####

Dmatrix <- cbind(rep(1,9),
                 rep(c(1,0,0),3),
                 rep(c(0,1,0),3))
# rownames(Dmatrix) <- samplename
Cmatrix <- cbind(
  c(0,1,0),
  c(0,0,1))


EP_logFC_1M <- EventPointer_RNASeq(Events = AllEvents_RNASeq_1M,Design = Dmatrix,Contrast = Cmatrix,Statistic = "LogFC",PSI = TRUE)
## even it doesn't take long time to run the above code, it is recommended the result in a .RData file.
# save(EP_logFC_1M,file="./suppa2_simulacion/ep_old/1M/EP_logFC_1M.RData")



