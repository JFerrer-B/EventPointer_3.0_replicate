#replicate_EP_BAM

### index steps:

#1. prepare bam files
#2. find all events
#3. compute statistic


#1. prepare bam files ####
pathBAM <- "./HVS_dataset/bam_files"
nameBAM <- dir(pathBAM,pattern = "*Coord.out.bam$")
ref <- "GTF"
pathref <- "./HVS_dataset/annotation/Homo_sapiens.GRCh37.65.gtf"
cores <- 1


SG_RNASeq_rMATS_data <- PrepareBam_EP(Samples=nameBAM,
                                      SamplePath=pathBAM,
                                      Ref_Transc=ref,
                                      fileTransc=pathref,
                                      cores=cores)

# save(SG_RNASeq_rMATS_data,file="./HVS_dataset/ep_bam/SG_RNASeq_rMATS_data.RData")

#2. find all events ####
AllEvents_RNASeq_rMATS_data <- EventDetection(SG_RNASeq_rMATS_data,cores=16,Path="../HVS_dataset/ep_bam/")
# save(AllEvents_RNASeq_rMATS_data,file = "./HVS_dataset/ep_bam/AllEvents_RNASeq_rMATS_data.RData")


#3. compute statistic ####
infosamples <- read.table(file = "./HVS_dataset/filereport_read_run_PRJNA172213_tsv.txt",stringsAsFactors = F,sep = "\t",header = T)
# View(infosamples)
samples_names <- gsub("Alig.*","",colnames(AllEvents_RNASeq_rMATS_data[[1]][[1]]$FPKM))
infosamples <- infosamples[match(samples_names,infosamples$Run),]
Dmatrix <- cbind(1,rep(c(0,1),each=3))
rownames(Dmatrix) <- infosamples$run_accession
colnames(AllEvents_RNASeq_rMATS_data[[1]][[1]]$FPKM)
gsub("Alig.*","",colnames(AllEvents_RNASeq_rMATS_data[[1]][[1]]$FPKM))

Cmatrix <- matrix(c(0,1),ncol=1)

EP_logFC_rMATS_data <- EventPointer_RNASeq(Events = AllEvents_RNASeq_rMATS_data,Design = Dmatrix,Contrast = Cmatrix,Statistic = "LogFC",PSI = TRUE)
# save(EP_logFC_rMATS_data,file="./HVS_dataset/ep_bam/EP_logFC_rMATS_data.RData")



