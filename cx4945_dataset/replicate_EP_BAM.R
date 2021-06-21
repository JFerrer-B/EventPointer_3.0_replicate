# replicate EP_BAM

### index steps:

#1. prepare bam files
#2. find all events
#3. compute statistic

#1. prepare bam files ####
pathBAM <- "./cx4945_dataset/bam_files/"
nameBAM <- dir(pathBAM,pattern = "*Coord.out.bam$")
ref <- "GTF"
pathref <- "./cx4945_dataset/annotation/Homo_sapiens.GRCh37.74.gtf"
cores <- 1

SG_RNASeq_JP_data <- PrepareBam_EP(Samples=nameBAM,
                                   SamplePath=pathBAM,
                                   Ref_Transc=ref,
                                   fileTransc=pathref,
                                   cores=cores)

save(SG_RNASeq_JP_data,file="./cx4945_dataset/ep_bam/SG_RNASeq_JP_data.RData")


#2. find all events ####
AllEvents_RNASeq_JP_data <- EventDetection(SG_RNASeq_JP_data,cores=16,Path="./cx4945_dataset/ep_bam/")
save(AllEvents_RNASeq_JP_data,file = "./cx4945_dataset/ep_bam/AllEvents_RNASeq_JP_data.RData")


#3. compute statistic ####
infosamples <- read.table(file = "./cx4945_dataset/SraRunTable.txt",stringsAsFactors = F,sep = ",",header = T)

samples_names <- gsub("Alig.*","",colnames(AllEvents_RNASeq_JP_data[[1]][[1]]$FPKM))
infosamples <- infosamples[match(samples_names,infosamples$Run),]


celltype <- factor(infosamples$cell_line,levels = c("SUM-149","MDA-MB-231","MDA-MB-468"))
Dmatrix <- model.matrix(~celltype)
colnames(Dmatrix) <- levels(celltype)
rownames(Dmatrix) <- paste0(infosamples$cell_line,"_",infosamples$agent)
Dmatrix <- cbind(Dmatrix,(infosamples$agent=="CX4945")+0)
colnames(Dmatrix)[4] <- "Treatment_CX4945"

Cmatrix <- matrix(c(0,0,0,1),ncol=1)
Dmatrix


EP_logFC_JP_data <- EventPointer_RNASeq(Events = AllEvents_RNASeq_JP_data,Design = Dmatrix,Contrast = Cmatrix,Statistic = "LogFC",PSI = TRUE)
save(EP_logFC_JP_data,file="./cx4945_dataset/ep_bam/EP_logFC_JP_data.RData")

