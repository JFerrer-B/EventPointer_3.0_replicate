#replicate_EP_ST_ML



# index of steps:

#1. detect all possible events in reference transcriptome
#2. get psi of all the events
#3. run statistic (both the ST and ML version)



library(EventPointer)

######1. detect all possible events in reference transcriptome ######


inputFile <- "./HVS_dataset/annotation/Homo_sapiens.GRCh37.65.gtf"
Transcriptome <- "Ensembl_GRCh37_v65"
Pathtxt <- "./HVS_dataset/EP_events/"
PathGTF <- "./HVS_dataset/EP_events/"
EventXtrans_grch37_v65 <- EventDetection_transcriptome(inputFile = inputFile,
                                                       Transcriptome = Transcriptome,
                                                       Pathtxt = Pathtxt,
                                                       cores = 16)

# save(EventXtrans_grch37_v65,file="./HVS_dataset/EP_events/EventXtrans_grch37_v65.RData")




load("./HVS_dataset/EP_events/EventXtrans_grch37_v65.RData")
EventsInfo <- read.delim(file="./HVS_dataset/EP_events/EventsFound_Ensembl_GRCh37_v65.txt",
                         header = TRUE,stringsAsFactors = FALSE)

EventsInfo$ID <- paste0(EventsInfo$GeneName,"_",EventsInfo$EventNumber)

######2. get psi of all the events ######

PathToSamples <- "./HVS_dataset/kallisto_output"
PathToSamples <- dir(PathToSamples,full.names = TRUE)
type="kallisto"
RNADATA_Kallisto <- getbootstrapdata(PathSamples = PathToSamples,type = type)

PSI_EXPRESION_K <- GetPSI_FromTranRef(Samples = RNADATA_Kallisto,PathsxTranscript = EventXtrans_grch37_v65,Bootstrap = TRUE,Filter = FALSE,Qn=0.25)
PSI_KALLISTO <- PSI_EXPRESION_K$PSI


######3. run statistic (both the ST and ML version) ######


infosamples <- read.table(file = "./HVS_dataset/filereport_read_run_PRJNA172213_tsv.txt",stringsAsFactors = F,sep = "\t",header = T)
# View(infosamples)
dim(infosamples)
infosamples <- infosamples[,1:3]
Dmatrix <- cbind(1,rep(c(0,1),each=3))
rownames(Dmatrix) <- infosamples$run_accession
Dmatrix
colnames(ML_isoform_expression)
Cmatrix <- matrix(c(0,1),ncol=1)

ResultBootstrap_kallisto_rmats_data <- EventPointer_Bootstraps(PSI = PSI_KALLISTO,
                                                               Design = Dmatrix,
                                                               Contrast = Cmatrix,
                                                               cores = 12,
                                                               ram = 20,
                                                               nBootstraps = 1000,
                                                               UsePseudoAligBootstrap = 1)




ResultBootstrap_kallisto_rmats_data_ML <- EventPointer_Bootstraps(PSI = PSI_KALLISTO,
                                                                  Design = Dmatrix,
                                                                  Contrast = Cmatrix,
                                                                  cores = 12,
                                                                  ram = 20,
                                                                  nBootstraps = 1000,
                                                                  UsePseudoAligBootstrap = 0)




