#replicate_EP_ST_ML

#index step:

# 1. find all possible events in reference transcriptome
# 2. get the value of PSI from salmon output (example for 1M is displeyed)
# 3. run EventPointer ST statistic (example for 1M is displeyed)
# 4. run EventPointer ML statistic (example for 1M is displeyed)

## 1. find all possible events in the reference transcriptome ######

## this step is common for both EP_ML and EP_ST:

library(EventPointer)

#the same used by SUPPA2
inputGTF <- "./simulation_data/annotation/refseq_hg19.formatted.gtf"
transcriptome <- "refseq_19_formatted"
Pathtxt <- "./simulation_data/EP_events/"

EventsxTran_refseq_19_formatted <- EventDetection_transcriptome(inputFile = inputGTF,
                                                                Transcriptome = transcriptome,
                                                                Pathtxt = Pathtxt,
                                                                cores = 16)

save(EventsxTran_refseq_19_formatted,file="./simulation_data/EP_events/EventsxTran_refseq_19_formatted.RData")

dim(EventsxTran_refseq_19_formatted$ExTP1)
dim(EventsxTran_refseq_19_formatted$ExTP2)



## 2. get the value of PSI from salmon output (example for 1M is displeyed) ######


read_1M <- "1M" 

PathSamples_1M <- dir(paste0("./simulated_data/salmon_output/",read_1M),full.names = TRUE,pattern = "SRR*")
RNA_seq_data_1M <- getbootstrapdata(PathSamples = PathSamples_1M,type = "salmon")
rownames(RNA_seq_data_1M[[1]]) <- gsub("hg19_refGene_","",rownames(RNA_seq_data_1M[[1]]))


miPSI_1M <- GetPSI_FromTranRef(Samples = RNA_seq_data_1M,
                               PathsxTranscript = EventsxTran_refseq_19_formatted,
                               Bootstrap = TRUE,
                               Filter = FALSE,
                               Qn = QN_ALL)

PSI_1M <- miPSI_1M$PSI
dim(PSI_1M)
## is an array with the boostrap. in the main manuscript is depicted with PSI**. One of the matrices that builds this array is the PSI_ML




## 3. run EventPointer ST statistic (example for 1M is displeyed)  ######

ML_1M <- sapply(RNA_seq_data_1M,function(X) X[,1])
yyx <- which(EventsxTran_refseq_19_formatted$transcritnames %in% rownames(ML_1M))
ppx <- match(EventsxTran_refseq_19_formatted$transcritnames[yyx],rownames(ML_1M))
any(is.na(ppx))
ref_1M <- as.matrix(EventsxTran_refseq_19_formatted$ExTPRef[,yyx] %*% ML_1M[ppx,])
mirowmeans_1M <- rowMeans2(ref_1M)
names(mirowmeans_1M) <- rownames(ref_1M)

samplename <- colnames(PSI_1M[1:5,1,])

######## EP_ST
ResultBootstrap_1M <- EventPointer_Bootstraps(PSI = PSI_1M,
                                              Design = Dmatrix,
                                              Contrast = Cmatrix,
                                              cores = 16,
                                              ram = 20,
                                              nBootstraps = 1000,
                                              UsePseudoAligBootstrap = 1,
                                              Threshold = 0)

###### EP_ST with threshold
Result_WITH_BOOTSTRAPS_THRESHOLD_1M <- EventPointer_Bootstraps(PSI = PSI_1M,
                                                               Design = Dmatrix,
                                                               Contrast = Cmatrix,
                                                               cores = 16,
                                                               ram = 20,
                                                               nBootstraps = 1000,
                                                               UsePseudoAligBootstrap = 1,
                                                               Threshold = 0.1)


## 4. run EventPointer ML statistic (example for 1M is displeyed)  ######
###### EP_ML
Result_NO_BOOTSTRAPS_1M <- EventPointer_Bootstraps(PSI = PSI_1M,
                                                   Design = Dmatrix,
                                                   Contrast = Cmatrix,
                                                   cores = 16,
                                                   ram = 20,
                                                   nBootstraps = 1000,
                                                   UsePseudoAligBootstrap = 0,
                                                   Threshold = 0)

###### EP_ML with threshold
Result_NO_BOOTSTRAPS_THRESHOLD_1M <- EventPointer_Bootstraps(PSI = PSI_1M,
                                                             Design = Dmatrix,
                                                             Contrast = Cmatrix,
                                                             cores = 16,
                                                             ram = 20,
                                                             nBootstraps = 1000,
                                                             UsePseudoAligBootstrap = 0,
                                                             Threshold = 0.1)



















































