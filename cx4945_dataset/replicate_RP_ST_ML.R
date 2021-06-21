# replicate EP_ST and EP_ML

# index of steps:

#1. detect all possible events in reference transcriptome
#2. get psi of all the events
#3. run statistic (both the ST and ML version)



######1. detect all possible events in reference transcriptome ######

library(EventPointer)

Eventsxisof_grch37_v74 <- EventDetection_transcriptome(inputFile = "./cx4945/annotation/Homo_sapiens.GRCh37.74.gtf",Transcriptome = "grch37_v74",
                                                       Pathtxt = "./cx4945/EP_events/",cores = 16)
save(Eventsxisof_grch37_v74,file="./cx4945/EP_events/Eventsxisof_grch37_v74.RData")



######2. get psi of all the events ######

load("./cx4945/EP_events/Eventsxisof_grch37_v74.RData")
Table_EventPointer <- read.table("./cx4945/EP_events/EventsFound_grch37_v74.txt",header = TRUE,stringsAsFactors = FALSE,sep = "\t",quote = "")
rownames(Table_EventPointer) <- paste0(Table_EventPointer$GeneName,"_",Table_EventPointer$EventNumber)


pathfiles <- dir("./cx4945/kallisto_output/",full.names = TRUE)
RNA_seq <- getbootstrapdata(PathSamples = pathfiles,type = "kallisto")

Eventsxisof_grch37_v74$transcritnames[1:5]
rownames(RNA_seq[[1]])[1:5]


PSI_evexp <- GetPSI_FromTranRef(Samples = RNA_seq,
                                PathsxTranscript = Eventsxisof_grch37_v74,
                                Bootstrap = T,
                                Filter = T,
                                Qn = 0.001)
PSI <- PSI_evexp$PSI




######3. run statistic (both the ST and ML version) ######

infosamples <- read.table(file = "./cx4945/SraRunTable.txt",stringsAsFactors = F,sep = ",",header = T)
# View(infosamples)
infosamples <- infosamples[match(colnames(PSI[,1,]),infosamples$Run),]

infosamples$source_name
infosamples$agent
infosamples$cell_line

##each cell type
celltype <- factor(infosamples$cell_line,levels = c("SUM-149","MDA-MB-231","MDA-MB-468"))
Dmatrix <- model.matrix(~celltype)
colnames(Dmatrix) <- levels(celltype)
rownames(Dmatrix) <- paste0(infosamples$cell_line,"_",infosamples$agent)
Dmatrix <- cbind(Dmatrix,(infosamples$agent=="CX4945")+0)
colnames(Dmatrix)[4] <- "Treatment_CX4945"

Cmatrix <- matrix(c(0,0,0,1),ncol=1)



## EP_ST 
Fit_Bootstrap <- EventPointer_Bootstraps_print(PSI = PSI,
                                               Design = Dmatrix,
                                               Contrast = Cmatrix,
                                               cores = 6,
                                               ram = 50,
                                               nBootstraps = 1000,
                                               UsePseudoAligBootstrap = 1,
                                               Threshold = 0)


## EP_ML
Fit_Bootstrap_nb <- EventPointer_Bootstraps_print(PSI = PSI,
                                                  Design = Dmatrix,
                                                  Contrast = Cmatrix,
                                                  cores = 6,
                                                  ram = 50,
                                                  nBootstraps = 1000,
                                                  UsePseudoAligBootstrap = 0,
                                                  Threshold = 0)




