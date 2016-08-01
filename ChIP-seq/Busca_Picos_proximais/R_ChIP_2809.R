library ("ChIPpeakAnno")
library(dataframes2xls)
library(ggplot2)

##command: R CMD BATCH <args1 = GFF> <args2 = Peaks> <args3 = Output_TAG>
args = commandArgs(trailingOnly = TRUE)

print(args)

#=============================Dados do GFF===========================================================

GFF_Esm_rgData <- GFF2RangedData(args[1], header=FALSE)

#====================================================================================================

PeakData_RD <-  BED2RangedData(args[2], header=FALSE)

#=========================ChIPpeakAnno - busca por genes proximos a regiões enriquecidas================================

#ChIPpeakAnno
AnnoPeak = annotatePeakInBatch(PeakData_RD, AnnotationData=GFF_Esm_rgData)

write.table((as.data.frame(AnnoPeak)) ,file=args[3] , sep="\t")

#Dados estatísticos da distancia do fragmento (menor distancia)
##Summary_shortDistance <- summary(as.data.frame(AnnoPeak)[13])
##write.table((as.data.frame(Summary_shortDistance)[3]) ,file="/media/criph/Analysis/FastQ/ChIP_1_Rep1_jan2015/Analysis/Dataset/Ty_Esm_peaks_MACS14.aPeak_short_distance", sep="\t", row.names=FALSE, col.names=FALSE)


#histograma da distancia da ponta do fragmento ao gene (distance feature)
##data_hist <- as.data.frame(AnnoPeak)
##pdf("/media/criph/Analysis/FastQ/ChIP_1_Rep1_jan2015/Analysis/Dataset/Ty_Esm_peaks_MACS14_graph.pdf", width=8, height = 6)
##ggplot(data_hist, aes(x=data_hist[,12])) + geom_histogram(binwidth=200) + xlab("Frag_Gene_Distance")
##dev.off()
##rm (data_hist)

#Imprime Arquivos com informações:
##Ty_Esm_peaks_MACS14_summary <- summary(as.data.frame(Ty_Esm_peaks_MACS14.RangedData)[4])
##Ty_Esm_peaks_MACS14 <- dim(Ty_Esm_peaks_MACS14.RangedData)[1]
##write.table((as.data.frame(Ty_Esm_peaks_MACS14.RangedData)) ,file="/media/criph/Analysis/FastQ/ChIP_1_Rep1_jan2015/Analysis/Dataset/Ty_Esm_peaks_MACS14.RangedData", sep="\t")
##write.table((as.data.frame(Ty_Esm_peaks_MACS14_summary)[3]) ,file="/media/criph/Analysis/FastQ/ChIP_1_Rep1_jan2015/Analysis/Dataset/Ty_Esm_peaks_MACS14_summary", sep="\t", row.names=FALSE, col.names=FALSE)
##write.table((as.data.frame(Ty_Esm_peaks_MACS14)) ,file="/media/criph/Analysis/FastQ/ChIP_1_Rep1_jan2015/Analysis/Dataset/Ty_Esm_peaks_MACS14", sep="\t")

