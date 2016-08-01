##command: R CMD BATCH <args1 = GFF> <args2 = Peaks> <args3 = Output_TAG>

#options(echo=FALSE)
args = commandArgs(trailingOnly = TRUE)
print(args)


#if (exists(args[1] & args[2] & args[3])) {

	library(ChIPpeakAnno)
	library(dataframes2xls)
	library(ggplot2)


	#=============================Dados do GFF===========================================================

	GFF_Esm_rgData <- GFF2RangedData(args[1], header=FALSE)

	#====================================================================================================

	PeakData_RD <-  BED2RangedData(args[2], header=FALSE)

	#=========================ChIPpeakAnno - busca por genes proximos a regiões enriquecidas================================

	#ChIPpeakAnno
	AnnoPeak = annotatePeakInBatch(PeakData_RD, AnnotationData=GFF_Esm_rgData)

	write.table((as.data.frame(AnnoPeak)) ,file=paste0(args[3],".aPeaks"), sep="\t")

	#Dados estatísticos da distancia do fragmento (menor distancia)
	Summary_shortDistance <- summary(as.data.frame(AnnoPeak)[13])
	write.table((as.data.frame(Summary_shortDistance)[3]) ,file=paste0(args[3], "_ShortDist_summary.out"), sep="\t", row.names=FALSE, col.names=FALSE)


	#histograma da distancia da ponta do fragmento ao gene (distance feature)
	data_hist <- as.data.frame(AnnoPeak)
	pdf(paste0(args[3] , "_Histogram.pdf") , width=8, height = 6)
	ggplot(data_hist, aes(x=data_hist[,13])) + geom_histogram(binwidth=200) + xlab("Frag_Gene_Distance")
	dev.off()

	#graficos de posicionamento
	pdf(paste0(args[3] , "_pos_graph.pdf") , width=8, height = 6)
	ggplot(data_hist, aes(x=data_hist[,12])) + geom_histogram(binwidth=200) + xlab("Position_Gene")
	dev.off()
	rm (data_hist)

	#Imprime Arquivos com informações do Pico enriquecido:
	Peak_summary <- summary(as.data.frame(PeakData_RD)[4])
	NunPeaks <- dim(PeakData_RD)[1]
	write.table((as.data.frame(Peak_summary)[3]) , file=paste0(args[3], "_Peak_summary") , sep="\t", row.names=FALSE, col.names=FALSE)
	write.table((as.data.frame(NunPeaks)) ,file=paste0(args[3], "_NumPeaks"), sep="\t")


#}

#else {print "Falta de argumento"}
