## - chamada de argumentos
args <- commandArgs(trailingOnly = TRUE)
print(args)

#args[1] = GFF
#args[2] = tabDiffBind
#args[3] = tag file output
#args[4] = summit (distancia esperado para o pico de cada fragmento)

print(paste0("ARG1=",args[1]))
print(paste0("ARG2=",args[2]))
print(paste0("ARG3=",args[3]))
print(paste0("ARG4=",args[4]))

cat("\nDiffBind analysis...\n")

#if (exists(args[1] & args[2] & args[3] & args[4]) {
#1 - chamar a biblioteca do DiffBind

library("DiffBind")
library(dataframes2xls)
library(ggplot2)

dir_D = "DiffBind_D_DESEQ2/" 
dir.create(dir_D)

#2 - Carregar as tabelas para do DiffBind de Narrow e Broad(Opcional)
samples <- read.csv(args[2])
write.table((as.data.frame(samples)) ,file=paste0(dir_D,args[3],"_samples.Tab"), sep="\t")

#3 - leitura dos picos - verificar como imprimir o plot como pdf.
AnA <- dba(sampleSheet=args[2])

#4 - Contar as reads - veririfar como estabelecer um summit no momento da chamada do script
AnA <- dba.count(AnA, summits=as.numeric(args[4]))

#5 - Estabelecendo contraste - deve-se usar minMembers para estabelecer contrates de 2 replicatas (default 3 ou mais).
AnA <- dba.contrast(AnA, categories=DBA_CONDITION, minMembers=2)

#6 - Análise diferencial
AnA <- dba.analyze(AnA)
#Enviando dados do objeto dba para um arquivo (importante para se ter dados de análise Dseq e Edger).
sink(file=paste0(dir_D,args[3],"_DB_analisys.txt"))
AnA
sink()

#grafico mostrando a subdivisão amostral
pdf(paste0(dir_D , args[3] , "_contrastPlot.pdf") , width=8, height = 6)
plot(AnA, contrast=1)
dev.off()

#7 - Recuperação do sítios de interação diferenciais.

#AnA.DB <- dba.report(AnA)


AnA.DB <- dba.report(AnA)
write.table((as.data.frame(AnA.DB)) ,file=paste0( dir_D, args[3] , "_DB_analisys.out"), sep="\t")


#8 - MAplot
print("MA_Plot ...")
pdf(paste0(dir_D, args[3] , "_MAplot.pdf") , width=8, height = 6)
dba.plotMA(AnA)
dev.off()

#9 - BoxPlot

pdf(paste0(dir_D, args[3] , "_Boxplot.pdf") , width=8, height = 6)
BPlotN = dba.plotBox(AnA)
dev.off()

## movendo arquivos criados para diretório DiffBind
#system("mkdir -p DiffBind_D_DESEQ2")
#system("find . -maxdepth 1 -type f  -not -wholename '*tabdiff_*' -exec mv '{}' './DiffBind_D' ';' ")
file.rename(from = "Rplots.pdf", to = paste0(dir_D,"Rplots.pdf"))

#10 - ChIPpeakAnno - inicio

cat("\nChiIPpeakAnno analysis...\n")

        library(ChIPpeakAnno)

	dir_C = "ChIPpeakAnno_D_DESEQ2/"
	dir.create(dir_C)

        #=============================Dados do GFF===========================================================

        GFF_Esm_rgData <- GFF2RangedData(args[1], header=FALSE)

        #====================================================================================================

	###PARA NARROW

        PeakData_RD <-  AnA.DB

        #=========================ChIPpeakAnno - busca por genes proximos a regiões enriquecidas================================

        #ChIPpeakAnno
        AnnoPeak = annotatePeakInBatch(PeakData_RD, AnnotationData=GFF_Esm_rgData)

        write.table((as.data.frame(AnnoPeak)) ,file=paste0(dir_C, args[3],".out"), sep="\t")

        #Dados estatísticos da distancia do fragmento (menor distancia)
        Summary_shortDistance <- summary(as.data.frame(AnnoPeak)[13])
        write.table((as.data.frame(Summary_shortDistance)[3]) ,file=paste0(dir_C,args[3], "_ShortDist_summary.out"), sep="\t", row.names=FALSE, col.names=FALSE)


        #histograma da distancia da ponta do fragmento ao gene (distance feature)
        data_hist <- as.data.frame(AnnoPeak)
        pdf(paste0(dir_C,args[3] , "_Histogram.pdf") , width=8, height = 6)
        ggplot(data_hist, aes(x=data_hist[,18])) + geom_histogram(binwidth=200) + xlab("Frag_Gene_Distance")
        dev.off()

        #graficos de posicionamento
        pdf(paste0(dir_C,args[3] , "_pos_graph.pdf") , width=8, height = 6)
        ggplot(data_hist, aes(x=data_hist[,19])) + geom_histogram(binwidth=200) + xlab("Position_Gene")
        dev.off()
        rm (data_hist)

        #Imprime Arquivos com informações do Pico enriquecido:
        Peak_summary <- summary(as.data.frame(PeakData_RD)[4])
        NunPeaks <- dim(PeakData_RD)[1]
        write.table((as.data.frame(Peak_summary)[3]) , file=paste0(dir_C, args[3], "_Peak_summary") , sep="\t", row.names=FALSE, col.names=FALSE)
        write.table((as.data.frame(NunPeaks)) ,file=paste0(dir_C, args[3], "_NumPeaks"), sep="\t")

	#movendo arquivos para diretório especifico ChIPpeakAnno
#	system("mkdir -p ChIPpeakAnno_D_DESEQ2")
#	system("find . -maxdepth 1 -type f  -not -wholename '*tabdiff_*' -exec mv '{}' './ChIPpeakAnno_D' ';' ")

#}

#else {print "Falta de argumento"}


