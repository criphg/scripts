#Ref:
#https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html
#https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
#https://jianhong.github.io/trackViewerBiocAsia2020Workshop/articles/trackViewer.html

#Mutações compartilhadas - referencia: https://covariants.org/shared-mutations
setwd("/media/sf_D_DRIVE/B_1_1_529/")


library(trackViewer)
#library(ggplot2)
df_S = read.delim("DF_S_mutation_count.tab")
#features <- GRanges("chr1", IRanges(c(1,501,1001), 
#                                    width=c(120,400,405), 
#                                    names=paste0("block", 1:3)),
#                   fill = c("#FF8833", "#51C6E6","#DFA320"),
#                   height = c(0.02,0.05,0.08))
#GRanges("chr1", IRanges(c(1,501,1001), width = ))

#SNP = c(10,100,105,108,400,410,420,600,700,805, 840, 1400,1402)

#sample.gr = GRanges("chr1", IRanges(SNP, width =1 , names=paste0("snp", SNP)),
#                    color = sample.int(6, length(SNP), replace = TRUE),
#                    score = sample.int(5, length(SNP), replace = TRUE))

#lolliplot(sample.gr, features)


#df_S


#features2 = GRanges("spike", IRanges(c(1),
#                   width=c(1276),
#                   names=c("Spike")), 
#                   fill = c("red"), 
#                   height= c(0.05))
#Spike = GRanges("spike", IRanges(df_S[,1], width = 1, names = df_S[,4]),
#                color = df_S[,3])
#                #score = sample.int(1, length(df_S), replace = TRUE))
#lolliplot(Spike, features2)


#gr = GRanges(
#  seqnames = Rle(c("Spike","teste"), c(1,1)), 
 # ranges = IRanges(c(1,1), width = c(1276,1276), 
#                   names = c("Spike","Teste"), 
#                   fill = c("red","green"), height = (c(0.05, 0.05)))
#)
#gr_mut = GRanges(
#  seqnames = Rle(c("Spike","teste"), c(length(df_S[,1]),1)),
#  ranges = IRanges(c(df_S[,1],100), width = 1, names = c(df_S[,4],"Ateste"),
#                   color = 1:40)
#)
#lolliplot(gr_mut, gr)



#region_dt_S = data.frame("coord" = c("1-1273","22-852","853-1273","319-541"),
#                         "name" = c("Spike","S1","S2","RBD"))


#gr_spike = GRanges(
#  "Spike", 
#  ranges = IRanges(c(1,22,853,319), end = c(1273,852,1273,541), 
#                   names = c("Spike","S1","S2","RDB"), 
#                   fill = c("red","green","grey","blue"), height = (c(0.08,0.05,0.05,0.03))))
  


SNP = df_S[,1]
MUT_name = df_S[,4]

{c_v1 = 2
c_v2 = 7
color_v = df_S[,5]
color_v[color_v == 0] = c_v1
color_v[color_v == 1] = c_v2}
#color_v = df_S[,5]+c_v

{df_S = transform(df_S, Freq=ave(seq(nrow(df_S)), coord, FUN=length))

gr_spike = GRanges(
  "Spike", 
  ranges = IRanges(c(1,22,853,319), end = c(1273,852,1273,541), 
                   names = c("Spike","S1","S2","RDB"), 
                   fill = c("red","grey","grey40","black"), height = (c(0.08,0.05,0.05,0.03))))

Omicron_Mut = GRanges("Spike", IRanges(SNP, width =1 , names = MUT_name),
                       color = color_v,
                       score = ifelse(color_v == c_v1, 6, 15))
Omicron_Mut$label.parameter.rot = 45
Omicron_Mut$label = as.character(1:length(df_S[,1]))
Omicron_Mut$label.col = "white"
Omicron_Mut$border = ifelse(color_v == c_v1, "black", "gray80")
#Omicron_Mut$border.width = 20

legend = list(labels=c("Shared","Exclusive"),col=c(0,0),fill=c(c_v1,c_v2))

xaxis <- c(1, 200, 400, 600, 800 ,1000, 1200, 1273)

lolliplot(Omicron_Mut, gr_spike, legend = legend, ylab = "Ômicron Mutations", type = "pin", xaxis = xaxis, yaxis = FALSE)
}

#yaxis = c(0,5,15)

jpeg(file = "Omicron_mutations_20211130_S.jpg", width = 2684, height = 894, pointsize = 40, units = "px")
lolliplot(Omicron_Mut, gr_spike, legend = legend, 
          ylab = "Ômicron Mutations", 
          type = "pin", xaxis = xaxis, yaxis = FALSE, cex = 0.6)
dev.off()

pdf(file = "Omicron_mutations_20211130_S.pdf", width = 84, height = 44, pointsize = 82)
lolliplot(Omicron_Mut, gr_spike, legend = legend, 
          ylab = "Ômicron Mutations", 
          type = "pin", xaxis = xaxis, yaxis = FALSE, cex = 0.7)
dev.off()

###########################################################################
#####N


df_N = read.delim("DF_N_mutation_count.tab")

SNP = df_N[,1]
MUT_name = df_N[,4]

{c_v1 = 2
  c_v2 = 7
  color_v = df_N[,5]
  color_v[color_v == 0] = c_v1
  color_v[color_v == 1] = c_v2}
#color_v = df_S[,5]+c_v

{df_N = transform(df_N, Freq=ave(seq(nrow(df_N)), coord, FUN=length))
  
  gr_N = GRanges(
    "N", 
    ranges = IRanges(c(1), end = c(419), 
                     names = c("Nucleocapsid"), 
                     fill = c("royalblue3"), height = (c(0.08))))
  
  Omicron_Mut2 = GRanges("N", IRanges(SNP, width =1 , names = MUT_name),
                        color = color_v,
                        #score = ifelse(color_v == c_v1, 6, 15))
                        score = ifelse(color_v == c_v1, 3, 8))
  Omicron_Mut2$label.parameter.rot = 45
  Omicron_Mut2$label = as.character(1:length(df_N[,1]))
  Omicron_Mut2$label.col = "white"
  Omicron_Mut2$border = ifelse(color_v == c_v1, "black", "gray80")
  #Omicron_Mut$border.width = 20
  
  legend = list(labels=c("Shared","Exclusive"),col=c(0,0),fill=c(c_v1,c_v2))
  
  xaxis <- c(1, 200, 400, 419)
  
  lolliplot(Omicron_Mut2, gr_N, legend = legend, ylab = "Ômicron Mutations", type = "pin", xaxis = xaxis, yaxis = FALSE)
}

#yaxis = c(0,5,15)

jpeg(file = "Omicron_mutations_20211130_N.jpg", width = 2684, height = 894, pointsize = 40, units = "px")
lolliplot(Omicron_Mut2, gr_N, legend = legend, 
          ylab = "Ômicron Mutations", 
          type = "pin", xaxis = xaxis, yaxis = FALSE, cex = 0.7)
dev.off()

pdf(file = "Omicron_mutations_20211130_N.pdf", width = 84, height = 44, pointsize = 82)
lolliplot(Omicron_Mut2, gr_N, legend = legend, 
          ylab = "Ômicron Mutations", 
          type = "pin", xaxis = xaxis, yaxis = FALSE, cex = 0.7)
dev.off()


pdf(file = "Omicron_mutations_20211130_N_2.pdf", width = 34, height = 44, pointsize = 82)
lolliplot(Omicron_Mut2, gr_N, legend = legend, 
          ylab = "Ômicron Mutations", 
          type = "pin", xaxis = xaxis, yaxis = FALSE, cex = 0.7)
dev.off()

