a = "/media/criph/STO2/Circlize/"
setwd(a)

library(circlize)

##############TABELAS - DADOS #######################################
#####################################################################
#coordenadas dos chromossomos
x_c = read.table("data/chrm_size_Lbr.tab", sep = "\t")
x_c2 = data.frame(V1 = x_c$V1, V2 = rep(0,36))
x_c3= rbind(x_c,x_c2)

#coordenadas dos CDS
tab_cds = read.table(file = "data/Tabela_geneLocation_geneID.tab", sep = "\t")

#removendo os scafolds
id_exclude_cds = !grepl("_SCAF", tab_cds$V1)
tab_cds = tab_cds[id_exclude_cds,]
#removendo os scaffolds do factor.
fact_1 = factor(as.character(tab_cds$V1))
levels(fact_1)
tab_cds$V1 = fact_1 


#coordenads dos genes DE - PRO_META.
# tab_DEG_PM = read.csv("data/PRO_META_FCp.csv", header = F)
# tab_DEG_MA = read.csv("data/META_AMA_FCp.csv", header = F)
# tab_DEG_AP = read.csv("data/AMA_PRO_FCp.csv", header = F)
# #recuperando coordenadas das DEG
# tab_cood_all = read.table("data/Tabela_geneLocation_geneID.tab", sep = "\t")
# #Função para recuperar as coordenadas dos CDS e eliminar os Scaffolds
# coord_recover_Chrm = function(tab_DEG, tab_cood_all){
#   index = match(tab_DEG[,1], tab_cood_all[,5])
#   tab_cood_DEG = cbind(tab_DEG, tab_cood_all[index,])
#   #removendo os scafolds
#   id_exclude = !grepl("_SCAF", tab_cood_DEG[,4])
#   tab_chr_DEG = tab_cood_DEG[id_exclude,]
#   #nomeando colunas para evitar dois V1
#   colnames(tab_chr_DEG) = paste0("c", seq(1,dim(tab_chr_DEG)[2]))
#   #removendo os scaffolds do factor.
#   fact = factor(as.character(tab_chr_DEG[,4]))
#   levels(fact)
#   tab_chr_DEG[,4] = fact
#   return(tab_chr_DEG)
# }  

tab = read.csv("data/20180417_DEncRNAs_FromTable20170816.csv", header = T)
tab$PROvsMETA.UP = as.numeric(gsub("-",0,tab$PROvsMETA.UP))
tab$PROvsMETA.DOWN = as.numeric(gsub("-",0,tab$PROvsMETA.DOWN))
tab$METAvsAMA.UP = as.numeric(gsub("-",0,tab$METAvsAMA.UP))
tab$METAvsAMA.DOWN = as.numeric(gsub("-",0,tab$METAvsAMA.DOWN))
tab$AMAvsPRO.UP = as.numeric(gsub("-",0,tab$AMAvsPRO.UP))
tab$AMAvsPRO.DOWN = as.numeric(gsub("-",0,tab$AMAvsPRO.DOWN))

#Divisão de 1 pelas colunas DOWM para obter valor corredo de FC "negativo"
tab[,14] = ifelse(tab[,14] != 0, 1/tab[,14], 0)
tab[,16] = ifelse(tab[,16] != 0, 1/tab[,16], 0)
tab[,18] = ifelse(tab[,18] != 0, 1/tab[,18], 0)

#Acrescimo de ultima coluna com soma de FC UP+DW (PM MA AP)
tab = cbind(tab, FC_PM = rowSums(tab[,c(13,14)]))
tab = cbind(tab, FC_MA = rowSums(tab[,c(15,16)]))
tab = cbind(tab, FC_AP = rowSums(tab[,c(17,18)]))

#Obtenção dos valores de logFC
tab = cbind(tab, Log_FC_PM = ifelse(tab$FC_PM != 0, log(tab$FC_PM, 2), 0))
tab = cbind(tab, Log_FC_MA = ifelse(tab$FC_MA != 0, log(tab$FC_MA, 2), 0))
tab = cbind(tab, Log_FC_AP = ifelse(tab$FC_AP != 0, log(tab$FC_AP, 2), 0))

#View(tab_ncRNADE_Pred2_ALL_15[grepl("_SCAF",tab_ncRNADE_Pred2_ALL_15[,2]),])

#fusão da tabela com classes
tab_class = read.csv("data/ncRNA_characteristics_result--after_Pfam_V2.csv", skip = 1, header = T)
idx_tab_class = match(tab$ncRNA.ID, tab_class$ncRNA.ID)
tab = cbind(tab, tab_class[idx_tab_class, 7:ncol(tab_class)])

# #ncRNA DE - FC 1.5 - PRED >= 2
# tab_DEG_PM = tab[((tab$Log_FC_PM >= log(1.5, 2)) | (tab$Log_FC_PM <= log(1/1.5, 2))) & (tab$Prediction.Score >= 2),]
# tab_DEG_MA = tab[((tab$Log_FC_MA >= log(1.5, 2)) | (tab$Log_FC_MA <= log(1/1.5, 2))) & (tab$Prediction.Score >= 2),]
# tab_DEG_AP = tab[((tab$Log_FC_AP >= log(1.5, 2)) | (tab$Log_FC_AP <= log(1/1.5, 2))) & (tab$Prediction.Score >= 2),]

#ncRNA DE - FC 1.5 (não foi feito filtro de preditores)
tab_DEG_PM = tab[((tab$Log_FC_PM >= log(1.5, 2)) | (tab$Log_FC_PM <= log(1/1.5, 2))),]
tab_DEG_MA = tab[((tab$Log_FC_MA >= log(1.5, 2)) | (tab$Log_FC_MA <= log(1/1.5, 2))),]
tab_DEG_AP = tab[((tab$Log_FC_AP >= log(1.5, 2)) | (tab$Log_FC_AP <= log(1/1.5, 2))),]



# ####GRUPOS  - FC 1.5 - PRED >= 2 (não usado!)
# 
# tab_ncRNADE_Pred2_ALL_15  = tab[((tab$Log_FC_PM >= log(1.5, 2)) | (tab$Log_FC_PM <= log(1/1.5, 2))) & ((tab$Log_FC_MA >= log(1.5, 2)) | (tab$Log_FC_MA <= log(1/1.5, 2))) & ((tab$Log_FC_AP >= log(1.5, 2)) | (tab$Log_FC_AP <= log(1/1.5, 2))) & (tab$Prediction.Score >= 2),]
# tab_ncRNADE_ALL_15  = tab[((tab$Log_FC_PM >= log(1.5, 2)) | (tab$Log_FC_PM <= log(1/1.5, 2))) & ((tab$Log_FC_MA >= log(1.5, 2)) | (tab$Log_FC_MA <= log(1/1.5, 2))) & ((tab$Log_FC_AP >= log(1.5, 2)) | (tab$Log_FC_AP <= log(1/1.5, 2))) ,]


#Função para eliminar registros encontrados em Scaffolds
remove_scaf = function(tab_DEG){
  #removendo os scafolds
  id_include = !grepl("_SCAF", tab_DEG[,2])
  tab_chr_DEG = tab_DEG[id_include,]
  #removendo os scaffolds do factor.
  fact = factor(as.character(tab_chr_DEG[,2]))
  levels(fact)
  tab_chr_DEG[,2] = fact
  return(tab_chr_DEG)
}
#tabelas dos DEG não scaffolds com suas coordenadas (usando função acima)
tab_chr_DEG_PM = remove_scaf(tab_DEG_PM)
tab_chr_DEG_MA = remove_scaf(tab_DEG_MA)
tab_chr_DEG_AP = remove_scaf(tab_DEG_AP)

# ###FALTA SEPARAR OS GRUPOS - não usado
# tab_Groups = remove_scaf(tab_ncRNADE_Pred2_ALL_15)
# table(grepl("_SCAF",tab_ncRNADE_Pred2_ALL_15[,2]))
#tab_ncRNADE_Pred2_least1_15  = tab[(((tab$PROvsMETA.UP >= 1.5) | (tab$PROvsMETA.DOWN >= 1.5)) | ((tab$METAvsAMA.UP > 2) | (tab$METAvsAMA.DOWN >= 1.5)) | ((tab$AMAvsPRO.UP >= 1.5) | (tab$AMAvsPRO.DOWN >= 1.5))) & (tab$Prediction.Score >= 2),]
### A SER FEITO!!!! COM OS ARQUIVOS OBTIDOS DO PROGRAMA DE SEPARACAO DE GRUPOS


####Ja feito acima para a tabela completa.
# #Divisão de 1 pelas colunas DOWM para obter valor corredo de FC "negativo"
# tab_chr_DEG_PM[,14] = ifelse(tab_chr_DEG_PM[,14] != 0, 1/tab_chr_DEG_PM[,14], 0)
# tab_chr_DEG_MA[,16] = ifelse(tab_chr_DEG_MA[,16] != 0, 1/tab_chr_DEG_MA[,16], 0)
# tab_chr_DEG_AP[,18] = ifelse(tab_chr_DEG_AP[,18] != 0, 1/tab_chr_DEG_AP[,18], 0)
# 
# #Acrescimo de ultima coluna com soma de FC UP+DW (PM MA AP)
# tab_chr_DEG_PM = cbind(tab_chr_DEG_PM, FC_PM = rowSums(tab_chr_DEG_PM[,c(13,14)]))
# tab_chr_DEG_MA = cbind(tab_chr_DEG_MA, FC_MA = rowSums(tab_chr_DEG_MA[,c(15,16)]))
# tab_chr_DEG_AP = cbind(tab_chr_DEG_AP, FC_AP = rowSums(tab_chr_DEG_AP[,c(17,18)]))
# 
# #View(tab_ncRNADE_Pred2_ALL_15[grepl("_SCAF",tab_ncRNADE_Pred2_ALL_15[,2]),])

#bed_DEG - funcao para criar o bed e acrescentar codigo de location -> 1-5UTR 2-Undetermined 3-3UTR 4-CDS
bed_func_DEG = function(tab_chr_DEG_PM, col) {
  bed_PM = tab_chr_DEG_PM[,c(2:5,col,7,6)]
  #Ordenar Chr, Start
  bed_PM = bed_PM[order(bed_PM[,1] , bed_PM[,2]),]
  bed_PM[,4] = gsub('\\+', 1, bed_PM[,4])
  bed_PM[,4] = gsub('\\-', -1, bed_PM[,4])
  bed_PM[,4] = as.numeric(bed_PM[,4])
  #Transformando a localizacao de factor para char
  bed_PM[,6] = as.character(bed_PM[,6])
  #substituindo localizacao por numero.
  for (i in 1:nrow(bed_PM)) {
    if(bed_PM[i,6] == "5UTR") {
      bed_PM[i,6] = 1
    } else if (bed_PM[i,6] == "undetermined") {
      bed_PM[i,6] = 2
    } else if (bed_PM[i,6] == "3UTR") {
      bed_PM[i,6] = 3
    } else if (bed_PM[i,6] == "CDS") {
      bed_PM[i,6] = 4
    } else {
      cat(paste0("error >>",bed_PM[i,6],"<< "))
    }
  }
  bed_PM[,6] = as.numeric(bed_PM[,6])
  bed_PM[,7] = as.numeric(bed_PM[,7])
  #Ja feita anteriormente.
  #transformação fc para log2FC
  #bed_PM[,5] = log(bed_PM[,5],2)
  return(bed_PM)
}

bed_PM = bed_func_DEG(tab_chr_DEG_PM, 22)
bed_MA = bed_func_DEG(tab_chr_DEG_MA, 23)
bed_AP = bed_func_DEG(tab_chr_DEG_AP, 24)

#bed expressos em all contrastes
bed_PM_loc = bed_PM[, c(1:4,6,7)]
bed_MA_loc = bed_MA[, c(1:4,6,7)]
bed_AP_loc = bed_AP[, c(1:4,6,7)]

#cbind das colunas de interesse
bed_all_loc = rbind(bed_PM_loc, bed_MA_loc, bed_AP_loc)
bed_all_loc = unique(bed_all_loc[order(bed_all_loc$Chromosome, bed_all_loc$Smaller.coord),])
#separaão de subgrupos por tamanho
id_length_4000 = bed_all_loc$Length <= 4000
table(id_length_4000)
bed_all_loc_4k = bed_all_loc[id_length_4000,]
#decidi usar 2000 pois so sao excluidos 83 ncRNA
id_length_2000 = bed_all_loc$Length <= 2000
table(id_length_2000)
bed_all_loc_2k = bed_all_loc[id_length_2000,]
#calculando media e mediana dos genes removendo scaffods
median(bed_all_loc$Length)
mean(bed_all_loc$Length)
#calculando media e mediana dos genes não removendo scaffods
temp = rbind(tab_DEG_PM, tab_DEG_MA, tab_DEG_AP)
temp_uq = unique(temp)
median(temp_uq$Length)
mean(temp_uq$Length)
table(temp_uq$Length <= 2000)


###RANGE DE Log2FC para DEG
range_DEG = range(bed_PM[,5], bed_MA[,5], bed_AP[,5])

# test_1 = NULL
# for (i in 1:nrow(bed_PM)){
#   if(bed_PM$c7[i] == "+") {
#     test_1[i] = 1
#   } else { test_1[i] = -1 }
# }

####GROUPS!!!
library(openxlsx)
G1_tab = read.xlsx("data/table_ncRNA_G1.xlsx",startRow = 2)
G2_tab = read.xlsx("data/table_ncRNA_G2.xlsx",startRow = 2)
G3_tab = read.xlsx("data/table_ncRNA_G3.xlsx",startRow = 2)
G4_tab = read.xlsx("data/table_ncRNA_G4.xlsx",startRow = 2)
G5_tab = read.xlsx("data/table_ncRNA_G5.xlsx",startRow = 2)
G6_tab = read.xlsx("data/table_ncRNA_G6.xlsx",startRow = 2)
#usando somente as colunas de interesse
bedG1_temp = cbind(G1_tab[,5:8],group = rep(1,nrow(G1_tab)))
bedG2_temp = cbind(G2_tab[,5:8],group = rep(2,nrow(G2_tab)))
bedG3_temp = cbind(G3_tab[,5:8],group = rep(3,nrow(G3_tab)))
bedG4_temp = cbind(G4_tab[,5:8],group = rep(4,nrow(G4_tab)))
bedG5_temp = cbind(G5_tab[,5:8],group = rep(5,nrow(G5_tab)))
bedG6_temp = cbind(G6_tab[,5:8],group = rep(6,nrow(G6_tab)))
#juntando os arquivos
bed_group_temp = rbind(bedG1_temp, bedG2_temp, bedG3_temp, bedG4_temp, bedG5_temp, bedG6_temp)
#removendo Scaffold genes
id_include_bed_g = !grepl("_SCAF", bed_group_temp[,1])
bed_group = bed_group_temp[id_include_bed_g,]
#ordenando e substituindo +/- por 1/-1
bed_group = bed_group[order(bed_group[,1] , bed_group[,2]),]
bed_group[,4] = gsub('\\+', 1, bed_group[,4])
bed_group[,4] = gsub('\\-', -1, bed_group[,4])
bed_group[,4] = as.numeric(bed_group[,4])
bed_group[,2] = as.integer(bed_group[,2])
bed_group[,3] = as.integer(bed_group[,3])


######
#bed CDS 
bed_cds = tab_cds[,1:4]
bed_cds[,4] = gsub('\\+', 1, bed_cds[,4])
bed_cds[,4] = gsub('\\-', -1, bed_cds[,4])
bed_cds[,4] = as.numeric(bed_cds[,4])
#ordenando
bed_cds = bed_cds[order(bed_cds$V1, bed_cds$V2),]
bed_cds = cbind(bed_cds, rep(1:3, nrow(tab_cds)/3)) 
colnames(bed_cds)[5] = "V5"


####Making bed_int - selecionando os genes de interesse em bed_group

ids_search = c("LbrM2903_34_lncRNA380","LbrM2903_33_lncRNA289","LbrM2903_32_lncRNA243","LbrM2903_33_lncRNA177","LbrM2903_30_lncRNA54","LbrM2903_23_lncRNA173","LbrM2903_20.1_lncRNA201","LbrM2903_20.1_lncRNA234")
ids_search = paste(ids_search, collapse = "|")
tab_int_v = tab[grep(ids_search, tab$ncRNA.ID),]

bed_int = tab_int_v[,2:5]
bed_int[,4] = gsub('\\+', 1, bed_int[,4])
bed_int[,4] = gsub('\\-', -1, bed_int[,4])
bed_int[,4] = as.numeric(bed_int[,4])






#Funcão para adicionar cores...
col_fun_CDS = colorRamp2(breaks = c(1, 2, 3), colors = c("grey0", "grey40", "grey80"))
col_fun_DEG = colorRamp2(breaks = c(range_DEG[1], 0, range_DEG[2]), colors = c("dark green", "white", "dark red"))
#1 - G1, 2- G2 ...
col_fun_group = colorRamp2(breaks = c(1,2,3,4,5,6), colors = c("blue2","darkmagenta","green3","red","yellow","cyan"))
#1-5UTR 2-Undetermined 3-3UTR 4-CDS
col_fun_location = colorRamp2(breaks = c(1,2,3,4), colors = c("sienna3","grey","yellow3","coral4"))


################################################################################
################################################################################

#CIRCOS ########################################################################
#test = data.frame(c1 = x_c$V1, c2 = rep(0, 36), c3 = x_c$V2)
#tiff("imageTeste.tiff", height = 500, width = 500, units = "cm", res = 100)


plot_circlize = function(){
pdf("Circos_with_int_19032019.pdf", 1000, 1000)
circos.clear()
circos.par(cell.padding = c(0,0,0,0), track.height = 0.1, 
           track.margin = c(0.001,0.001), 
           canvas.xlim = c(-1.2,1.2), canvas.ylim = c(-1.2,1.2), gap.after = c(rep(1,35), 5))

circos.initialize(factors = x_c3$V1, x = x_c3$V2)
#circos.initialize(factors = bed_cds$V1, x = bed_cds$V2)
#circos.genomicInitialize(bed_cds, plot = F)

#, major.by = NULL, plotType = c("ideogram", "axis", "labels"))
#circos.initializeWithIdeogram(test, plotType = c("axis","labels"))

#text(0,0, "plotType = c('axis','labels')", cex = 1)
#circos.genomicTrack(bed_cds, stack = T,
#                    panel.fun = function(region, value, ...) {
#                      circos.genomicPoints(region, value, pch = 16, cex = 0.3)
#                    })
#col_cds = data.frame(C1 = rep(c("red","blue","green"), nrow(bed_cds)/3))
#col_cds = rep(c("red","blue","green"), 96/3)

#bed_cds = cbind(bed_cds, col_cds)
#colnames(col_cds) = c("V1", "V2")
##TRACK CDS
circos.genomicTrack(bed_cds, ylim = c(-1,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0, col = col_fun_CDS(value[[2]]) , border = col_fun_CDS(value[[2]]))
                      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
                                  CELL_META$sector.index, cex = 60, 
                                  facing = "clockwise", niceFacing = T, adj = c(0,0.5))
                      circos.axis(h = "top", labels.cex = 60, major.tick.length = uy(8,"cm"))
                    })
circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
            "CDS", cex = 110, 
            facing = "clockwise", niceFacing = T)

##TRACK DEG 
#Função TRACK DEG
track_circos_DEG = function(bed_PM, col_fun_DEG,text_v) {
      circos.genomicTrack(bed_PM, ylim = c(-1,1), bg.border = "black", bg.col = "grey95",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0, 
                                         col = col_fun_DEG(value[[2]]) , border = col_fun_DEG(value[[2]]))
                      
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      #                      panel.fun = function(x,y) {
                      #circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
                      #            CELL_META$sector.index, cex = 80, las = 3)
                      #circos.axis(labels.cex = 60)
                    })
  circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
              text_v , cex = 90, 
              facing = "clockwise", niceFacing = T)
}
track_circos_DEG(bed_PM, col_fun_DEG, "P/M")
track_circos_DEG(bed_MA, col_fun_DEG, "M/A") 
track_circos_DEG(bed_AP, col_fun_DEG, "A/P")


############int


circos.par(cell.padding = c(0,0,0,0), track.height = 0.04, 
           track.margin = c(0,0.002))

circos.genomicTrack(bed_int, ylim = c(-1,1), bg.border = NA,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0, 
                                         col = "black" , border = "black")
 #                     circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      #                      panel.fun = function(x,y) {
                      #circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
                      #            CELL_META$sector.index, cex = 80, las = 3)
                      #circos.axis(labels.cex = 60)
                    })


for (i in c(1:19,20.1,20.2,21:35)){
  track_idx = 5
  i = formatC(i, width = 2, flag = "0")
  xlim_v = get.cell.meta.data("xlim", sector.index = paste0("LbrM2903_",i),track.index = track_idx)
  ylim_v = get.cell.meta.data("ylim", sector.index = paste0("LbrM2903_",i),track.index = track_idx)
  circos.lines(c(0,0), c(ylim_v[1]-uy(4,"cm"),ylim_v[2]+uy(4,"cm")) , sector.index = paste0("LbrM2903_",i),track.index = track_idx, col = "black")
  circos.lines(c(xlim_v[2],xlim_v[2]), c(ylim_v[1]-uy(4,"cm"),ylim_v[2]+uy(4,"cm")), sector.index = paste0("LbrM2903_",i),track.index = track_idx , col = "black")
  # circos.lines(xlim_v, c(1,1) , sector.index = paste0("LbrM2903_",i),track.index = track_idx, col = "black")
  circos.lines(xlim_v, c(0, 0), sector.index = paste0("LbrM2903_",i),track.index = track_idx, lty = 2, col = "#00000040")
  }


circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
            "**" , cex = 90, 
            facing = "clockwise", niceFacing = T)




circos.par(cell.padding = c(0,0,0,0), track.height = 0.1, 
           track.margin = c(0.001,0.001))



#####GROUPS
circos.genomicTrack(bed_group, ylim = c(-1,1), bg.col = "grey80",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0,
                                         col = col_fun_group(value[[2]]) , border = col_fun_group(value[[2]]))
#                      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
#                                  CELL_META$sector.index, cex = 60, 
#                                  facing = "clockwise", niceFacing = T, adj = c(0,0.5))
#                      circos.axis(h = "top", labels.cex = 60, major.tick.length = uy(60))
                    })
circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
            "3xDE", cex = 80, 
            facing = "clockwise", niceFacing = T)

####Circos Point
circos.genomicTrack(bed_all_loc_2k , ylim = range(bed_all_loc_2k$Length),
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, numeric.column = 3, 
                                           pch = 16, cex = 10, col = (col_fun_location(value[[2]])))
                      circos.lines(CELL_META$cell.xlim, c(200, 200), lty = 2, col = "black")
                      ####PAREI AQUI FAZER O GRAFICO DE PONTOS>
                      ####
                                            
           #           circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
          #                        CELL_META$sector.index, cex = 60, 
          #                        facing = "clockwise", niceFacing = T, adj = c(0,0.5))
           #           circos.axis(h = "top", labels.cex = 60, major.tick.length = uy(60))
                    })
circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
            "L", cex = 70, 
            facing = "clockwise", niceFacing = T)


####LEGEND
#circos.initialize(factors = x_c3$V1, x = x_c3$V2)
#circos.genomicTrack(bed_group, ylim = c(-1,1), bg.col = "grey80",
#                    panel.fun = function(region, value, ...) {
#                      circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0,
#                                         col = col_fun_group(value[[2]]) , border = col_fun_group(value[[2]]))
                      #                      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
                      #                                  CELL_META$sector.index, cex = 60, 
                      #                                  facing = "clockwise", niceFacing = T, adj = c(0,0.5))
                      #                      circos.axis(h = "top", labels.cex = 60, major.tick.length = uy(60))
#                    })

dev.off()
}
plot_print = plot_circlize()



#LEGENDA
library(ComplexHeatmap)
library(gridBase)

legend_func = function(){
tiff(filename = "Legend_19032019.tiff")
lgd_CDS = Legend(at = c(" "," "," "), type = "lines", #labels_gp = gpar(fontsize= 1000),
                    legend_gp = gpar(col = c("grey0", "grey40", "grey80"), cex = 60, cex.pt = 20) , title_position = "topleft", 
                    title = "Track 1\nCDS", nrow = 1)#, grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))

 lgd_DEG = Legend(at = c(round(range_DEG[1], digits = 2), 0, round(range_DEG[2], digits = 2)), col_fun = col_fun_DEG, #legend_gp = gpar(cex = 20), #labels_gp = gpar(fontsize= 1000),
                  title_position = "topleft", title = "Track 2,3,4\nLog2FC")#, grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))
 
 lgd_ncRNA =  Legend(at = c(" "), type = "lines", #labels_gp = gpar(fontsize= 1000),
                     legend_gp = gpar(col = c("black"), cex = 60, cex.pt = 20) , title_position = "topleft", 
                     title = "Track 5\nncRNA of Interest(**)", nrow = 1)#, grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))
 
 
 lgd_group = Legend(at = c("G1","G2","G3","G4","G5","G6"), type = "lines", #labels_gp = gpar(fontsize= 1000),
                    legend_gp = gpar(col = c("blue2","darkmagenta","green3","red","yellow","cyan"), cex = 20) , title_position = "topleft", 
                    title = "Track 6\n3xDE - Groups") #, grid_height = unit(0.6, "npc"))#, grid_width = unit (0.1, "npc"))
 
 lgd_length = Legend(at = c("5`UTR","Undetermined","3`UTR","CDS"), type = "points", #labels_gp = gpar(fontsize= 1000),
                     legend_gp = gpar(col = c("sienna3","grey","yellow3","coral4"), cex = 20) , title_position = "topleft", 
                     title = "Track 7\nL (Location)") #,grid_height = unit(0.6, "npc"))#, grid_width = unit (0.1, "npc"))



 lgd_all = packLegend(lgd_CDS, lgd_DEG,lgd_ncRNA, lgd_group, lgd_length)

# pushViewport(viewport(x = unit(0.8, "npc"), y = unit(0.8, "npc"), 
#                       width = grobWidth(lgd_CDS), 
#                       height = grobHeight(lgd_CDS), 
#                       just = c("right", "bottom")))
par(cex = 10)
grid.draw(lgd_all)
#grid.draw(lgd_CDS)
#upViewport()
dev.off()
}
legend_print = legend_func()

############################################################################
#############################################################################







##################################
##################################




# pdf("imageTeste3.pdf", 1000, 1000)
# 
# plot.new()
# circle_size = unit(1, "snpc") # snpc unit gives you a square region
# pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
#                       just = c("center", "top")))
# par(omi = gridOMI(), new = TRUE)
# 
# circos.clear()
# circos.par(cell.padding = c(0,0,0,0), track.height = 0.1, 
#            track.margin = c(0.001,0.001), 
#            canvas.xlim = c(-1.2,1.2), canvas.ylim = c(-1.2,1.2), gap.after = c(rep(1,35), 5))
# 
# circos.initialize(factors = x_c3$V1, x = x_c3$V2)
# 
# #####GROUPS
# circos.genomicTrack(bed_group, ylim = c(-1,1), bg.col = "grey80",
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicRect(region, value, ytop.column = 1 , ybottom = 0,
#                                          col = col_fun_group(value[[2]]) , border = col_fun_group(value[[2]]))
#                       #                      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(50,"cm"),
#                       #                                  CELL_META$sector.index, cex = 60, 
#                       #                                  facing = "clockwise", niceFacing = T, adj = c(0,0.5))
#                       #                      circos.axis(h = "top", labels.cex = 60, major.tick.length = uy(60))
#                     })
# circos.text(CELL_META$xlim[2] + 350000 , CELL_META$ycenter,
#             "3xDE", cex = 80, 
#             facing = "clockwise", niceFacing = T)
# 
# upViewport()
# 
# 
# lgd_CDS = Legend(at = c(" "," "," "), type = "lines", #labels_gp = gpar(fontsize= 1000),
#                  legend_gp = gpar(col = c("grey0", "grey40", "grey80"), cex = 60) , title_position = "topleft", 
#                  title = "Track 1\nCDS", nrow = 1, grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"), size = unit(0.05, "npc"))
# 
# lgd_DEG = Legend(at = c(round(range_DEG[1], digits = 2), 0, round(range_DEG[2], digits = 2)), col_fun = col_fun_DEG,legend_gp = gpar(cex = 60), #labels_gp = gpar(fontsize= 1000),
#                  title_position = "topleft", title = "Track 2,3,4\nLogFC", grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))
# 
# lgd_group = Legend(at = c("G1","G2","G3","G4","G5","G6"), type = "lines", #labels_gp = gpar(fontsize= 1000),
#                    legend_gp = gpar(col = c("blue2","darkmagenta","green3","red","yellow","cyan"), cex = 60) , title_position = "topleft", 
#                    title = "Track 4\nGroups", grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))
# 
# lgd_length = Legend(at = c("5`UTR","Undetermined","3`UTR","CDS"), type = "points", #labels_gp = gpar(fontsize= 1000),
#                     legend_gp = gpar(col = c("sienna3","grey","yellow3","coral4"), cex = 60) , title_position = "topleft", 
#                     title = "Track 5\nLength", grid_height = unit(0.6, "npc"), grid_width = unit (0.1, "npc"))
# 
# 
# 
# lgd_all = packLegend(lgd_CDS, lgd_DEG, lgd_group, lgd_length)
# 
# 
# pushViewport(viewport(x = 0.5, y = unit(1, "npc") - circle_size, 
#                       width = grobWidth(lgd_all), height = grobHeight(lgd_all), 
#                       just = c("center", "top")))
# grid.draw(lgd_all)
# upViewport()
# 
# 
# 
# 
# 
# pushViewport(viewport(x = unit(0.8, "npc"), y = unit(0.2, "npc"), 
#                       width = grobWidth(lgd_list_vertical), 
#                       height = grobHeight(lgd_list_vertical), 
#                       just = c("right", "bottom")))
# grid.draw(lgd_all)
# upViewport()
# 
# dev.off()
# 
# 
# # #PM
# # index_PM = match(tab_DEG_PM$V1, tab_cood_all$V5)
# # tab_cood_DEG = cbind(tab_DEG, tab_cood_all[index_PM,])
# # #removendo os scafolds
# # id_exclude = !grepl("_SCAF", tab_cood_DEG$V1)
# # tab_chr_DEG = tab_cood_DEG[id_exclude,]
# # #removendo os scaffolds do factor.
# # head(tab_chr_DEG$V1) 
# # levels_v1 = x_c3$V1
# # fact = factor(as.character(tab_chr_DEG$V1))
# # levels(fact)
# # tab_chr_DEG$V1 = fact 
# 
# #iniciando o circos com o tamanho de cada chrm
# #circos.clear()
# #circos.par( "track.height" = 0.05)
# # circos.initialize(factors = x_c3$V1, x = x_c3$V2)
# # 
# # #circos.initialize(factors = tab_chr_DEG$V1, tab_chr_DEG$V2)
# # # circos.track(factors = tab_chr_DEG$V1, x = tab_chr_DEG$V2, y = tab_chr_DEG$FoldChange,  
# # #              panel.fun = function(x,y) {
# # #                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(3,"mm"), 
# # #                            CELL_META$sector.index, cex = 0.4, las = 3)
# # #                circos.axis(labels.cex = 0.2)
# # #              })
# # circos.par( "track.height" = 0.1)
# # 
# # ##Plot CDS
# # circos.track(factors = tab_cds$V1, x = tab,ylim = c(0,2),  
# #              panel.fun = function(x,y) {
# #                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(3,"mm"), 
# #                            CELL_META$sector.index, cex = 0.4, las = 3)
# #                circos.axis(labels.cex = 0.2)
# # 
# #                circos.rect()
# #              })  
# # circos.trackPoints(tab_cds$V1, tab_cds$V2 ,rep(1, length(tab_cds$V1)) , pch = 16, cex = 0.5)
# # circos.genomicRect(tab_cds$V2, 0, tab_cds$V3, 1, tab_cds$V1)
# # 
# # circos.rect(tab_cds$V2, rep(1, length(tab_cds$V1)), tab_cds$V3, rep(2, length(tab_cds$V1)), sector.index = tab_cds$V1)
# # 
# # 
# # ###Plot DEG
# # #rangeFC global
# # rangeUP_FC = ceiling(max(tab_chr_DEG_PM$c3, tab_chr_DEG_MA$c3, tab_chr_DEG_AP$c3))
# # 
# # ###Plot P_M DEG
# # circos.track(factors = x_c3$V1, ylim = c(0,5))
# # # ,  
# # #              panel.fun = function(x,y) {
# # #                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(3,"mm"), 
# # #                            CELL_META$sector.index, cex = 0.4, las = 3)
# # #                circos.axis(labels.cex = 0.2)
# # #              })  
# # 
# # ###########consigo simular o tamanho dos chrm e controlar o tamanho das tracks.
# # #col = rep(c("#FF0000", "#00FF00"), 32)
# # #por enquanto sem cor.
# # circos.trackPoints(tab_chr_DEG_PM[,4], tab_chr_DEG_PM[,5], tab_chr_DEG_PM[,3] , pch = 16, cex = 0.5)
# # 
# # ###Plot M_A DEG
# # circos.track(factors = x_c3$V1, ylim = c(0,5))  
# # circos.trackPoints(tab_chr_DEG_MA[,4], tab_chr_DEG_MA[,5], tab_chr_DEG_MA[,3] , pch = 16, cex = 0.5)
# # 
# # ###Plot A_P DEG
# # circos.track(factors = x_c3$V1, ylim = c(0,5))  
# # circos.trackPoints(tab_chr_DEG_AP[,4], tab_chr_DEG_AP[,5], tab_chr_DEG_AP[,3] , pch = 16, cex = 0.5)
# # 
# # 
# # 
# # 
# # ###########Tutorial
# # 
# # circos.track(factors = x_c3$V1, ylim = c(0,20), track.index =5, track.bg = "red")
# # 
# # 
# # circos.track(factors, ylim = c(0, 1), track.index = 1, ...)
# # 
# # 
# # #circos.trackPoints(tab_chr_DEG$V1, tab_chr_DEG$V2, tab_chr_DEG$FoldChange , pch = 16, cex = 0.5)
# # 
# # # circos.trackLines(df$factors, df$x, df$y, col = col, pch = 16, cex = 0.5)
# # # 
# # # range (rep(0, 36), x_c$V2)
# # # n = c(2, 3, 5) 
# # #  s = c("aa", "bb", "cc", "dd", "ee") 
# # #  b = c(TRUE, FALSE, TRUE, FALSE, FALSE) 
# # #  x = list(n, s, b, 3)   # x contains copies of n, s, b 
# # #  
# # #  x_c$V1[1]
# # # i = 2
# # #   
# # #  
# # # a = list()
# # #  for (i in 1:36) {
# # #    
# # #    a[i] = list(c(0,x_c$V2[i]))
# # #  }
