#Rscript <program> <S1_ChIP> <S2_ChIP> <S1_SHU> <S2_SHU>

args = commandArgs(trailingOnly = TRUE)
print(args)


Tab_ChIP_Ty = read.table(file=args[1], header=T, sep="\t")
Tab_ChIP_TyM = read.table(file=args[2], header=T, sep="\t")

Tab_Shu_Ty = read.table(file=args[3], header=T, sep="\t")
Tab_Shu_TyM = read.table(file=args[4], header=T, sep="\t")

#boxplot S1xS1SHU
pdf(paste0("S1xS1SHU_Histogram.pdf") , width=8, height = 6)
boxplot(Tab_ChIP_Ty$distancetoFeature, Tab_Shu_Ty$distancetoFeature, names=c("Ty_ChIP","Ty_Shuffle"), main= "Box-plot Distance Feature", ylim=c(-50000,50000))
dev.off()

#boxplot S2xS2SHU
pdf(paste0("S2xS2SHU_Histogram.pdf") , width=8, height = 6)
boxplot(Tab_ChIP_TyM$distancetoFeature, Tab_Shu_TyM$distancetoFeature, names=c("TyM_ChIP","TyM_Shuffle"), main= "Box-plot Distance Feature", ylim=c(-50000,50000))
dev.off()


kstesteS1 = ks.test(Tab_ChIP_Ty$distancetoFeature, Tab_Shu_Ty$distancetoFeature, alternative = c("two.sided"), exact = NULL)
kstesteS2 = ks.test(Tab_ChIP_TyM$distancetoFeature, Tab_Shu_TyM$distancetoFeature, alternative = c("two.sided"), exact = NULL)

wiltesteS1 = wilcox.test(Tab_ChIP_Ty$distancetoFeature, Tab_Shu_Ty$distancetoFeature)
wiltesteS2 = wilcox.test(Tab_ChIP_TyM$distancetoFeature, Tab_Shu_TyM$distancetoFeature)


