#Rscript Rename_files_2021_04_22.R xlsx_PlateMap dir_seq_ab1 dir_seq_out Sheet_number(optional) complement_to_outname(ex sample"_2", optional)
#setwd("/media/sf_F_DRIVE/Dropbox/CTvacinas/Results/Plates_Remote/")
args = commandArgs(trailingOnly = T)

#args[1] = "Placa Renata.xlsx"
#args[2] = "CTV_2021.08.26_2021-08-26HELENA"
#args[3] = "CTV_2021.08.26_2021-08-26HELENA_renamed"
#args[4] = 1

if (length(args) < 3) {
  stop("At least three arguments must be supplied (input: xlsx_PlateMap dir_seq_ab1 dir_seq_out)", call.=FALSE)
}
if(length(args) < 4){
  args[4] = 1
  args[5] = ""
}
if(length(args) < 5){
  args[5] = ""
}

library(openxlsx)
library(glue)
pl_design = read.xlsx(args[1], sheet = as.numeric(args[4]))

if(!exists("pl_design")){
  stop("This plate table is not in correct format.")
}
cat("\nArguments:\n")
print(args)
#View(pl_design)

ct_v = 1
plat_v = NULL
nam_v = NULL
for (j in 2:13){
  for (i in 2:9){
    if(is.null(pl_design[i,j])){
      stop("This plate table is not in correct format or this is not a correct table sheet(please specify the number of correct sheet in the last argument).")
    }
    #if(grepl("^[0-9][0-9][0-9][0-9]$", pl_design[i,j])){
    if(!is.na(pl_design[i,j])){
      plate_c = glue(toupper(letters[i-1]),sprintf("%02d",j-1))
      nam_v = c(nam_v, plate_c)
      plat_v = c(plat_v, paste0(pl_design[i,j],colnames(pl_design)[j],"_",plate_c,args[5]))
      ct_v = ct_v + 1
      #print(plate_id)
      #print(name_v)
    } 
  }
}
nam_v
plat_v
names(plat_v) = nam_v

folder_v = args[2]
folder_out_v = args[3]

list_f = list.files(path = folder_v, pattern = "*.ab1")
#list_f


idx_list = NULL
for (i in 1:length(nam_v)){
#  test_v = "x"
  test_v = grep(nam_v[i],list_f)
#  ct = i
#  print(test_v)
  if (length(test_v) == 0) {
      idx_list = c(idx_list,NA) 
   }
   else { idx_list = c(idx_list,grep(nam_v[i],list_f))}
}

#idx_list

file_match = list_f[idx_list]

#file_match[90]
#nam_v[90]
#plat_v[90]

dir.create(folder_out_v)
for(j in 1:length(plat_v)){
  file.copy(from = file.path(folder_v,file_match[j]),to = file.path(folder_out_v,paste0(plat_v[j],".ab1")))
}
