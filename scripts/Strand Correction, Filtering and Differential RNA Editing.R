### Master script for RNA editing meta-analysis

# .csv files produced by table_annovar.pl script

# converted files that include:
# .vcf files converted with convert2annovar.pl (do not forget --includeinfo flag!)
# outSigTable files converted with R for ANNOVAR

## Set variables:

# Define "human" or "mouse" samples for required files
# input_RG is assigned with the .csv file, extracted from ENSEMBLE .gtf file, with strand info
# RepeatMask for RepeatMasket.txt table with repeat sequences info (Alu elements, etc)

genome = "mouse"

library(data.table)

if (genome == "mouse"){
  
  input_RG = read.csv("./required_files/mouse/Mus_musculus_ensembl_dataset.csv")
  input_RG = input_RG[!duplicated(input_RG$gene_name),]
  
  RepeatMask = read.table("./required_files/mouse/mm10_RepeatMasker.txt")
  
  cDNA_snp = fread("./required_files/mouse/mm10_cDNA_snp142.txt")
  
  radar = fread("./required_files/mouse/Converted_with_R_mm10_Mouse_AG_all_v2.txt")
  
  C_to_U = fread("./required_files/mouse/Non-redundant_known_C-to-U_RNA_editing_positions")

} else if (genome == "human"){
  
  input_RG = read.csv("./required_files/human/Homo_sapiens.GRCh38.83.csv")
  input_RG = input_RG[!duplicated(input_RG$gene_name),]
  
  cDNA_snp = fread("./required_files/mouse/hg38_cDNA_snp144.txt")
  
  RepeatMask = read.table("./required_files/human/hg38_RepeatMasker.txt")
  
  radar = fread("./required_files/human/Converted_with_R_hg38_Human_AG_all_v2.txt")
}

# Define number of parallel jobs (for the step: # add RepeatMasker file info)

jobs = 3

# Define threshold variable

threshold_var = 2

# Define minimum threshold for variation frequency

var_freq = 0

# Dual phenotypes for concatenated tables: t.test

dual_phenotypes = list(c("control","cjd"))

# Dual stages for concatenated tables: t.test

dual_stages = list(c("120d","180d"))

# Define tools used for variant detection

tools = c("REDItools", "VarScan")

# Define sample phenotypes (control, AD, CJD, etc)

phenotypes = c("control", "cjd")

# Define sample stages (120d, 180d, etc)

stages = c("120d","180d")

# Define sample names (B1, B2, etc)
# Sample name MUST exist within the original file

sample_names = c("B1","B2","B3","B4",
                 "C1", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4",
                 "E1", "E2", "E3", "E4", "E5")

#

## Scripts:

# Start: rename files

dir.create("./renamed_files/csv_files/unprocessed", recursive = TRUE)
dir.create("./renamed_files/csv_files/with_strand", recursive = TRUE)
dir.create("./renamed_files/csv_files/with_RepeatMask", recursive = TRUE)
dir.create("./renamed_files/converted_files/unprocessed", recursive = TRUE)

for(i1 in tools){
  for(i2 in phenotypes){
    for(i3 in stages){
      for(i4 in sample_names){
        
        # current directory defined by i1, i2, i3 and i4
        
        # list of current .csv files to be moved and renamed
        d1_csv = list.files(paste(".","original_files/csv_files",i1,i2,i3, sep = "/"))
        # list of current .csv files to be moved and renamed (full path names!)
        d2_csv = list.files(full.names = TRUE ,paste(".","original_files/csv_files",i1,i2,i3, sep = "/"))
        # list of current converted files to be moved and renamed
        d1_convert = list.files(paste(".","original_files/converted_files",i1,i2,i3, sep = "/"))
        # list of current converted files to be moved and renamed (full path names!)
        d2_convert = list.files(full.names = TRUE ,paste(".","original_files/converted_files",i1,i2,i3, sep = "/"))
        
        # pick with grep function the file containing the sample_name
        # if loops for picking only the samples that exist in the current directory
        
        j1 = grep(i4, d1_csv)
        
        if(length(j1) == 0){print(paste0("sample ",paste(i1,i2,i3,i4,sep=".")," does not exist"))
        } else if (length(j1) == 1){
          file.copy(from = d2_csv[j1], to = paste(sep = "/",".","renamed_files","csv_files","unprocessed",paste(sep=".",i1,i2,i3,i4,"csv")))
        }
        
        j2 = grep(i4, d1_convert) 
        
        if(length(j2) == 0){print(paste0("sample ",paste(i1,i2,i3,i4,sep=".")," does not exist"))
        } else if (length(j2) == 1){
          file.copy(from = d2_convert[j2], to = paste(sep = "/",".","renamed_files","converted_files","unprocessed",paste(sep=".",i1,i2,i3,i4)))
        }
      }
    }
  }
}

# End: rename files

# Start: add strand info

source(file = "./R_Scripts/add_strand_info.R")

l_unprocessed = length(list.files("./renamed_files/csv_files/unprocessed"))

for(i in seq(l_unprocessed)){
  add_strand_info(TA = list.files(full.names = TRUE,"./renamed_files/csv_files/unprocessed")[i],
                  RG = input_RG,
                  output = paste0("./renamed_files/csv_files/with_strand/",list.files("./renamed_files/csv_files/unprocessed")[i]))
}

# End: add strand info

# Start: add RepeatMasker file info (Alu elements, etc)

source(file = "./R_Scripts/add_RepeatMask_info.R")

library(parallel)

# create a list with empty vectors, which number is equal to jobs
tasks = as.list(seq(jobs))

# find the number of files in /with_strand/ directory
l_with_strand = length(list.files("./renamed_files/csv_files/with_strand"))

# create limits of numerical vectors based on the value of the variable jobs
ind_with_strand = round(seq(1,l_with_strand, len = jobs+1))

# create subsetting numerical vectors (1:x, etc) for diving the files
# assign these vectors to the empty vectors of tasks list
for(j in seq(jobs)){
  if(j == 1){
    tasks[[j]] = 1:ind_with_strand[j+1]
  } else {
    tasks[[j]] = (ind_with_strand[j]+1):ind_with_strand[j+1]
  }
}

# Use the mclapply function of parallel package to run simultaneously
# the custom function add_RepeatMask for the divided files
mclapply(
  tasks,
  function(x){
    for(i in x){
      add_RepeatMask(TA = list.files(full.names = TRUE,"./renamed_files/csv_files/with_strand")[i],
                     RM = RepeatMask,
                     output = paste0("./renamed_files/csv_files/with_RepeatMask/",list.files("./renamed_files/csv_files/with_strand")[i]))
    }},
  mc.cores = length(tasks)
)

# remove RepeatMasker table after extracting the information
rm(RepeatMask)

# End: add RepeatMasker file info (Alu elements, etc)

# Start: Filter for RNA editing positions

# Create directories with filtered RNA editing positions
for(i in c("ADAR_all_positions/xlsx",
           "ADAR_all_positions/tables",
           "APOBEC_all_positions/xlsx",
           "APOBEC_all_positions/tables")){
  dir.create(recursive = TRUE, paste0("./renamed_files/",i))
}

source("./R_Scripts/RNA_editing_filter.R")

# store the lists with the .csv file and the converted file names
# create two lists, one with full path and one without
a1 = list.files("./renamed_files/csv_files/with_RepeatMask")
a2 = list.files(full.names = TRUE,"./renamed_files/csv_files/with_RepeatMask")
b1 = list.files("./renamed_files/converted_files/unprocessed")
b2 = list.files(full.names = TRUE,"./renamed_files/converted_files/unprocessed")

# call custom function RNA_editing_filter to extract all possible A-to-I RNA editing positions
for(i in seq(length(a1))){
  RNA_editing_filter(table_annovar = a2[i],
                     index = b2[i],
                     pipeline = (
                       if(length(grep("REDItools",a1[i]))==1){
                         "REDItools"
                       } else if(length(grep("VarScan",a1[i]))==1){
                         "VarScan"
                       }
                     ),
                     freq = var_freq,
                     enzyme = "ADAR",
                     output = b1[i])
}

# call custom function RNA_editing_filter to extract all possible C-to-U RNA editing positions
for(i in seq(length(a1))){
  RNA_editing_filter(table_annovar = a2[i],
                     index = b2[i],
                     pipeline = (
                       if(length(grep("REDItools",a1[i]))==1){
                         "REDItools"
                       } else if(length(grep("VarScan",a1[i]))==1){
                         "VarScan"
                       }
                     ),
                     freq = var_freq,
                     enzyme = "APOBEC",
                     output = b1[i])
}

# remove these variables after extracting their respective information
rm(cDNA_snp)
rm(radar)
rm(C_to_U)
rm(a1)
rm(a2)
rm(b1)
rm(b2)

# End: Filter for RNA editing positions

# Start: Merging samples with same enzyme, tool, phenotype, stage

dir.create(recursive = TRUE,"./merged_samples/ADAR/xlsx")
dir.create(recursive = TRUE,"./merged_samples/ADAR/tables")
dir.create(recursive = TRUE,"./merged_samples/APOBEC/xlsx")
dir.create(recursive = TRUE,"./merged_samples/APOBEC/tables")

# Create lists with the files containing the RNA editing positions for each sample (full path included)
a1 = list.files(full.names = TRUE, "./renamed_files/ADAR_all_positions/tables")
b1 = list.files(full.names = TRUE, "./renamed_files/APOBEC_all_positions/tables")

source("./R_Scripts/Merging_tables.R")

for(i1 in tools){
  for(i2 in phenotypes){
    for(i3 in stages){
      
      # ADAR samples
      l1 = a1 %>%
        (function(x){x[grep(i1,x)]})%>%
        (function(x){x[grep(i2,x)]})%>%
        (function(x){x[grep(i3,x)]})
      # APOBEC samples
      l2 = b1 %>%
        (function(x){x[grep(i1,x)]})%>%
        (function(x){x[grep(i2,x)]})%>%
        (function(x){x[grep(i3,x)]})
      
      # Very careful in this step: the names of the samples should exist
      # only within the names of the files and only once!
      
      # Sample names for ADAR
      existing_samples.1 = c()
      for(j in seq(length(sample_names))){
        if(length(grep(sample_names[j],l1)) == 1){
          existing_samples.1 = c(existing_samples.1,
                                 sample_names[j])
        }
      }
      # Sample names for APOBEC
      existing_samples.2 = c()
      for(j in seq(length(sample_names))){
        if(length(grep(sample_names[j],l2)) == 1){
          existing_samples.2 = c(existing_samples.2,
                                 sample_names[j])
        }
      }
      #Merge ADAR tables
      merge_tables(
        no.tables = length(l1),
        table_names = l1,
        index = paste0(".",existing_samples.1),
        mean_name = paste("Mean_Var_freq",i1,i2,i3,sep="."),
        output = paste(i1,i2,i3,sep=".")
      )
      #Merge APOBEC tables
      merge_tables(
        no.tables = length(l2),
        table_names = l2,
        index = paste0(".",existing_samples.2),
        mean_name = paste("Mean_Var_freq",i1,i2,i3,sep="."),
        output = paste(i1,i2,i3,sep=".")
      )
    }
  }
}

rm(a1)
rm(b1)

# End: Merging samples with same enzyme, tool, phenotype, stage

# Start: Merging REDItools and VarScan tables

for(i1 in phenotypes){
  for(i2 in stages){
    for(i3 in c("ADAR","APOBEC")){
      ind_tools.tables = grep(paste0(".",i1,".",i2),
                              list.files(paste0("./merged_samples/",i3,"/tables")))
      tools.tables = list.files(full.names = TRUE,
                                paste0("./merged_samples/",i3,"/tables"))[ind_tools.tables]
      
      REDItool_table = read.table(tools.tables[grep(".REDItools",tools.tables)],
                                  header = TRUE)
      VarScan_table = read.table(tools.tables[grep(".VarScan",tools.tables)],
                                 header = TRUE)
      if(i3 == "ADAR"){
        col_end.1 = which(colnames(REDItool_table) == "RADAR")
        col_end.2 = which(colnames(VarScan_table) == "RADAR")
        
        merged_tool_table = base::merge(REDItool_table,VarScan_table, 
                                        by.x = colnames(REDItool_table)[1:col_end.1],
                                        by.y = colnames(VarScan_table)[1:col_end.2],
                                        suffixes = paste0(".",tools))
        
        
        write.table(merged_tool_table,
                    paste0("./merged_samples/",i3,"/tables/Merged_",i3,
                           ".REDItools.VarScan.",i1,".",i2),
                    row.names = FALSE)
        write.xlsx(merged_tool_table,
                   paste0("./merged_samples/",i3,"/xlsx/Merged_",i3,
                          ".REDItools.VarScan.",i1,".",i2,".xlsx"),
                   row.names = FALSE)
        
        l_reditools = sapply(seq(nrow(REDItool_table)),
                             function(y){paste0(REDItool_table$Chr[y], "-", REDItool_table$Start[y])})
        
        l_varscan = sapply(seq(nrow(VarScan_table)),
                           function(y){paste0(VarScan_table$Chr[y], "-", VarScan_table$Start[y])})
        
        library(VennDiagram)
        
        overlap = calculate.overlap(
          x = list(l_reditools,
                   l_varscan)
        )
        
        venn.plot=draw.pairwise.venn(length(overlap[[1]]),
                                     length(overlap[[2]]), 
                                     length(overlap[[3]]),
                                     lty = "blank",
                                     fill = c("blue","red"),
                                     category = c("REDItools","VarScan"),
                                     cex=1.5);
        dev.copy(png,paste0("./merged_samples/",i3,"/xlsx/","Venn_plot_",i3,
                            ".REDItools.VarScan.",i1,".",i2,".png"));dev.off();dev.off()
      } 
      if(i3 == "APOBEC") {
        if(genome == "mouse"){
          col_end.1 = which(colnames(REDItool_table) == "known_C_to_U")
          col_end.2 = which(colnames(VarScan_table) == "known_C_to_U")
        } else if (genome == "human"){
          col_end.1 = which(colnames(REDItool_table) == "repFamily.2")
          col_end.2 = which(colnames(VarScan_table) == "repFamily.2")
        }
        
        merged_tool_table = base::merge(REDItool_table,VarScan_table, 
                                        by.x = colnames(REDItool_table)[1:col_end.1],
                                        by.y = colnames(VarScan_table)[1:col_end.2],
                                        suffixes = paste0(".",tools))
        
        
        write.table(merged_tool_table,
                    paste0("./merged_samples/",i3,"/tables/Merged_",i3,
                           ".REDItools.VarScan.",i1,".",i2),
                    row.names = FALSE)
        write.xlsx(merged_tool_table,
                   paste0("./merged_samples/",i3,"/xlsx/Merged_",i3,
                          ".REDItools.VarScan.",i1,".",i2,".xlsx"),
                   row.names = FALSE)
        
        l_reditools = sapply(seq(nrow(REDItool_table)),
                             function(y){paste0(REDItool_table$Chr[y], "-", REDItool_table$Start[y])})
        
        l_varscan = sapply(seq(nrow(VarScan_table)),
                           function(y){paste0(VarScan_table$Chr[y], "-", VarScan_table$Start[y])})
        
        library(VennDiagram)
        
        overlap = calculate.overlap(
          x = list(l_reditools,
                   l_varscan)
        )
        
        venn.plot=draw.pairwise.venn(length(overlap[[1]]),
                                     length(overlap[[2]]), 
                                     length(overlap[[3]]),
                                     lty = "blank",
                                     fill = c("blue","red"),
                                     category = c("REDItools","VarScan"),
                                     cex=1.5);
        dev.copy(png,paste0("./merged_samples/",i3,"/xlsx/","Venn_plot_",i3,
                            ".REDItools.VarScan.",i1,".",i2,".png"));dev.off();dev.off()
      }
    }
  }
}

# End: Merging REDItools and VarScan tables

# Start: Concatenate all samples with the same tool

source("./R_Scripts/concatenate_tables.R")

for (i1 in c("ADAR", "APOBEC")){
  
  ind1 = grep(i1, list.files("./renamed_files"))
  
  path1 = paste0(list.files(full.names = TRUE,"./renamed_files")[ind1],
                 "/tables")
  
  for(i2 in tools){
    
    path2 = paste0("./concatenated_samples/",i1,"/",i2, "/tables")
    path3 = paste0("./concatenated_samples/",i1,"/",i2, "/xlsx")
    dir.create(recursive = TRUE, path2)
    dir.create(recursive = TRUE, path3)
    
    ind_tables=c()
    for(i3 in phenotypes){
      for(i4 in stages){
        
        ind2 = grep(paste(i2,i3,i4, sep = "."),
                    list.files(path1))
        
        ind_tables = c(ind_tables, ind2)
        
      }
    }
    
    table_names = list.files(full.names = TRUE,
                             path1)[ind_tables]
    
    index = c()
    for(j in 1:length(strsplit(table_names,"[.]"))){
      a=strsplit(table_names,"[.]")[[j]]
      index = c(index, paste0(".",a[length(a)-2],".",a[length(a)-1],".",a[length(a)]))
    }
    
    output = "."
    for(j1 in 1:length(phenotypes)){
      if(j1 == 1){
        output = paste0(output, phenotypes[j1])
      } else {
        output = paste(output, phenotypes[j1], sep = ".")
      }
    }
    for(j2 in 1:length(stages)){
      if(j1 == 1){
        output = paste0(output, stages[j2])
      } else {
        output = paste(output, stages[j2], sep = ".")
      }
    }
    output = paste0(i1,output)
    
    
    if (i1 == "ADAR"){
      ind_end = "RADAR"
    } else if(i1 == "APOBEC"){
      if(genome == "mouse"){
        ind_end = "known_C_to_U"
      } else {
        ind_end = c()
      }
    }
    #
    concatenate_tables(no.tables = length(table_names),
                       table_names = table_names,
                       index = index,
                       output = output)
    
  }
}

# End: Concatenate all samples with the same tool

# Start: Concatenate all samples from both tools

source("./R_Scripts/concatenate_tables.R")

library(dplyr)

for (enzyme in c("ADAR", "APOBEC")){
  
  dir.create(recursive = TRUE, paste0("./concatenated_samples/", enzyme, "/REDItools.VarScan/tables"))
  dir.create(recursive = TRUE, paste0("./concatenated_samples/", enzyme, "/REDItools.VarScan/xlsx"))
  
  path2 = paste0("./concatenated_samples/", enzyme, "/REDItools.VarScan/tables")
  path3 = paste0("./concatenated_samples/", enzyme, "/REDItools.VarScan/xlsx")
  
  if (enzyme == "ADAR"){
    ind_end = "RADAR"
  } else if(enzyme == "APOBEC"){
    if(genome == "mouse"){
      ind_end = "known_C_to_U"
    } else {
      ind_end = c()
    }
  }
  
  ind1 = grep(enzyme, list.files("./renamed_files"))
  
  path1 = paste0(list.files(full.names = TRUE,"./renamed_files")[ind1],
                 "/tables")
  
  ind_sort = c()
  for(i2 in phenotypes){
    for(i3 in stages){
      ind_sort = c(ind_sort, grep(paste0(i2, ".",i3), list.files(path1)))
    }
  }
  
  table_names = list.files(full.names=TRUE, path1)[ind_sort]
  
  index = c()
  for(j in 1:length(strsplit(table_names,"[.]"))){
    a=strsplit(table_names,"[.]")[[j]]
    index = c(index, paste0(".",a[length(a)-2],".",a[length(a)-1],".",a[length(a)]))
  }
  
  output = "."
  for(j1 in 1:length(phenotypes)){
    if(j1 == 1){
      output = paste0(output, phenotypes[j1])
    } else {
      output = paste(output, phenotypes[j1], sep = ".")
    }
  }
  for(j2 in 1:length(stages)){
    if(j1 == 1){
      output = paste0(output, stages[j2])
    } else {
      output = paste(output, stages[j2], sep = ".")
    }
  }
  output = paste0(enzyme,output)
  
  concatenate_tables(no.tables = length(table_names),
                     table_names = table_names,
                     index = index,
                     output = output,
                     both.tools = TRUE)
  
}

# End: Concatenate all samples of both tools

# Start:
# Common positions and t.test (dual comparisons) for samples with the same
# tool and stage, and different phenotype

library(genefilter)
library(dplyr)
library(xlsx)

for(d1 in c("ADAR","APOBEC")){
  for(d2 in c(tools, "REDItools.VarScan")){
    for(d3 in c("tables","xlsx")){
      dir.create(recursive = TRUE, paste0("./dual_phenotype_comparisons/",
                                          d1, "/", d2, "/", d3))
    }
  }
}

for(enzyme in c("ADAR", "APOBEC")){
  for(d in c(tools, "REDItools.VarScan")){
    
    path1 = paste0("./concatenated_samples/",enzyme,"/",d,"/tables")
    
    if (d == "REDItools.VarScan"){
      tools_2 = tools
    } else if(d == "REDItools"){
      tools_2 = "REDItools"
    } else if (d == "VarScan"){
      tools_2 = "VarScan"
    }
    for( i1 in tools_2){
      for(i2 in stages){
        t = read.table(header = TRUE, list.files(full.names = TRUE,
                                                 path1))
        if(enzyme == "ADAR"){
          col_end = which(colnames(t)=="RADAR")
        } else if (enzyme == "APOBEC"){
          if(genome == "mouse"){
            col_end = which(colnames(t)=="known_C_to_U")
          } else if (genome == "human"){
            col_end = which(colnames(t)=="repFamily.2")
          }
        }
        
        t %>%
          (function(x){
            if(d == "REDItools.VarScan"){
              x[,c(1:col_end, grep(i1, colnames(x)))]
            } else {x}
          }) %>%
          (function(x){
            x[,c(1:col_end, grep(i2, colnames(x)))]
          }) %>%
          (function(x){
            for(j1 in 1:length(dual_phenotypes)){
              f = c()
              ind_dual_phenotypes = c()
              for(j2 in dual_phenotypes[[j1]]){
                ind_dual_phenotypes = c(ind_dual_phenotypes,grep(j2, colnames(x)))
                f = c(f, rep((length(f)+1), length(grep(j2, colnames(x)))))
              }
              x[,c(1:col_end.1,ind_dual_phenotypes)] %>%
                (function(y){
                  mm = as.matrix(y[,grep("Variant_allele_frequency", colnames(y))])
                  mutate(y, p_val = as.numeric(sub(NaN, 1, rowttests(mm, fac = as.factor(f))$p.val)))
                }) %>%
                (function(y){
                  
                  if(d == "REDItools.VarScan"){
                    
                    ind1 = paste0(".",i1,
                                  ".", dual_phenotypes[[j1]][1],
                                  ".", i2)
                    ind2 = paste0(".",i1,
                                  ".", dual_phenotypes[[j1]][2],
                                  ".", i2)
                    ind3 = paste0(".",i1,
                                  ".", dual_phenotypes[[j1]][1],
                                  ".", dual_phenotypes[[j1]][2],
                                  ".", i2)
                    
                  } else {
                    
                    ind1 = paste0(".", dual_phenotypes[[j1]][1],
                                  ".", i2)
                    
                    ind2 = paste0(".", dual_phenotypes[[j1]][2],
                                  ".", i2)
                    
                    ind3 = paste0(".", dual_phenotypes[[j1]][1],
                                  ".", dual_phenotypes[[j1]][2],
                                  ".", i2)
                  }
                  
                  m1 = rowMeans(y[,grep(ind1, colnames(y))])
                  m2 = rowMeans(y[,grep(ind2, colnames(y))])
                  st1 = apply(y[,grep(ind1, colnames(y))], 1, sd)
                  st2 = apply(y[,grep(ind2, colnames(y))], 1, sd)
                  
                  h = cbind(y[,1:col_end.1], 
                            y[,grep(ind1, colnames(y))],
                            as.data.frame(m1),
                            as.data.frame(st1),
                            y[,grep(ind2, colnames(y))],
                            as.data.frame(m2),
                            as.data.frame(st2),
                            y[,which(colnames(y)=="p_val")])
                  
                  colnames(h)[which(colnames(h) == "m1")] = paste0("Mean_Var_freq",ind1)
                  colnames(h)[which(colnames(h) == "m2")] = paste0("Mean_Var_freq",ind2)
                  colnames(h)[which(colnames(h) == "st1")] = paste0("stdev",ind1)
                  colnames(h)[which(colnames(h) == "st2")] = paste0("stdev",ind2)
                  colnames(h)[grep("p_val", colnames(h))] = paste0("p_val",ind3)
                  
                  h = h[order(h$p_val),]
                  
                  h = h[!rowSums(h[,grep("Variant_allele_frequency",colnames(h))]) == 0,]
                  
                  write.table(row.names = FALSE, h, file = paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/tables/",enzyme,ind3))
                  write.xlsx(row.names = FALSE, h, file = paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/xlsx/",enzyme,paste0(ind3, ".xlsx")))
                  
                  # Venn plot a. general and b. for p_val <= 0.05
                  # a.
                  l_phenotype.1 = h[!h[,grep(paste0("Mean_Var_freq",ind1),colnames(h))] == 0,] %>%
                    (function(k){
                      l1=c()
                      for(j4 in seq(nrow(k))){
                        l1=c(l1,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l1
                    })
                  l_phenotype.2 = h[!h[,grep(paste0("Mean_Var_freq",ind2),colnames(h))] == 0,] %>%
                    (function(k){
                      l2=c()
                      for(j4 in seq(nrow(k))){
                        l2=c(l2,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l2
                    })
                  overlap.1 = calculate.overlap(
                    x = list(l_phenotype.1,l_phenotype.2)
                  )
                  venn.plot=draw.pairwise.venn(length(overlap.1[[1]]),
                                               length(overlap.1[[2]]), 
                                               length(overlap.1[[3]]),
                                               lty = "blank",
                                               fill = c("blue","red"),
                                               category = c(paste0(dual_phenotypes[[j1]][1],".", i2),
                                                            paste0(dual_phenotypes[[j1]][2],".", i2)),
                                               cex=1.5);dev.copy(png,paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/xlsx/Venn_plot",ind3,".png"));dev.off();dev.off()
                  
                  # b.
                  h2 = h[h[,grep("p_val",colnames(h))]<=0.05,]
                  
                  write.table(row.names = FALSE, h2, file = paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/tables/",enzyme,"_below_0.05_p_val",ind3))
                  write.xlsx(row.names = FALSE, h2, file = paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/xlsx/",enzyme,"_below_0.05_p_val",paste0(ind3, ".xlsx")))
                  
                  l_phenotype.1 = h2[!h2[,grep(paste0("Mean_Var_freq",ind1),colnames(h2))] == 0,] %>%
                    (function(k){
                      l1=c()
                      for(j4 in seq(nrow(k))){
                        l1=c(l1,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l1
                    })
                  l_phenotype.2 = h2[!h2[,grep(paste0("Mean_Var_freq",ind2),colnames(h2))] == 0,] %>%
                    (function(k){
                      l2=c()
                      for(j4 in seq(nrow(k))){
                        l2=c(l2,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l2
                    })
                  overlap.1 = calculate.overlap(
                    x = list(l_phenotype.1,l_phenotype.2)
                  )
                  venn.plot=draw.pairwise.venn(length(overlap.1[[1]]),
                                               length(overlap.1[[2]]), 
                                               length(overlap.1[[3]]),
                                               lty = "blank",
                                               fill = c("blue","red"),
                                               category = c(paste0(dual_phenotypes[[j1]][1],".", i2),
                                                            paste0(dual_phenotypes[[j1]][2],".", i2)),
                                               cex=1.5);dev.copy(png,paste0("./dual_phenotype_comparisons/",enzyme,"/",d,"/xlsx/Venn_plot_below_0.05_p_val",ind3,".png"));dev.off();dev.off()
                  
                })
            }
          })
      }
    }
  }
}

# End:
# Common positions and t.test (dual comparisons) for samples with the same
# tool and stage, and different phenotype

# Start:
# Common positions and t.test (dual comparisons) for samples with the same
# tool and phenotype, and different stage

library(genefilter)
library(dplyr)
library(xlsx)

for(d1 in c("ADAR","APOBEC")){
  for(d2 in c(tools, "REDItools.VarScan")){
    for(d3 in c("tables","xlsx")){
      dir.create(recursive = TRUE, paste0("./dual_stage_comparisons/",
                                          d1, "/", d2, "/", d3))
    }
  }
}

for(enzyme in c("ADAR", "APOBEC")){
  for(d in c(tools, "REDItools.VarScan")){
    
    path1 = paste0("./concatenated_samples/",enzyme,"/",d,"/tables")
    
    if (d == "REDItools.VarScan"){
      tools_2 = tools
    } else if(d == "REDItools"){
      tools_2 = "REDItools"
    } else if (d == "VarScan"){
      tools_2 = "VarScan"
    }
    for( i1 in tools_2){
      for(i2 in phenotypes){
        t = read.table(header = TRUE, list.files(full.names = TRUE,
                                                 path1))
        if(enzyme == "ADAR"){
          col_end = which(colnames(t)=="RADAR")
        } else if (enzyme == "APOBEC"){
          if(genome == "mouse"){
            col_end = which(colnames(t)=="known_C_to_U")
          } else if (genome == "human"){
            col_end = which(colnames(t)=="repFamily.2")
          }
        }
        
        t %>%
          (function(x){
            if(d == "REDItools.VarScan"){
              x[,c(1:col_end, grep(i1, colnames(x)))]
            } else {x}
          }) %>%
          (function(x){
            x[,c(1:col_end, grep(i2, colnames(x)))]
          }) %>%
          (function(x){
            for(j1 in 1:length(dual_stages)){
              f = c()
              ind_dual_stages = c()
              for(j2 in dual_stages[[j1]]){
                ind_dual_stages = c(ind_dual_stages,grep(j2, colnames(x)))
                f = c(f, rep((length(f)+1), length(grep(j2, colnames(x)))))
              }
              x[,c(1:col_end.1,ind_dual_stages)] %>%
                (function(y){
                  mm = as.matrix(y[,grep("Variant_allele_frequency", colnames(y))])
                  mutate(y, p_val = as.numeric(sub(NaN, 1, rowttests(mm, fac = as.factor(f))$p.val)))
                }) %>%
                (function(y){
                  
                  if(d == "REDItools.VarScan"){
                    
                    ind1 = paste0(".",i1, ".", i2,
                                  ".", dual_stages[[j1]][1]
                    )
                    ind2 = paste0(".",i1, ".", i2,
                                  ".", dual_stages[[j1]][2]
                    )
                    ind3 = paste0(".",i1, ".", i2,
                                  ".", dual_stages[[j1]][1],
                                  ".", dual_stages[[j1]][2]
                    )
                    
                  } else {
                    
                    ind1 = paste0(".", i2,
                                  ".", dual_stages[[j1]][1]
                    )
                    ind2 = paste0(".", i2,
                                  ".", dual_stages[[j1]][2]
                    )
                    ind3 = paste0(".", i2,
                                  ".", dual_stages[[j1]][1],
                                  ".", dual_stages[[j1]][2]
                    )
                  }
                  
                  m1 = rowMeans(y[,grep(ind1, colnames(y))])
                  m2 = rowMeans(y[,grep(ind2, colnames(y))])
                  st1 = apply(y[,grep(ind1, colnames(y))], 1, sd)
                  st2 = apply(y[,grep(ind2, colnames(y))], 1, sd)
                  
                  h = cbind(y[,1:col_end.1], 
                            y[,grep(ind1, colnames(y))],
                            as.data.frame(m1),
                            as.data.frame(st1),
                            y[,grep(ind2, colnames(y))],
                            as.data.frame(m2),
                            as.data.frame(st2),
                            y[,which(colnames(y)=="p_val")])
                  
                  colnames(h)[which(colnames(h) == "m1")] = paste0("Mean_Var_freq",ind1)
                  colnames(h)[which(colnames(h) == "m2")] = paste0("Mean_Var_freq",ind2)
                  colnames(h)[which(colnames(h) == "st1")] = paste0("stdev",ind1)
                  colnames(h)[which(colnames(h) == "st2")] = paste0("stdev",ind2)
                  colnames(h)[grep("p_val", colnames(h))] = paste0("p_val",ind3)
                  
                  h = h[order(h$p_val),]
                  
                  h = h[!rowSums(h[,grep("Variant_allele_frequency",colnames(h))]) == 0,]
                  
                  write.table(row.names = FALSE, h, file = paste0("./dual_stage_comparisons/",enzyme,"/",d,"/tables/",enzyme,ind3))
                  write.xlsx(row.names = FALSE, h, file = paste0("./dual_stage_comparisons/",enzyme,"/",d,"/xlsx/",enzyme,paste0(ind3, ".xlsx")))
                  
                  # Venn plot a. general and b. for p_val <= 0.05
                  # a.
                  l_stage.1 = h[!h[,grep(paste0("Mean_Var_freq",ind1),colnames(h))] == 0,] %>%
                    (function(k){
                      l1=c()
                      for(j4 in seq(nrow(k))){
                        l1=c(l1,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l1
                    })
                  l_stage.2 = h[!h[,grep(paste0("Mean_Var_freq",ind2),colnames(h))] == 0,] %>%
                    (function(k){
                      l2=c()
                      for(j4 in seq(nrow(k))){
                        l2=c(l2,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l2
                    })
                  overlap.1 = calculate.overlap(
                    x = list(l_stage.1,l_stage.2)
                  )
                  venn.plot=draw.pairwise.venn(length(overlap.1[[1]]),
                                               length(overlap.1[[2]]), 
                                               length(overlap.1[[3]]),
                                               lty = "blank",
                                               fill = c("blue","red"),
                                               category = c(paste0(dual_stages[[j1]][1],".", i2),
                                                            paste0(dual_stages[[j1]][2],".", i2)),
                                               cex=1.5);dev.copy(png,paste0("./dual_stage_comparisons/",enzyme,"/",d,"/xlsx/Venn_plot",ind3,".png"));dev.off();dev.off()
                  
                  # b.
                  h2 = h[h[,grep("p_val",colnames(h))]<=0.05,]
                  
                  write.table(row.names = FALSE, h2, file = paste0("./dual_stage_comparisons/",enzyme,"/",d,"/tables/",enzyme,"_below_0.05_p_val",ind3))
                  write.xlsx(row.names = FALSE, h2, file = paste0("./dual_stage_comparisons/",enzyme,"/",d,"/xlsx/",enzyme,"_below_0.05_p_val",paste0(ind3, ".xlsx")))
                  
                  l_stage.1 = h2[!h2[,grep(paste0("Mean_Var_freq",ind1),colnames(h2))] == 0,] %>%
                    (function(k){
                      l1=c()
                      for(j4 in seq(nrow(k))){
                        l1=c(l1,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l1
                    })
                  l_stage.2 = h2[!h2[,grep(paste0("Mean_Var_freq",ind2),colnames(h2))] == 0,] %>%
                    (function(k){
                      l2=c()
                      for(j4 in seq(nrow(k))){
                        l2=c(l2,paste0(k$Chr[j4],"-",k$Start[j4]))
                      }
                      l2
                    })
                  overlap.1 = calculate.overlap(
                    x = list(l_stage.1,l_stage.2)
                  )
                  venn.plot=draw.pairwise.venn(length(overlap.1[[1]]),
                                               length(overlap.1[[2]]), 
                                               length(overlap.1[[3]]),
                                               lty = "blank",
                                               fill = c("blue","red"),
                                               category = c(paste0(dual_stages[[j1]][1],".", i2),
                                                            paste0(dual_stages[[j1]][2],".", i2)),
                                               cex=1.5);dev.copy(png,paste0("./dual_stage_comparisons/",enzyme,"/",d,"/xlsx/Venn_plot_below_0.05_p_val",ind3,".png"));dev.off();dev.off()
                  
                })
            }
          })
      }
    }
  }
}

# End:
# Common positions and t.test (dual comparisons) for samples with the same
# tool and phenotype, and different stage

# Venn_plots for total positions by both tools in separate script
# ggplot - barplots for individual tools in separate script

# Start:
# ggplot - barplots for position distribution between different regions
# and chisq.test to check if position distribution dependents upon 
# phenotype and stage

library(dplyr)
library(ggplot2)
library(xlsx)
library(reshape)
library(scales)

for(enzyme in c("ADAR","APOBEC")){
  for(d in c(tools,"REDItools.VarScan")){
    
    path1 = paste0("./concatenated_samples/",enzyme, "/",d, "/tables")
    
    t1 = read.table(header = TRUE ,list.files(full.names = TRUE, path1))
    
    # Different phenotype dual comparison
    
    for(i1 in stages){
      t2 = t1 %>%
        (function(x){
          x[,grep("Variant_allele_frequency", colnames(x))]
        }) %>%
        (function(x){
          x[,grep(i1, colnames(x))]
        }) %>%
        (function(x){
          for(j1 in 1:length(dual_phenotypes)){
            
            s = as.list(seq(length(dual_phenotypes[[j1]])))
            
            for(j2 in 1:length(dual_phenotypes[[j1]])){
              
              t3 = x[,grep(dual_phenotypes[[j1]][j2], colnames(x))]
              
              t4 = as.data.frame(table(t1[!rowSums(t3) == 0,]$Func.refGene))
              
              t4=cbind(t4,rep(paste0(d,".",
                                     dual_phenotypes[[j1]][j2],".",
                                     i1), nrow(t4)))
              
              colnames(t4)=c("Region","Positions","Sample")
              
              s[[j2]]=t4
            }
            
            for(j2 in 1:length(dual_phenotypes[[j1]])){
              ind = s[[j2]]$Region %in% c("exonic",
                                          "UTR3",
                                          "UTR5",
                                          "intronic",
                                          "intergenic",
                                          "ncRNA_exonic",
                                          "ncRNA_intronic")
              
              t5=data.frame(Region = "other", 
                            Positions = sum(s[[j2]][!ind,]$Positions),
                            Sample = paste0(d,".",
                                            dual_phenotypes[[j1]][j2],".",
                                            i1))
              s[[j2]] = rbind(s[[j2]][ind,],t5)
            }
            
            # chi-squared contingency table test
            
            d3_names = c()
            
            for(j2 in 1:length(dual_phenotypes[[j1]])){
              if(j2 == 1){
                d3 = merge(s[[j2]][,1:2],s[[j2+1]][,1:2],
                           by.x = colnames(s[[j2]])[1],
                           by.y = colnames(s[[j2+1]])[1])
              }
              if(j2 == 2){print("skipped 2nd")}
              if(j2 > 2){
                d3 = merge(d3,s[[j2]][,1:2],
                           by.x = colnames(d3)[1],
                           by.y = colnames(s[[j2]])[1])
              }
              d3_names = c(d3_names, paste0(d,".",
                                            dual_phenotypes[[j1]][j2],".",
                                            i1))
            }
            
            rownames(d3) = d3[,1]
            
            d3 = d3[,2:ncol(d3)]
            
            colnames(d3) = d3_names
            
            chisq.test(d3)$p.val
            
            write.xlsx(d3, file = paste0("./dual_phenotype_comparisons/",
                                         enzyme,"/",d,"/xlsx/","Chisq.table_",chisq.test(d3)$p.val,"_",enzyme,".",d,".",
                                         dual_phenotypes[[j1]][j2],".",
                                         i1,".xlsx"))
            
            # ggplot - barplots
            
            for(j2 in 1:length(dual_phenotypes[[j1]])){
              if(j2 == 1){
                d4=rbind(s[[j2]],s[[j2+1]])
              }
              if(j2 == 2){print("skipped 2nd")}
              if(j2 > 2){
                d4=rbind(d4,s[[j2]])
              }
            }
            
            plot_1 = d4 %>%
              ggplot() +
              aes(x = Sample, y = Positions,
                  fill = Region) +
              theme(plot.title = element_text(lineheight=.8, face="bold"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1)) +
              geom_bar(position = "fill", stat = "identity") +
              scale_y_continuous(labels = percent)
            
            print(plot_1);dev.copy(png,paste0("./dual_phenotype_comparisons/",
                                              enzyme,"/",d,"/xlsx/","Percentage_barplot_",enzyme,".",d,".",
                                              dual_phenotypes[[j1]][j2],".",
                                              i1,".png"));dev.off();dev.off()
            
            plot_2 = d4 %>%
              ggplot() +
              aes(x = Sample, y = Positions,
                  fill = Region) +
              theme(plot.title = element_text(lineheight=.8, face="bold"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1)) +
              geom_bar(stat = "identity")
            
            print(plot_2);dev.copy(png,paste0("./dual_phenotype_comparisons/",
                                              enzyme,"/",d,"/xlsx/","Positions_barplot_",enzyme,".",d,".",
                                              dual_phenotypes[[j1]][j2],".",
                                              i1,".png"));dev.off();dev.off()
            
          }
        })
    }
    
    # Different stage dual comparison
    
    for(i1 in phenotypes){
      t2 = t1 %>%
        (function(x){
          x[,grep("Variant_allele_frequency", colnames(x))]
        }) %>%
        (function(x){
          x[,grep(i1, colnames(x))]
        }) %>%
        (function(x){
          for(j1 in 1:length(dual_stages)){
            
            s = as.list(seq(length(dual_stages[[j1]])))
            
            for(j2 in 1:length(dual_stages[[j1]])){
              
              t3 = x[,grep(dual_stages[[j1]][j2], colnames(x))]
              
              t4 = as.data.frame(table(t1[!rowSums(t3) == 0,]$Func.refGene))
              
              t4=cbind(t4,rep(paste0(d,".", 
                                     i1, ".",
                                     dual_stages[[j1]][j2]), nrow(t4)))
              
              colnames(t4)=c("Region","Positions","Sample")
              
              s[[j2]]=t4
            }
            
            for(j2 in 1:length(dual_stages[[j1]])){
              ind = s[[j2]]$Region %in% c("exonic",
                                          "UTR3",
                                          "UTR5",
                                          "intronic",
                                          "intergenic",
                                          "ncRNA_exonic",
                                          "ncRNA_intronic")
              
              t5=data.frame(Region = "other", 
                            Positions = sum(s[[j2]][!ind,]$Positions),
                            Sample = paste0(d,".",
                                            i1,".",
                                            dual_stages[[j1]][j2]))
              s[[j2]] = rbind(s[[j2]][ind,],t5)
            }
            
            # chi-squared contingency table test
            
            d3_names = c()
            
            for(j2 in 1:length(dual_stages[[j1]])){
              if(j2 == 1){
                d3 = merge(s[[j2]][,1:2],s[[j2+1]][,1:2],
                           by.x = colnames(s[[j2]])[1],
                           by.y = colnames(s[[j2+1]])[1])
              }
              if(j2 == 2){print("skipped 2nd")}
              if(j2 > 2){
                d3 = merge(d3,s[[j2]][,1:2],
                           by.x = colnames(d3)[1],
                           by.y = colnames(s[[j2]])[1])
              }
              d3_names = c(d3_names, paste0(d,".",i1,".",
                                            dual_stages[[j1]][j2]))
            }
            
            rownames(d3) = d3[,1]
            
            d3 = d3[,2:ncol(d3)]
            
            colnames(d3) = d3_names
            
            chisq.test(d3)$p.val
            
            write.xlsx(d3, file = paste0("./dual_stage_comparisons/",
                                         enzyme,"/",d,"/xlsx/","Chisq.table_",chisq.test(d3)$p.val,"_",enzyme,".",d,".",
                                         i1,".",dual_stages[[j1]][j2],".xlsx"))
            
            # ggplot - barplots
            
            for(j2 in 1:length(dual_stages[[j1]])){
              if(j2 == 1){
                d4=rbind(s[[j2]],s[[j2+1]])
              }
              if(j2 == 2){print("skipped 2nd")}
              if(j2 > 2){
                d4=rbind(d4,s[[j2]])
              }
            }
            
            plot_1 = d4 %>%
              ggplot() +
              aes(x = Sample, y = Positions,
                  fill = Region) +
              theme(plot.title = element_text(lineheight=.8, face="bold"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1)) +
              geom_bar(position = "fill", stat = "identity") +
              scale_y_continuous(labels = percent)
            
            print(plot_1);dev.copy(png,paste0("./dual_stage_comparisons/",
                                              enzyme,"/",d,"/xlsx/","Percentage_barplot_",enzyme,".",d,".",
                                              i1,".",dual_stages[[j1]][j2],
                                              ".png"));dev.off();dev.off()
            
            plot_2 = d4 %>%
              ggplot() +
              aes(x = Sample, y = Positions,
                  fill = Region) +
              theme(plot.title = element_text(lineheight=.8, face="bold"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1)) +
              geom_bar(stat = "identity")
            
            print(plot_2);dev.copy(png,paste0("./dual_stage_comparisons/",
                                              enzyme,"/",d,"/xlsx/","Positions_barplot_",enzyme,".",d,".",
                                              i1,".",dual_stages[[j1]][j2],
                                              ".png"));dev.off();dev.off()
            
          }
        })
    }
    
  }
}
