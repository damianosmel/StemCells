library(cgdsr)
library(readxl)

###
# Code to calculate CNA (copy number alternations) and mutation allele frequency 
# for all mutations of each clone of a donor in xlsx file
###

cat("\n###\n")
cat("# Step 0: Connect to cBioPortal\n")
cat("###\n")
# connect to cBioPortal
cgds <- CGDS("https://www.cbioportal.org/")
test(cgds)

setwd("/home/damian/Documents/L3S/projects/stemCells/maike/analysis/code")

cat("\n###\n")
cat("# Step 1: Read curated studies.\n")
cat("###\n")
#get list of curated studies
curated_studies_file = "cna_mutationAlleleFreq/curated-study-list-20190522.txt"

get_curated_studies <- function(file_path){
  cat("Reading list of curated studies, please wait..\n")
  curated_studies <- readLines(file_path)
  cat("Success: ",length(curated_studies)," studies are loaded.\n")
  return(curated_studies)
}
curated_studies <- get_curated_studies(curated_studies_file)

cat("\n###\n")
cat("# Step 2: Save ids for cna + mutation information ids.\n")
cat("###\n")
#for each study find cna and mutation information ids
#todo: make the matrix dataframe with unterstable column names
cna_mut_ids_in_studies <- matrix("notAvailable",length(curated_studies),2)
cat(cna_mut_ids_in_studies,"\n")
study_index <- 1
for(curated_study in curated_studies){
  cat("\n--- ", curated_study, " ---\n")
  genetic_profile_study <- getGeneticProfiles(cgds,curated_study)
  if(nrow(genetic_profile_study)==0){
    cat("Study",curated_study,"does not contain CNA and mutation info.\n")
    study_index <- study_index + 1
    next
  }
  profile_ids <- genetic_profile_study[,"genetic_profile_id"]
  profile_names <- genetic_profile_study[,"genetic_profile_name"]
  cat("Profile names:\n")
  cat(profile_names,"\n")
  cna_index <- grep("GISTIC",profile_names,ignore.case=TRUE,value=FALSE)
  mut_index <- grep("Mutations",profile_names,ignore.case=TRUE,value=FALSE)
  if(length(cna_index)==0){
    cat("Study",curated_study,"does not contain CNA info.\n")
  }
  else{
    cna_index <- cna_index
    cat("Study",curated_study,"contains CNA at ", cna_index,"\n")
    cna_mut_ids_in_studies[study_index,1] <- profile_ids[cna_index]
  }
  if(length(mut_index)==0){
    cat("Study",curated_study,"does not contain mutation info.\n")
  }
  else{
    if(length(mut_index)>1){#get the first occurence of mutations keyword
      mut_index <- mut_index[1]
    }
    cat("Study",curated_study,"contains mutations at ", mut_index,"\n")
    cna_mut_ids_in_studies[study_index,2] <- profile_ids[mut_index]
  }
  study_index <- study_index + 1
}


### Only cancer variants (up July 2019) ###
# cat("\n###\n")
# cat("# Step 3: Parse xlsx and get mutated genes for each clone of each donor\n")
# cat("###\n")
# get_mutated_genes_donor_clone <- function(sheet_name){
#   data_sheet <- read_excel(path="/home/damian/Documents/L3S/projects/stemCells/maike/genes/All_donors_variants_noSlash.xlsx", sheet=sheet_name)
#   donor_clone_mutated_genes <- list()
#   if(sheet_name == "Table_S7a"){
#     cat("Creating mutated genes for each clone of donors in",sheet_name,"\n")
#     donor_clone <- data_sheet$Donor_clone
#     print(donor_clone)
#     unique_donor_clone <- unique(donor_clone)
#     print(unique_donor_clone)
#     for (donor_clone in unique_donor_clone){
#       if(!is.na(donor_clone)){
#         rows_index = data_sheet$Donor_clone == donor_clone
#         mutated_genes <- unique(data_sheet[rows_index,]$Gene_symbol)
#         print(mutated_genes)
#         donor_clone_mutated_genes[[donor_clone]] <- mutated_genes[!is.na(mutated_genes)]
#       }
#       cat("---\n")
#     }
#     
#   }else if(sheet_name == "Table_S7b"){
#     cat("Creating mutated genes for each clone of donors in",sheet_name,"\n")
#     donor_clone <- data_sheet$Donor_and_iPSC_clone
#     print(donor_clone)
#     unique_donor_clone <- unique(donor_clone)
#     print(unique_donor_clone)
#     for(donor_clone in unique_donor_clone){
#       if(!is.na(donor_clone)){
#         rows_index = data_sheet$Donor_and_iPSC_clone == donor_clone
#         mutated_genes <- unique(data_sheet[rows_index,]$Gene_symbol)
#         donor_clone_mutated_genes[[donor_clone]] <- mutated_genes[!is.na(mutated_genes)]
#       }
#       cat("---\n")
#     }
#   }
#   cat("Loaded mutated genes per clone of donor:\n")
#   print(donor_clone_mutated_genes)
#   return(donor_clone_mutated_genes)
# }
# 
# sheet_name <- "Table_S7a"
# mutated_genes_S7a <- get_mutated_genes_donor_clone(sheet_name)
# 
# sheet_name <- "Table_S7b"
# mutated_genes_S7b <- get_mutated_genes_donor_clone(sheet_name)
# 
# #get the total list of mutated genes
# mutated_genes_all <- c(mutated_genes_S7a,mutated_genes_S7b)



### All variants of clone of donor (from Aug. 2019) ###
# cat("\n###\n")
# cat("# Step 3: Parse directory with tabular files and get mutated genes for each clone of each donor.\n")
# cat("###\n")
parse_variant_table <- function(file){
  # credits: https://camnugent.wordpress.com/2018/11/08/returning-multiple-values-from-a-function-in-r/
  table <- read.table(file, header=TRUE,colClasses = c("character","character","integer","character","character","character")) # load file
  print(typeof(table$DONOR_ID))
  return(list(donor_clone=unique(table$DONOR_ID),variants=unique(table$GENE_ID)))
}

parse_all_variants <- function(dir_path){
  ###
  # Read list of tabular files and get list of variant genes for each clone of donor
  ###
  #credits: https://www.r-bloggers.com/looping-through-files/
  donor_clone2vars <- list()
  
  file.names <- dir(dir_path, pattern =".tabular",full.names =TRUE)
  for(i in 1:length(file.names)){
    donor_clone_variants <- parse_variant_table(file.names[i])
    cat("Getting genes for donor_clone: ",donor_clone_variants$donor_clone)
    donor_clone2vars[[donor_clone_variants$donor_clone]] <- donor_clone_variants$variants
  }
  return(donor_clone2vars)
}

#dir_path = "/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants" #correct variants from first time
dir_path = "/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_2check_maike" #corrected by Maike
mutated_genes_all <- parse_all_variants(dir_path)

cat("\n###\n")
cat("# Step 4: Create cna + allele frequency csv for mutated genes on each clone of donor.\n")
cat("###\n")

#create names of allele frequency for each mutated gene
create_genes_allelefreq_names <- function(mutated_genes){
  genes_allele_freq <- c()
  #make all gene names having a minus character (-) to have dot character instead (.)
  #and get the lower case of them
  mutated_genes <- tolower(gsub("-",".",mutated_genes))
  cat("---\n")
  for(mutated_gene in mutated_genes){
    genes_allele_freq <- c(genes_allele_freq,paste0(mutated_gene, "_allele_freq"))
  }
  cat("Success: list of ",length(genes_allele_freq)," mutated genes are appended with allele_frequency.\n")
  return(genes_allele_freq)
}

###
# cna + allele freq
###
donor_clone_index <- 1
for(donor_clone in names(mutated_genes_all)){
  cat("\n----- ----- -----\n")
  cat("Processing:",donor_clone,"\n")
  
  mutated_genes <- sort(mutated_genes_all[[donor_clone]])
  mutated_genes_allele_freq <- create_genes_allelefreq_names(mutated_genes)
  
  total_rows_count <- 0
  for(study_index in c(1:length(curated_studies))){
    current_rows_count <- 0
    cat("\n--- ", curated_studies[study_index]," ---\n")
    cat("ids: ", cna_mut_ids_in_studies[study_index,],"\n")
    cat("study: ",study_index," of ",length(curated_studies),"\n")
    all_cases <- getCaseLists(cgds, curated_studies[study_index])[1,1]
    
    #get mutation data (var_read_count/(ref_read_count+var_read_count))
    cat("\nMutations: \n")
    if(cna_mut_ids_in_studies[study_index,2] !="notAvailable"){
      cat("Mutation available\n")
      mutation_profile <- getMutationData(cgds,all_cases,cna_mut_ids_in_studies[study_index,2],mutated_genes)
      #create mutations (+allele frequency) dataframe
      mutation_profile$allele_freq <- mutation_profile$variant_read_count_tumor/(mutation_profile$reference_read_count_tumor+mutation_profile$variant_read_count_tumor)
    }else{
      cat("Mutation not available\n")
    }
    #cat("mutation_profile names:",names(mutation_profile),"\n")
    
    ### Get CNA data ###
    #out = studies cases x genes
    cat("\nCNA: \n")
    if(cna_mut_ids_in_studies[study_index,1] !="notAvailable"){
      cat("CNA available\n")
      genomic_profile <- getProfileData(cgds,mutated_genes, cna_mut_ids_in_studies[study_index,1],all_cases)
      #cat(names(genomic_profile),"\n")
      #cat("===\n")
      #for a reason the getProfileData returns genes that we didn't query
      #credits: https://www.listendata.com/2015/06/r-keep-drop-columns-from-data-frame.html
      # https://stackoverflow.com/questions/17598134/compare-two-character-vectors-in-r
      # https://stackoverflow.com/questions/25576096/subset-with-vector-specifying-columns-to-drop
      #cat("diff:\n")
      #print(setdiff(names(genomic_profile),mutated_genes))
      #cat("===\n")
      
      genomic_profile <- genomic_profile[, which(names(genomic_profile) %in% mutated_genes)]
      #print(names(genomic_profile))
      #cat("===\n")
      for(gene in mutated_genes){
        if(!(gene %in% names(genomic_profile))){
          genomic_profile[gene] <-NaN
        }
      }
      #sort column names
      genomic_profile <- genomic_profile[ , order(names(genomic_profile))]
      #cat("= sort =\n")
      #print(names(genomic_profile))
      #cat("===\n")
      genomic_profile <- cbind(Studies=curated_studies[study_index],genomic_profile)
      
      #append genomic profile with columns for allele frequency
      for(gene_allele_freq in mutated_genes_allele_freq){
        genomic_profile[gene_allele_freq] <- NaN
      }
      #print(genomic_profile)
      current_rows_count <- current_rows_count + nrow(genomic_profile)
    }else{
      cat("CNA not available\n")
      #create empty data frame if cna info is not available
      column_names <- c("studies",mutated_genes,mutated_genes_allele_freq)
      genomic_profile <- setNames(data.frame(matrix(ncol = length(mutated_genes)+length(mutated_genes_allele_freq)+1, nrow = 0)), column_names)
    }
    #cat("genomic_profile names:",names(genomic_profile),"\n")
    
    ### write allele freq information into cna genomic profile ###
    cat("\nAppend CNA with allele frequency:\n")
    if(cna_mut_ids_in_studies[study_index,2] != "notAvailable" & nrow(mutation_profile)>0){ #covered case no mutation profile available
      #loop over mutation_profile data
      #check if study case is in genomic_profile
      #if yes, update genomic_profile column for the specific mutation
      for(mutation in 1:nrow(mutation_profile)){
        #print(mutation_profile[mutation,])
        gene <- mutation_profile[mutation,"gene_symbol"]
        if(gene %in% mutated_genes){
          gene <- tolower(gsub("-",".",gene))
          #cat("Gene: ",gene,"\n")
          case <- mutation_profile[mutation,"case_id"]
          case <- gsub("-",".",case)
          #cat("Case: ",case,"\n")
          if(case %in% row.names(genomic_profile) & cna_mut_ids_in_studies[study_index,1] !="notAvailable"){
            #cat(case,"in row\n")
            allele_freq_gene_column <- paste0(gene,"_allele_freq")
            #cat("Old value: ", genomic_profile[case,allele_freq_gene_column],"\n")
            genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
            #cat("New value: ", genomic_profile[case,allele_freq_gene_column],"\n")
          }else if(cna_mut_ids_in_studies[study_index,1] =="notAvailable"){
            #cat("Not available cna: update allele freq\n")
            allele_freq_gene_column <- paste0(gene,"_allele_freq")
            if(case %in% rownames(genomic_profile)){
              #cat(case," exists, update row\n")
              genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
            }else{
              genomic_profile[case,] <-c(curated_studies[study_index],rep(NaN,times=length(mutated_genes)),rep(NaN,times=length(mutated_genes_allele_freq)))
              #cat(case,"is new, create row!\n")
              #cat("Old value: ", genomic_profile[case,allele_freq_gene_column],"\n")
              genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
              #cat("New value: ", genomic_profile[case,allele_freq_gene_column],"\n")
              current_rows_count <- current_rows_count + 1  
            }
        }
      }#if gene in mutated genes
      }
      cat("---\n")
    }else{
      cat(curated_studies[study_index],"does not provide or contains zero mutation information.\n")
    }
    
    #normalize resulted data frame to dot as separator character
    colnames(genomic_profile) <- tolower(gsub("-",".",colnames(genomic_profile)))
    #cat("\n### ###\n")
    if(study_index==1){
      cat("Creating cna_mut_studies.\n")
      #cat("genomic_profile names:",names(genomic_profile),"\n")
      cna_mut_studies <- genomic_profile
    }else{
      cat("Appending cna_mut_studies.\n")
      #cat("cna_mut_studies names:",names(cna_mut_studies),"\n")
      #cat("===\n")
      #cat("genomic_profile names:",names(genomic_profile),"\n")
      cna_mut_studies <- rbind(cna_mut_studies,genomic_profile)
      #test
    }
    total_rows_count <- total_rows_count + current_rows_count
    #assert that the size of actual total dataframe equals the sum of cna and mutation information over all studies
    if(nrow(cna_mut_studies) != total_rows_count){
      cat(curated_studies[study_index],"creates more rows than computed.\n")
      cat("assert nrow(cna_mut_studies) != total_rows_count: fail => break!")
      break
    }
  }
  
  cat("CNA and allele frequency data frame for all studies:\n")
  head(cna_mut_studies)
  cat("\nNumber of total rows:",total_rows_count,"\n")
  
  ### Only cancer variants (up July 2019) ###
  #csv_name <- paste0(gsub(", ","_",donor_clone),".csv")
  #csv_name <- paste0("cna_mutationAlleleFreq/",csv_name)
  ### All variants (from Aug 2019) ###
  csv_name <- paste0(donor_clone,".csv")
  #csv_name <- paste0("/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_cna_mutAllFreq/",csv_name) #correct from first time
  csv_name <- paste0("/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_2check_maike_cna_mutAllFreq/",csv_name) #corrected by Maike
  write.table(data.frame("samples"=rownames(cna_mut_studies),cna_mut_studies),csv_name, row.names=FALSE,sep=",")
  cat("CNA + allele frequency for mutations are written in",csv_name,"\n")
  cat("----- ----- -----\n")
  
  #remove current genomic_profile and update donor_clone index
  rm(genomic_profile,mutation_profile,cna_mut_studies)
  donor_clone_index <- donor_clone_index + 1
}

