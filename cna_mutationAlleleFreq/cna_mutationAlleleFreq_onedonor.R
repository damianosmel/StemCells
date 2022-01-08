###
# Code to calculate CNA (copy number alternations) and mutation allele frequency 
# for all mutations of a donor
###

library(cgdsr)
library(readxl)
# connect to cBioPortal
cgds <- CGDS("https://www.cbioportal.org/")
test(cgds)

setwd("/home/damian/Documents/L3S/projects/stemCells/maike/analysis/code")

#get list of curated studies
curated_studies_file = "cna_mutationAlleleFreq/curated-study-list-20190522.txt"

get_curated_studies <- function(file_path){
  cat("Reading list of curated studies, please wait..\n")
  curated_studies <- readLines(file_path)
  cat("Success: ",length(curated_studies)," studies are loaded.\n")
  return(curated_studies)
}
curated_studies <- get_curated_studies(curated_studies_file)

#get mutated genes
#test
mutated_genes <- c("HPCAL1","HNRNPF","PPP1CC","REXO1","FBXW9","JAK1","SLC25A44","CEP350","ERC2","CLDN4","KPNA7","TMEM168","NAIF1","MYO7A","SCEL","HIRIP3","PLA2G3","XIAP")
mutated_genes

#create names of allele frequency for each mutated gene
create_genes_allelefreq_names <- function(mutated_genes){
  genes_allele_freq <- c()
  for(mutated_gene in mutated_genes){
    genes_allele_freq <- c(genes_allele_freq,paste0(mutated_gene, "_allele_freq"))
  }
  cat("Sucess: list of ",length(genes_allele_freq)," mutated genes ared appended with allele_frequency.\n")
  return(genes_allele_freq)
}

mutated_genes_allele_freq <- create_genes_allelefreq_names(mutated_genes)
mutated_genes_allele_freq

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

cat("\n Cna_mut_ids:")
cna_mut_ids_in_studies

###
# cna + allele freq
###
total_rows_count <- 0
for(study_index in c(1:length(curated_studies))){
  current_rows_count <- 0
  cat("\n--- ", curated_studies[study_index]," ---\n")
  cat("ids: ", cna_mut_ids_in_studies[study_index,],"\n")
  all_cases <- getCaseLists(cgds, curated_studies[study_index])[1,1]
  #get mutation data
  #cat("\nMutations: \n")
  if(cna_mut_ids_in_studies[study_index,2] !="notAvailable"){
    mutation_profile <- getMutationData(cgds,all_cases,cna_mut_ids_in_studies[study_index,2],mutated_genes)
    #create mutations (+allele frequency) dataframe
    mutation_profile$allele_freq <- mutation_profile$variant_read_count_tumor/(mutation_profile$reference_read_count_tumor+mutation_profile$variant_read_count_tumor)
    #print(mutation_profile)
  }else{
    cat("Current study does not contain mutations.")
  }
  
  #get cna data
  #cat("\nCNA: \n")
  if(cna_mut_ids_in_studies[study_index,1] !="notAvailable"){
    genomic_profile <- getProfileData(cgds,mutated_genes, cna_mut_ids_in_studies[study_index,1],all_cases)
    genomic_profile <- cbind(Studies=curated_studies[study_index],genomic_profile)
    #append genomic profile with columns for allele frequency
    for(gene_allele_freq in mutated_genes_allele_freq){
      genomic_profile[gene_allele_freq] <- NA
    }
    #print(genomic_profile)
    current_rows_count <- current_rows_count + nrow(genomic_profile)
  }else{
    cat("CNA info not available\n")
    #create empty data frame if cna info is not available
    column_names <- c("Studies",mutated_genes,mutated_genes_allele_freq)
    genomic_profile <- setNames(data.frame(matrix(ncol = length(mutated_genes)+length(mutated_genes_allele_freq)+1, nrow = 0)), column_names)
  }
  
  #write allele freq information into cna genomic profile
  cat("\nAppend CNA with allele frequency:\n")
  if(cna_mut_ids_in_studies[study_index,2] != "notAvailable" & nrow(mutation_profile)>0){ #covered case no mutation profile available
    #loop over mutation_profile data
    #check if case is in genomic_profile
    #if yes, update genomic_profile column for the specific mutation
    for(mutation in 1:nrow(mutation_profile)){
      #print(mutation_profile[mutation,])
      gene <- mutation_profile[mutation,"gene_symbol"]
      cat("Gene: ",gene,"\n")
      case <- mutation_profile[mutation,"case_id"]
      case <- gsub("-",".",case)
      cat("Case: ",case,"\n")
      if(case %in% row.names(genomic_profile) & cna_mut_ids_in_studies[study_index,1] !="notAvailable"){
        #cat(case,"in row\n")
        allele_freq_gene_column <- paste0(gene,"_allele_freq")
        #cat("Old value: ", genomic_profile[case,allele_freq_gene_column],"\n")
        genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
        #cat("New value: ", genomic_profile[case,allele_freq_gene_column],"\n")
      }else if(cna_mut_ids_in_studies[study_index,1] =="notAvailable"){
        cat("Not available cna: update allele freq\n")
        allele_freq_gene_column <- paste0(gene,"_allele_freq")
        if(case %in% rownames(genomic_profile)){
          #cat(case," exists, update row\n")
          genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
        }else{
          genomic_profile[case,] <-c(curated_studies[study_index],rep(NaN,times=length(mutated_genes)),rep(NA,times=length(mutated_genes_allele_freq)))
          #cat(case,"is new, create row!\n")
          #cat("Old value: ", genomic_profile[case,allele_freq_gene_column],"\n")
          genomic_profile[case,allele_freq_gene_column] <- mutation_profile[mutation,"allele_freq"]
          #cat("New value: ", genomic_profile[case,allele_freq_gene_column],"\n")
          current_rows_count <- current_rows_count + 1  
        }
      }
    }
    cat("---\n")
  }else{
    cat(curated_studies[study_index],"does not provide or contains zero mutation information.\n")
  }
  
  cat("\n### ###\n")
  if(study_index==1){
    cat("Creating cna_mut_studies.\n")
    cna_mut_studies <- genomic_profile
  }else{
    cat("Appending cna_mut_studies.\n")
    cna_mut_studies <- rbind(cna_mut_studies,genomic_profile)
  }
  total_rows_count <- total_rows_count + current_rows_count
  
  if(nrow(cna_mut_studies) != total_rows_count){
    cat(curated_studies[study_index],"creates more rows than computed.\n")
    cat("test: break!")
    break
  }
}

cat("CNA and allele frequency data frame for all studies:\n")
head(cna_mut_studies)
cat("\nNumber of total rows:",total_rows_count,"\n")
write.table(data.frame("Samples"=rownames(cna_mut_studies),cna_mut_studies),"cna_mutationAlleleFreq/hSVEC38_clones23.csv", row.names=FALSE,sep=",")