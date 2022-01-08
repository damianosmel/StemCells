###
# Read xlsx with all donors and get list of mutated genes for each clone
# time reference: the input xlsx have been the first mutation files ever processed.
###

library(readxl)
# is.na(): http://www.cookbook-r.com/Manipulating_data/Comparing_vectors_or_factors_with_NA/
# list: https://www.datamentor.io/r-programming/list/
# indexing: http://www.r-tutor.com/r-introduction/data-frame/data-frame-row-slice
# gsub: https://stackoverflow.com/questions/11936339/replace-specific-characters-within-strings
data_S7a <- read_excel(path="/home/damian/Documents/L3S/projects/stemCells/maike/genes/All_donors_variants.xlsx", sheet="Table_S7a")
#data_S7b <- read_excel(path="/home/damian/Documents/L3S/projects/stemCells/maike/genes/All_donors_variants.xlsx", sheet="Table_S7b")

donor_clone_S7a <- data_S7a$Donor_clone
unique_donor_clone_S7a <- unique(donor_clone_S7a)

donor_clone_mutated_genes <- list()

###
# loop over each unique combination of donor abd clone
# and save all genes found for that combination.
# credits: remove all NA from the list: 
# https://stackoverflow.com/questions/7706876/remove-na-values-from-a-vector
###
for (donor_clone in unique_donor_clone_S7a){
  if(!is.na(donor_clone)){
    rows_index = data_S7a$Donor_clone == donor_clone
    mutated_genes <- data_S7a[rows_index,]$Gene_symbol
    donor_clone_mutated_genes[[donor_clone]] <- mutated_genes[!is.na(mutated_genes)]
  }
}

get_mutated_genes_donor_clone <- function(sheet_name){
  data_sheet <- read_excel(path="/home/damian/Documents/L3S/projects/stemCells/maike/genes/All_donors_variants_noSlash.xlsx", sheet=sheet_name)
  donor_clone_mutated_genes <- list()
  if(sheet_name == "Table_S7a"){
    cat("Creating mutated genes for each clone of donors in",sheet_name,"\n")
    donor_clone <- data_sheet$Donor_clone
    #print(donor_clone)
    unique_donor_clone <- unique(donor_clone)
    #print(unique_donor_clone)
    for (donor_clone in unique_donor_clone){
      if(!is.na(donor_clone)){
        rows_index = data_sheet$Donor_clone == donor_clone
        mutated_genes <- data_sheet[rows_index,]$Gene_symbol
        donor_clone_mutated_genes[[donor_clone]] <- mutated_genes[!is.na(mutated_genes)]
      }
    }
    
  }else if(sheet_name == "Table_S7b"){
    cat("Creating mutated genes for each clone of donors in",sheet_name,"\n")
    donor_clone <- data_sheet$Donor_and_iPSC_clone
    unique_donor_clone <- unique(donor_clone)
    for(donor_clone in unique_donor_clone){
      if(!is.na(donor_clone)){
        rows_index = data_sheet$Donor_and_iPSC_clone == donor_clone
        mutated_genes <- data_sheet[rows_index,]$Gene_symbol
        donor_clone_mutated_genes[[donor_clone]] <- mutated_genes[!is.na(mutated_genes)]
      }
    }
  }
  cat("Loaded mutated genes per clone of donor:\n")
  print(donor_clone_mutated_genes)
  return(donor_clone_mutated_genes)
}

sheet_name <- "Table_S7a"
mutated_genes_S7a <- get_mutated_genes_donor_clone(sheet_name)

sheet_name <- "Table_S7b"
mutated_genes_S7b <- get_mutated_genes_donor_clone(sheet_name)
