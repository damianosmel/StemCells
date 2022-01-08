### ###
# Function to read Maike's vcf files, formatted in .tabular and create a list of genes per donor's clone 
# Having such list cna_mutationAlleleFreq_donor_clone_all.R can run
### ###
#setwd('/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_example')
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
    print(file.names[i])
    #list[donor_clone,var_genes] <- parse_variant_table(file.names[i])
    donor_clone_variants <- parse_variant_table(file.names[i])
    cat("donor_clone: ")
    print(typeof(donor_clone_variants$donor_clone))
    #cat("var_genes: ")
    #print(donor_clone_variants$variants)
    donor_clone2vars[[donor_clone_variants$donor_clone]] <- donor_clone_variants$variants
  }
  return(donor_clone2vars)
}

dir_path = "/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_example"
donor_clone2vars <- parse_all_variants(dir_path)