library(discover)

### ###
# Run discover tool for a given variant excel file (with cna and allele freq per variant gene)
# Publication: Canisius, Sander, John WM Martens, and Lodewyk FA Wessels. "A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence."
# Genome biology 17.1 (2016): 261.
# Credits: https://ccb.nki.nl/software/discover/doc/r/discover-intro.html
### ###

# Set up directory, input sanple file, list of enriched genes for sample
dir_path = "/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_cna_mutAllFreq"
setwd(dir_path)
data_file <-"hSVEC38C5.csv" #"hSVEC31C4.csv" , "hSVEC31C4.csv", "discover_example.csv" #try with the smallest variant cna allele freq csv:
discover_data_file <- gsub(".csv","_4discover.csv",data_file)
# enriched_genes <- c("hivep3","armh1", "ptpn23", "csmd3", "sall1","ctdsp2","mamstr") #hsvec31c4
enriched_genes <- c("hpcal1","hnrnpf","ppp1cc","rexo1","fbxw9","jak1","slc25a44","cep350","erc2","cldn4","kpna7","tmem168","naif1","myo7a","scel","hirip3","pla2g3","xiap")

cat("\n###\n")
cat("Converting variant csv with cna and allele freq to cna presence for discover tool\n")
cat("###\n")
all_data <- read.csv(data_file)
all_data$studies<- NULL #remove studies column
variant_data <- all_data[, which(colMeans(!is.na(all_data)) > 0.25)]

# get only the CNA data and not the allele freq columns
last_cna_index = dim(variant_data)[2] / 2
cna_data <- variant_data[1:last_cna_index]

for(r in 1:dim(cna_data)[1]){
  cna_row <- c()
  for(c in 2:dim(cna_data)[2]){
      if(is.nan(cna_data[r,c])){ #if there is change in cna, save gene name for association rules
        cna_data[r,c] <- 0
      }
      else if (is.na(cna_data[r,c])){
        cna_data[r,c] <- 0
      }
      else if (cna_data[r,c] < 0 || cna_data[r,c] > 0){
        cna_data[r,c] <- 1
      }
  }
}

# write updated cna data to csv
cat("Saving processed CNA data for input to DISCOVER tool\n")
processed_cna <- as.matrix(t(cna_data))
write.table(processed_cna, file = discover_data_file, sep = ",", col.names = FALSE, row.names = TRUE, append=FALSE)

### ###
# Run DISCOVER
### ###
print("=== ===")
print("Read 0/1 mutations from csv file")
cna4discover <- read.csv(discover_data_file,header=TRUE, row.names=1)

# create background matrix
print("=====")
print("Create the DISCOVER background matrix")
events <- discover.matrix(cna4discover)

# Pairwise test
print("=====")
print("Apply DISCOVER for pairwise tests on the subset genes")
# set the threshold for pairwise mutation test to all genes that are mutated at least on 10% of the samples
cancer_mutation_threshold <- 0.1 * dim(cna4discover)[2]
subset_genes <- rowSums(cna4discover) >= cancer_mutation_threshold
result.mutex <- pairwise.discover.test(events[subset_genes, ])
result.mutex

print("=====")
print("Results with FDR=0.01: ")
print(result.mutex, fdr.threshold=0.01)

print("=====")
pairwise_file = gsub(".csv","_discover_pairwise_ascending_pvalue.csv",data_file)
cat("Save pairwise results in ", pairwise_file, "\n")
pairwise_result <- as.data.frame(result.mutex)
write.table(pairwise_result, file = pairwise_file, sep = ",", col.names = TRUE, row.names = FALSE, append=FALSE)

# Groupwise test
print("=====")
enriched_matched <- sapply(row.names(cna4discover), function(l) length(intersect(enriched_genes, l)) > 0)
groupwise.discover.test(events[enriched_matched, ])
#groupwise.discover.test(events[row.names(cna4discover), ])
print("=====")
plot(events[enriched_matched, ])
plot(events[row.names(cna4discover), ])
print("== * ==")
print("==***==")
