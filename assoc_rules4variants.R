library(arules)
library(arulesViz)

### ###
# Run association rules for a given variant excel file (with cna and allele freq per variant gene)
### ###
#credits: https://michael.hahsler.net/research/arules_RUG_2015/demo/

dir_path = "/home/damian/Documents/L3S/projects/stemCells/maike/dataset/all_variants_cna_mutAllFreq"
data_file <-"hSVEC38C6.csv" #"hSVEC38C5.csv", "hSVEC31C4.csv" #try with the smallest variant cna allele freq csv:
arules_data_file <- gsub(".csv","_4arules.csv",data_file)
setwd(dir_path)

### ####
cat("\n###\n")
cat("Converting variant csv with cna and allele freq to cna presence for arules.\n")
cat("###\n")
all_data <- read.csv(data_file)
data_dim <- dim(all_data)
variant_data <- all_data[,3:data_dim[2]] #exclude samples and studies

last_cna_index = dim(variant_data)[2] / 2
cna_data <- variant_data[3:last_cna_index]
#append vector, credits: https://stackoverflow.com/questions/8617347/append-a-vector-as-a-row-in-a-csv-file
for(r in 1:dim(cna_data)[1]){
  cna_row <- c()
  for(c in 1:dim(cna_data)[2]){
    if (is.finite(cna_data[r,c]) == TRUE){

      if(cna_data[r,c] == 0){ #if there is change in cna, save gene name for association rules
        cna_row <- append(cna_row,names(cna_data)[c]) 
      }
    }
  }
  if (length(cna_row)>= 1){
    cat("Converting a new study instance\n")
    cna <- as.matrix(t(cna_row))
    write.table(cna, file = arules_data_file, sep = ",", 
                col.names = FALSE, row.names = FALSE, append=TRUE)
  }
}

### ###
cat("\n###\n")
cat("Converting to transactions")
cat("Getting first summary of data\n")
cat("###\n")

cna_transactions <- read.transactions(arules_data_file, format="basket",sep = ",")
summary(cna_transactions)

cat("Number of transactions:")
cat(nrow(cna_transactions))
itemFrequencyPlot(cna_transactions, topN=dim(cna_transactions)[2], cex.names=.5,horiz=TRUE)

min_support <- 0.1 #2/5
#max_len 3rd quartile
max_len <- 3.0#, 162 for hSVEC38C5 crushes the rstudio so I have used 10 didn't work so I used 6 in the end, #79.0, for hSVEC31C4.csv
cna_itemsets <- apriori(cna_transactions, parameter = list(target = "frequent",supp=min_support, minlen = 2, maxlen=max_len))
inspect(head(sort(cna_itemsets),n=50))

#sort by interest
quality(cna_itemsets)$lift <- interestMeasure(cna_itemsets, measure="lift", trans = cna_transactions)
inspect(head(sort(cna_itemsets, by = "lift"), n=50))

#plot item sets in graph
max_num_freq_itemsets = length(cna_itemsets)
plot(head(sort(cna_itemsets, by = "lift"), n=max_num_freq_itemsets), method="graph",control=list(layout=igraph::with_graphopt(spring.const=5, mass=200),cex=0.4,arrowSize=0.25))

#save found cna itemsets
arules_save_file <- gsub(".csv","_arules.csv",data_file)
write(sort(cna_itemsets,by="lift"), file = arules_save_file, sep = ",", quote = TRUE, row.names = FALSE)
all_mutated_genes <- read.csv(arules_save_file)

### ###
# Function to get the presence/absence of interesting genes in a row of all mutated genes data-frame
### ###
find_interesting_genes_in_mutated_genes = function(x, output) {
  # access mutated genes  in first column
  itemsets = as.character(x[1])
  genes = unlist(strsplit(substr(itemsets,2,nchar(itemsets)-1), ','),use.names=FALSE)
  maike_genes <- c("hivep3","armh1","ptpn23","csmd3","sall1")
  presence_interest_genes <- maike_genes %in% genes
  return(is.element(TRUE,presence_interest_genes))
}

# write rules for Maike genes in file for manual check
arules_subset_genes_file <- gsub(".csv","_arules_maike_genes.csv",data_file)
subset_mutated_genes <- subset(all_mutated_genes, apply(all_mutated_genes,1,find_interesting_genes_in_mutated_genes), c("items","support","count","lift"), drop = FALSE)
write.csv(subset_mutated_genes, file=arules_subset_genes_file, row.names=FALSE)
