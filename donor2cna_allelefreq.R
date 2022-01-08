library(cgdsr)
#get connection object and test its functionality
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

#studies to include
studies_list = ["msk_impact_2017","tcga_",..]

#get the list of genes for one donor
genes_list = ["JAK1", "XIAP"]

cat("Print available genetic profiles for ", mycancerstudy, ":")
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)
mygeneticprofile[,c(1:3)]
cat("---\n")


#get all samples - update
cat("Case list for ", mycancerstudy,":")
mycaselist = getCaseLists(mycgds,mycancerstudy)
allcancertypes = mycaselist[,"case_list_id"][22:23]
cat("---\n")

#get cna and mutation allele freq
#save sampe_id \tab cna;alleleFreq 
cat("For many genes we can select only one genetic profiles: \n")
cat("For two genes JAK1, XIAP, get copy number alternation: \n")
for(cancertype in allcancertypes){
  cat("\nType: ", cancertype)
  #genomicprofile = getProfileData(mycgds,c("JAK1","XIAP"), "msk_impact_2017_cna",cancertype)
  genomicprofile = getProfileData(mycgds,c("JAK1","XIAP"), "msk_impact_2017_mutations",cancertype)
  mutationprofile = getMutationData(mycgds,cancertype,"msk_impact_2017_mutations",c("JAK1","XIAP"))
  cat("\nCNA: ")
  print(genomicprofile)
  cat("\nMUtations: ")
  print(mutationprofile)
}
cat("---\n")