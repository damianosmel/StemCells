library(cgdsr)
#get connection object and test its functionality
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

#get all cancer studies
#mycancerstudy = getCancerStudies(mycgds)
#get available genetic measurements for these cancer studies
mycancerstudy = "msk_impact_2017"
#mycancerstudy = "nsclc_tcga_broad_2016"
cat("Print available genetic profiles for ", mycancerstudy, ":")
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)
mygeneticprofile[,c(1:3)]
cat("---\n")

cat("Case list for ", mycancerstudy,":")
mycaselist = getCaseLists(mycgds,mycancerstudy)
allcancertypes = mycaselist[,"case_list_id"][22:23]
cat("---\n")

cat("For one gene, we can get multiple genetic profiles: \n")
cat("For one gene, JAK1, get copy number alternation and mutation information: \n")
for(cancertype in allcancertypes){
  cat("\nType: ", cancertype)
  genomicprofile = getProfileData(mycgds,"JAK1",c("msk_impact_2017_cna","msk_impact_2017_mutations"),cancertype)
  mutationprofile = getMutationData(mycgds,cancertype,"msk_impact_2017_mutations","JAK1")
  cat("\nCNA and Mutations: ")
  print(genomicprofile)
  cat("\nMutations: ")
  print(mutationprofile)
  cat("\n")
  # allele freq = variant read count / (reference read count + variant read count)
  for (mutationrow in 1:nrow(mutationprofile)){
    mutation <- mutationprofile[mutationrow,]
    allelefreq = mutation[1,"variant_read_count_tumor"] / (mutation[1,"reference_read_count_tumor"] + mutation[1,"variant_read_count_tumor"])
    #allelefreqtype = "amplification"
    # if(allelefreq >= 1.0){
    #   allelefreqtype = "amplification"
    # }
    # else{
    #   allelefreqtype = "deletion"
    # }
    #cat("Mutation",mutationrow,"allele frequency=",allelefreq,"(",allelefreqtype,")\n")
    cat("Mutation",mutationrow,"allele frequency=",allelefreq,"\n")
  }
  cat("======\n")
  break
}
cat("---\n")

# cat("For many genes we can select only one genetic profiles: \n")
# cat("For two genes JAK1, XIAP, get copy number alternation: \n")
# for(cancertype in allcancertypes){
#   cat("\nType: ", cancertype)
#   #genomicprofile = getProfileData(mycgds,c("JAK1","XIAP"), "msk_impact_2017_cna",cancertype)
#   #genomicprofile = getProfileData(mycgds,c("JAK1","XIAP"), "msk_impact_2017_mutations",cancertype)
#   #mutationprofile = getMutationData(mycgds,cancertype,"msk_impact_2017_mutations",c("JAK1","XIAP"))
#   mutationprofile = getMutationData(mycgds,cancertype,"msk_impact_2017_mutations",c("XIAP"))
#   #cat("\nCNA: ")
#   #print(genomicprofile)
#   cat("\nMUtations: ")
#   print(mutationprofile)
#   break
# }
# cat("---\n")