library(cgdsr)
### 
# Credits
# https://stackoverflow.com/questions/19324003/r-split-string-by-tab-or-nonprinting-character
# https://www.biostars.org/p/163635/
###
#get connection object and test its functionality
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

mycancerstudy = "laml_tcga_pan_can_atlas_2018"

cat("Print available genetic profiles for ", mycancerstudy, ":")
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)
mygeneticprofile[,c(1:3)]
cat("---\n")

cat("Case list for ", mycancerstudy,":")
mycaselist = getCaseLists(mycgds,mycancerstudy)
#print(mycaselist)
cases = mycaselist[7,"case_ids"]
cases = strsplit(cases,' ')[[1]]
print(cases)

secondcancerstudy = "laml_tcga"

cat("Print available genetic profiles for ", secondcancerstudy,":")
secondstudyprofile = getGeneticProfiles(mycgds,secondcancerstudy)
secondstudyprofile[,c(1:3)]
cat("---\n")

cat("Case list for ",secondcancerstudy,":")
secondcaselist = getCaseLists(mycgds,secondcancerstudy)
#print(secondcaselist)
secondcases = secondcaselist[10,"case_ids"]
secondcases = strsplit(secondcases,' ')[[1]]
print(secondcases)

commoncases = intersect(cases, secondcases)
print(commoncases)
cat(length(commoncases))
