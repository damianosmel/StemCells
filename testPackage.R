###
# Check if the mutations found on the site are also retrieved via the API
###

#in the site we get three mutations for 
#https://www.cbioportal.org/results/mutations?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=msk_impact_2017&case_set_id=msk_impact_2017_cnaseq&clinicallist=NUM_SAMPLES_PER_PATIENT&data_priority=0&gene_list=JAK1&geneset_list=%20&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=msk_impact_2017_cna&genetic_profile_ids_PROFILE_MUTATION_EXTENDED=msk_impact_2017_mutations&show_samples=false&tab_index=tab_visualize

#we can find the N187S mutation also by:
library(cgdsr)
#get connection object and test its functionality
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

profile = getProfileData(mycgds,"JAK1",c("msk_impact_2017_cna","msk_impact_2017_mutations"),"msk_impact_2017_Non-Small_Cell_Lung_Cancer")

#get rows that they have the N187S mutation:
subset(profile, msk_impact_2017_mutations == "N187S")

