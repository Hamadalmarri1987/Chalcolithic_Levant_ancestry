# This code is for to test which populations can be joined in iran chalcolithic by using qpwave
# R version of Admixtools is used
# Version of admixtools ‘2.0.10’
# R version 4.3.0 (2023-04-21 ucrt)

#Importing Library
library(admixtools)

#Set working directory and paths
aadr_path <- "path/to/aadr/folder"
prefix <- file.path(aadr_path, "v62.0_1240k_public")



pops_iran_chalcolithic_with_outgroup_set_1 <- c(
  
  "Iran_SehGabi_C.AG",
  "Iran_TepeHissar_C.DG",
  
  
  "Georgia_Kotias_Mesolithic.SG", 
  "Georgia_Satsurblia_LateUP.SG", 
  "Morocco_Iberomaurusian.AG", 
  "Russia_AfontovaGora3_UP.AG", 
  "Serbia_IronGates_Mesolithic.AG",
  "Mbuti.DG"
)
  
pops_iran_chalcolithic_with_outgroup_set_2 = c(
  "Iran_SehGabi_C.AG",
  "Iran_TepeHissar_C.DG",
  
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",                  
  "Han.DG",               
  "Papuan.DG",            
  "Ethiopia_4500BP.SG",              
  "Karitiana.DG",            
  "Mbuti.DG"
)


# f2 with outgroups set 1
f2_blocks_iran_chalcolithic_outgroups_set_1 <- f2_from_geno(prefix,
                          pops = pops_iran_chalcolithic_with_outgroup_set_1,
                          maxmiss = 0.1,
                          minmaf = 0.01, 
                          auto_only = TRUE)

# f2 with outgroups set 2
f2_blocks_iran_chalcolithic_outgroups_set_2 <- f2_from_geno(prefix,
                          pops = pops_iran_chalcolithic_with_outgroup_set_2,
                          maxmiss = 0.1,
                          minmaf = 0.01, 
                          auto_only = TRUE)

# QpWave 

right_set_1 = c(
  "Georgia_Kotias_Mesolithic.SG", 
  "Georgia_Satsurblia_LateUP.SG", 
  "Morocco_Iberomaurusian.AG", 
  "Russia_AfontovaGora3_UP.AG", 
  "Serbia_IronGates_Mesolithic.AG",
  "Mbuti.DG"
)

right_set_2 =c(
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",                  
  "Han.DG",               
  "Papuan.DG",            
  "Ethiopia_4500BP.SG",              
  "Karitiana.DG",            
  "Mbuti.DG"
)

# QpWave with all the left pops

left_iran_chalcolithic = c(
  
  "Iran_SehGabi_C.AG",
  "Iran_TepeHissar_C.DG"
  
)

# Performing QpWave with all pops with outgroup set 1


qp_levant_1 = qpwave(f2_blocks_iran_chalcolithic_outgroups_set_1, left=left_iran_chalcolithic, right=right_set_1)
qp_levant_1$rankdrop


# Performing QpWave with all pops with outgroup set 2


qp_levant_2 = qpwave(f2_blocks_iran_chalcolithic_outgroups_set_2, left=left_iran_chalcolithic, right=right_set_2)
qp_levant_2$rankdrop



