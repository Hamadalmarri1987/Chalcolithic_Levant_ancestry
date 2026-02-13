# This code is for to test which populations can be joined in levant by using qpwave
# R version of Admixtools is used
# Version of admixtools ‘2.0.10’
# R version 4.3.0 (2023-04-21 ucrt)

#Importing library
library(admixtools)

#Set working directory and paths
aadr_path <- "C:/Users/WAQAS/Documents/CG/1/22/AA"
prefix <- file.path(aadr_path, "v62.0_1240k_public")

#Preparing for f2
#There will be two f2 with one set of outgroup each

pops_levant_with_outgroups_set_1 <- c(
  "Israel_PPNB.AG",
  "Jordan_Baja_Late_PPNB.AG",
  "Jordan_PPNB.AG", 
  
  "Georgia_Kotias_Mesolithic.SG", 
  "Georgia_Satsurblia_LateUP.SG", 
  "Morocco_Iberomaurusian.AG", 
  "Russia_AfontovaGora3_UP.AG", 
  "Serbia_IronGates_Mesolithic.AG",
  "Mbuti.DG"
)

pops_levant_with_outgroups_set_2 = c(
  "Israel_PPNB.AG",
  "Jordan_Baja_Late_PPNB.AG",
  "Jordan_PPNB.AG",
  
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
f2_blocks_levant_outgroups_set_1 <- f2_from_geno(prefix,
                                                 pops = pops_levant_with_outgroups_set_1,
                                                 maxmiss = 0.1,
                                                 minmaf = 0.01, 
                                                 auto_only = TRUE)


# f2 with outgroups set 2
f2_blocks_levant_outgroups_set_2 <- f2_from_geno(prefix,
                                                 pops = pops_levant_with_outgroups_set_2,
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

left_levant = c(
  "Israel_PPNB.AG",
  "Jordan_Baja_Late_PPNB.AG",
  "Jordan_PPNB.AG" 
)

# Performing QpWave with all pops with outgroup set 1


qp_levant_1 = qpwave(f2_blocks_levant_outgroups_set_1, left=left_levant, right=right_set_1)
qp_levant_1$rankdrop


# Performing QpWave with all pops with outgroup set 2


qp_levant_2 = qpwave(f2_blocks_levant_outgroups_set_2, left=left_levant, right=right_set_2)
qp_levant_2$rankdrop

# QpWave with Israel_PPNB.AG + Jordan_Baja_Late_PPNB.AG left pops

left_levant = c(
  "Israel_PPNB.AG",
  "Jordan_Baja_Late_PPNB.AG"
)

# Performing QpWave with Israel_PPNB.AG + Jordan_Baja_Late_PPNB.AG pops with outgroup set 1


qp_levant_1 = qpwave(f2_blocks_levant_outgroups_set_1, left=left_levant, right=right_set_1)
qp_levant_1$rankdrop


# Performing QpWave with Israel_PPNB.AG + Jordan_Baja_Late_PPNB.AG pops with outgroup set 2


qp_levant_2 = qpwave(f2_blocks_levant_outgroups_set_2, left=left_levant, right=right_set_2)
qp_levant_2$rankdrop


# QpWave with Israel_PPNB.AG + Jordan_PPNB.AG left pops

left_levant = c(
  "Israel_PPNB.AG",
  "Jordan_PPNB.AG"
)

# Performing QpWave with Israel_PPNB.AG + Jordan_PPNB.AG pops with outgroup set 1


qp_levant_1 = qpwave(f2_blocks_levant_outgroups_set_1, left=left_levant, right=right_set_1)
qp_levant_1$rankdrop


# Performing QpWave with Israel_PPNB.AG + Jordan_PPNB.AG pops with outgroup set 2


qp_levant_2 = qpwave(f2_blocks_levant_outgroups_set_2, left=left_levant, right=right_set_2)
qp_levant_2$rankdrop


# QpWave with Jordan_Baja_Late_PPNB.AG + Jordan_PPNB.AG left pops

left_levant = c(
  "Jordan_Baja_Late_PPNB.AG",
  "Jordan_PPNB.AG"
)

# Performing QpWave with Jordan_Baja_Late_PPNB.AG + Jordan_PPNB.AG pops with outgroup set 1


qp_levant_1 = qpwave(f2_blocks_levant_outgroups_set_1, left=left_levant, right=right_set_1)
qp_levant_1$rankdrop


# Performing QpWave with Jordan_Baja_Late_PPNB.AG + Jordan_PPNB.AG pops with outgroup set 2


qp_levant_2 = qpwave(f2_blocks_levant_outgroups_set_2, left=left_levant, right=right_set_2)
qp_levant_2$rankdrop

