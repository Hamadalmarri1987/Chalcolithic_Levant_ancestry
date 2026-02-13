# This code is for to test which populations can be joined in Iran Neolethic by using qpwave
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

pops_iran_with_outgroups_set_1 <- c(
 
  "Iran_HajjiFiruz_N.AG",         
  
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG',
  
  
  "Georgia_Kotias_Mesolithic.SG", 
  "Georgia_Satsurblia_LateUP.SG", 
  "Morocco_Iberomaurusian.AG", 
  "Russia_AfontovaGora3_UP.AG", 
  "Serbia_IronGates_Mesolithic.AG",
  "Mbuti.DG"
)

pops_iran_with_outgroups_set_2 = c(
  
  "Iran_HajjiFiruz_N.AG",          
  
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG',
  
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
f2_blocks_iran_outgroups_set_1 <- f2_from_geno(prefix,
                                                   pops = pops_iran_with_outgroups_set_1,
                                                   maxmiss = 0.1,
                                                   minmaf = 0.01, 
                                                   auto_only = TRUE
)

# f2 with outgroups set 2
f2_blocks_iran_outgroups_set_2 <- f2_from_geno(prefix,
                                               pops = pops_iran_with_outgroups_set_2,
                                               maxmiss = 0.1,
                                               minmaf = 0.01, 
                                               auto_only = TRUE
)


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

left_iran = c(
  "Iran_HajjiFiruz_N.AG",         
  
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG'
)

# Performing QpWave with all pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with all pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop



left_iran = c(
  "Iran_HajjiFiruz_N.AG",         
  
  "Iran_Luristan_PPN.SG"
)

# Performing QpWave with HajjiFiruz + Luristan pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with HajjiFiruz + Luristan pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_HajjiFiruz_N.AG",         
  
  "Iran_GanjDareh_N.AG"
)

# Performing QpWave with HajjiFiruz + GanjDareh pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with HajjiFiruz + GanjDareh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_HajjiFiruz_N.AG",         
  
  "Iran_TepeAbdulHosein_N.SG"
)

# Performing QpWave with HajjiFiruz + TepeAbdulHosein pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with HajjiFiruz + TepeAbdulHosein pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_HajjiFiruz_N.AG",         
  
  'Iran_Wezmeh_N.SG'
)

# Performing QpWave with HajjiFiruz + Wezmeh  pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with HajjiFiruz + Wezmeh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG"
)

# Performing QpWave with Luristan + GanjDareh pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with Luristan + GanjDareh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_Luristan_PPN.SG",
  
  "Iran_TepeAbdulHosein_N.SG"
)

# Performing QpWave with Luristan + TepeAbdulhosein pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with Luristan + TepeAbdulHosein pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_Luristan_PPN.SG",
  
  'Iran_Wezmeh_N.SG'
)

# Performing QpWave with Luristan + Wezmeh pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with Luristan + Wezmeh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG"
)

# Performing QpWave with GanjDareh + TepeAbdulHosein pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with GanjDareh + TepeAbdulHosein pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop



left_iran = c(
  "Iran_GanjDareh_N.AG",
  
  'Iran_Wezmeh_N.SG'
)

# Performing QpWave with GankDareh + Wezmeh pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with GanjDareh + Wezmeh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


left_iran = c(
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG'
)

# Performing QpWave with TepeAbdulHosein + Wezmeh pops with outgroup set 1


qp_iran_1 = qpwave(f2_blocks_iran_outgroups_set_1, left=left_iran, right=right_set_1)
qp_iran_1$rankdrop


# Performing QpWave with TepeAbdulHosein + Wezmeh pops with outgroup set 2


qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop

