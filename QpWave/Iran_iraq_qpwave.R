# This code is for to test which populations can be joined in Iran Neolethic by using qpwave
# R version of Admixtools is used
# Version of admixtools ‘2.0.10’
# R version 4.3.0 (2023-04-21 ucrt)

#Importing library
library(admixtools)

#Set working directory and paths
aadr_path <- "path/to/aadr/folder"
prefix <- file.path(aadr_path, "v62.0_1240k_public")

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

#Preparing for f2
#There will be f2 with only outgroup set 2


#f2 with only iran
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


# f2 with outgroups set 2
f2_blocks_iran_outgroups_set_2 <- f2_from_geno(prefix,
                                               pops = pops_iran_with_outgroups_set_2,
                                               maxmiss = 0.1,
                                               minmaf = 0.01, 
                                               auto_only = TRUE
)

# Performing QpWave with all pops with outgroup set 2

left_iran = c(
  "Iran_HajjiFiruz_N.AG",          
  
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG'
)

qp_iran_2 = qpwave(f2_blocks_iran_outgroups_set_2, left=left_iran, right=right_set_2)
qp_iran_2$rankdrop


# f2 with only iraq

pops_iraq_with_outgroups_set_2 = c(
  
  "Iraq_Shanidar.AG", 
  "Iraq_PPNA.AG",
  
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",                  
  "Han.DG",               
  "Papuan.DG",            
  "Ethiopia_4500BP.SG",              
  "Karitiana.DG",            
  "Mbuti.DG"
)


# f2 with outgroups set 2
f2_blocks_iraq_outgroups_set_2 <- f2_from_geno(prefix,
                                               pops = pops_iraq_with_outgroups_set_2,
                                               maxmiss = 0.1,
                                               minmaf = 0.01, 
                                               auto_only = TRUE
)

# Performing QpWave with all pops with outgroup set 2

left_iraq = c(
  "Iraq_Shanidar.AG", 
  "Iraq_PPNA.AG"
)

qp_iraq_2 = qpwave(f2_blocks_iraq_outgroups_set_2, left=left_iraq, right=right_set_2)
qp_iraq_2$rankdrop



# F2 with all pops

pops_iraq_iran_with_outgroups_set_2 = c(
  
  "Iraq_Shanidar.AG", 
  "Iraq_PPNA.AG",
  
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


# f2 with outgroups set 2
f2_blocks_iraq_iran_outgroups_set_2 <- f2_from_geno(prefix,
                                               pops = pops_iraq_iran_with_outgroups_set_2,
                                               maxmiss = 0.1,
                                               minmaf = 0.01, 
                                               auto_only = TRUE
)

left_iraq_iran = c(
  "Iraq_Shanidar.AG", 
  "Iraq_PPNA.AG",
  
  "Iran_HajjiFiruz_N.AG",          
  
  "Iran_Luristan_PPN.SG",
  
  "Iran_GanjDareh_N.AG",
  
  "Iran_TepeAbdulHosein_N.SG",  
  
  'Iran_Wezmeh_N.SG'
)

qp_iraq_iran_2 = qpwave(f2_blocks_iraq_outgroups_set_2, left=left_iraq_iran, right=right_set_2)
qp_iraq_iran_2$rankdrop



outgroups <- c(
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",
  "Han.DG",               
  "Papuan.DG",
  "Ethiopia_4500BP.SG",
  "Karitiana.DG",
  "Mbuti.DG"
)

iran_pops <- c(
  "Iran_HajjiFiruz_N.AG",
  "Iran_Luristan_PPN.SG",
  "Iran_GanjDareh_N.AG",
  "Iran_TepeAbdulHosein_N.SG",
  "Iran_Wezmeh_N.SG"
)

iraq_pops <- c(
  "Iraq_Shanidar.AG",
  "Iraq_PPNA.AG"
)

all_pops <- c(iran_pops, iraq_pops)

# GENERATE VALID COMBINATIONS (only sizes 2,3,4)
valid_combos <- function(all_pops, sizes) {
  combos <- list()
  idx <- 1
  
  for(k in sizes){
    cmb <- combn(all_pops, k, simplify = FALSE)
    
    for(c in cmb){
      # keep only mixed Iran-Iraq sets
      if(any(c %in% iran_pops) && any(c %in% iraq_pops)){
        combos[[idx]] <- c
        idx <- idx + 1
      }
    }
  }
  return(combos)
}

combo_list <- valid_combos(all_pops, sizes = c(2,3,4))

length(combo_list)   # number of valid sets
combo_list            # inspect sets


# RUN f2 + qpWave FOR EACH COMBINATION 

results <- list()

for(i in seq_along(combo_list)){
  
  left_set <- combo_list[[i]]
  name <- paste(left_set, collapse="__")
  
  cat("\n=========================\n")
  cat("Running combination:", paste(left_set, collapse=", "), "\n")
  cat("=========================\n")
  
  pops_for_f2 <- c(left_set, outgroups)
  
  # ---- Capture the printed output from f2_from_geno ----
  f2_output <- capture.output(
    f2_blocks <- f2_from_geno(prefix,
                              pops = pops_for_f2,
                              maxmiss = 0.1,
                              minmaf = 0.01,
                              auto_only = TRUE)
  )
  
  # FIX: use f2_output instead of out 
  snps_total_line    <- grep("SNPs read in total", f2_output, value = TRUE)
  snps_filtered_line <- grep("remain after filtering", f2_output, value = TRUE)
  
  
  snps_total    <- if (length(snps_total_line) > 0) snps_total_line else NA
  snps_filtered <- if (length(snps_filtered_line) > 0) snps_filtered_line else NA
  
  
  # Print for user
  cat("✔", snps_total, "\n")
  cat("!", snps_filtered, "\n")
  
  
  # qpWave
  qp <- qpwave(f2_blocks, left = left_set, right = outgroups)
  
  # Save output 
  results[[name]] <- list(
    snps_total    = snps_total,
    snps_filtered = snps_filtered,
    
    rankdrop      = qp$rankdrop
  )
  
  print(qp$rankdrop)
}


results

