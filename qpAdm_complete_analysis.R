#Comprehensive qpAdm Analysis for Levant Chalcolithic Ancestry
#Testing 3-way and 4-way models with multiple outgroup sets
#Including merged populations based on qpWave results

library(admixtools)
library(tidyverse)
library(readxl)
library(writexl)
library(magrittr)

#Set working directory and paths
results_dir <- "results_aDNA/"
dir.create(results_dir, showWarnings = FALSE)

aadr_path <- "AADR_v62"
prefix <- file.path(aadr_path, "v62.0_1240k_public")

#Define target population: Levant Chalcolithic
target_population <- "Israel_C.AG"


#POPULATION DEFINITIONS


levant_pops <- c(
  "Jordan_PPNB.AG",
  "Israel_PPNB.AG", 
  "Jordan_Baja_Late_PPNB.AG"
)

anatolia_pops <- c(
  "Turkey_Marmara_Barcin_N.SG",
  "Turkey_Central_Catalhoyuk_N.SG",
  "Turkey_Southeast_Cayonu_PPN.SG",
  "Turkey_Central_Boncuklu_PPN.AG",
  "Turkey_Central_Boncuklu_PPN.WGC.SG"
)


iran_pops <- c(
  "Iran_HajjiFiruz_N.AG",
  "Iran_GanjDareh_N.AG",
  "Iran_TepeAbdulHosein_N.SG",
  "Iran_Wezmeh_N.SG",
  "Iran_Luristan_PPN.SG"
)

iran_chl_pops <- c(
  "Iran_SehGabi_C.AG"
)

#Core Iraq populations
iraq_pops <- c(
  "Iraq_Shanidar.AG",
  "Iraq_PPNA.AG"
)

#Southeast Anatolia populations = Upper Mesopotamia per Lazaridis
se_anatolia_pops <- c(
  "Turkey_Southeast_Cayonu_PPN.SG",
  "Turkey_Southeast_Mardin_PPN.AG"
)

#FULL MESOPOTAMIA per Lazaridis = SE Anatolia + Iraq
mesopotamia_full <- c(
  "Turkey_Southeast_Cayonu_PPN.SG",
  "Turkey_Southeast_Mardin_PPN.AG",
  "Iraq_Shanidar.AG",
  "Iraq_PPNA.AG"
)


#MERGED POPULATION DEFINITIONS (based on qpWave results)

#Levant merged (all three form a clade, p > 0.5 for all combinations)
levant_merged <- c("Jordan_PPNB.AG", "Israel_PPNB.AG", "Jordan_Baja_Late_PPNB.AG")

#Anatolia merges based on qpWave results
#Barcin + Catalhoyuk (p = 0.723)
anatolia_barcin_catalhoyuk <- c("Turkey_Marmara_Barcin_N.SG", "Turkey_Central_Catalhoyuk_N.SG")

#Barcin + Cayonu (p = 0.334)
anatolia_barcin_cayonu <- c("Turkey_Marmara_Barcin_N.SG", "Turkey_Southeast_Cayonu_PPN.SG")

#Catalhoyuk + Cayonu (p = 0.828)
anatolia_catalhoyuk_cayonu <- c("Turkey_Central_Catalhoyuk_N.SG", "Turkey_Southeast_Cayonu_PPN.SG")

#Catalhoyuk + Boncuklu WGC (p = 0.911)
anatolia_catalhoyuk_boncuklu <- c("Turkey_Central_Catalhoyuk_N.SG", "Turkey_Central_Boncuklu_PPN.WGC.SG")

#Barcin + Catalhoyuk + Cayonu (p = 0.537 for f4rank=0)
anatolia_barcin_catalhoyuk_cayonu <- c("Turkey_Marmara_Barcin_N.SG", "Turkey_Central_Catalhoyuk_N.SG", 
                                        "Turkey_Southeast_Cayonu_PPN.SG")

#Iran merges based on qpWave results
#GanjDareh + TepeAbdulHosein (p = 0.153 / 0.527)
iran_ganjdareh_tepe <- c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG")

#TepeAbdulHosein + Wezmeh (p = 0.506)
iran_tepe_wezmeh <- c("Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG")

#GanjDareh + Wezmeh (p = 0.185)
iran_ganjdareh_wezmeh <- c("Iran_GanjDareh_N.AG", "Iran_Wezmeh_N.SG")

#GanjDareh + TepeAbdulHosein + Wezmeh (p = 0.361 for f4rank=0)
iran_zagros_merged <- c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG")

#Mesopotamia merges based on qpWave results
#Iraq pops can be merged (p = 0.554)
iraq_merged <- c("Iraq_Shanidar.AG", "Iraq_PPNA.AG")

#SE Anatolia merged (Cayonu + Mardin, p = 0.937)
se_anatolia_merged <- c("Turkey_Southeast_Cayonu_PPN.SG", "Turkey_Southeast_Mardin_PPN.AG")

#Extended Mesopotamia combinations (Lazaridis approach)
#Iraq_PPNA + Turkey_Southeast_Mardin (p = 0.341)
mesopotamia_ppna_mardin <- c("Iraq_PPNA.AG", "Turkey_Southeast_Mardin_PPN.AG")

#Full Mesopotamia per Lazaridis (all 4 populations)
#May have low SNP overlap
mesopotamia_lazaridis_full <- c("Turkey_Southeast_Cayonu_PPN.SG", "Turkey_Southeast_Mardin_PPN.AG",
                                 "Iraq_Shanidar.AG", "Iraq_PPNA.AG")

#Iran + Iraq combinations that form clades
#HajjiFiruz + Iraq_Shanidar (p = 0.699)
iran_iraq_hajji_shanidar <- c("Iran_HajjiFiruz_N.AG", "Iraq_Shanidar.AG")

#GanjDareh + Iraq_Shanidar (p = 0.524)
iran_iraq_ganj_shanidar <- c("Iran_GanjDareh_N.AG", "Iraq_Shanidar.AG")

#TepeAbdulHosein + Iraq_Shanidar (p = 0.601)
iran_iraq_tepe_shanidar <- c("Iran_TepeAbdulHosein_N.SG", "Iraq_Shanidar.AG")

#Wezmeh + Iraq_Shanidar (p = 0.610)
iran_iraq_wezmeh_shanidar <- c("Iran_Wezmeh_N.SG", "Iraq_Shanidar.AG")

#GanjDareh + TepeAbdulHosein + Iraq_Shanidar (p = 0.868 for f4rank=0)
iran_iraq_zagros_shanidar <- c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG", "Iraq_Shanidar.AG")

#GanjDareh + TepeAbdulHosein + Wezmeh + Iraq_Shanidar (p = 0.617 for f4rank=0)
iran_iraq_full_zagros <- c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG", "Iraq_Shanidar.AG")


#OUTGROUP DEFINITIONS (6 sets)

out1 <- c(
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",                  
  "Han.DG",               
  "Papuan.DG",            
  "Chukchi.DG",              
  "Karitiana.DG",            
  "Mbuti.DG"
)

out2 <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Israel_Natufian.AG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Ethiopia_4500BP.SG",
  "Karitiana.DG",
  "Mbuti.DG"
)

out3 <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Ethiopia_4500BP.SG",
  "Karitiana.DG",
  "Mbuti.DG"
)

out4 <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Chukchi.DG",
  "Karitiana.DG",
  "Mbuti.DG"
)

out5 <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Israel_Natufian.AG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Chukchi.DG",
  "Karitiana.DG",
  "Mbuti.DG"
)

out6 <- c(
  "Russia_UstIshim_IUP.DG",         
  "Russia_Kostenki14_UP.SG",          
  "Russia_MA1_UP.SG",                  
  "Han.DG",               
  "Papuan.DG",            
  "Ethiopia_4500BP.SG",              
  "Karitiana.DG",            
  "Mbuti.DG"
)

#Outgroup set from qpWave document
out_qpwave1 <- c(
  "Georgia_Kotias_Mesolithic.SG",
  "Georgia_Satsurblia_LateUP.SG",
  "Morocco_Iberomaurusian.AG",
  "Russia_AfontovaGora3_UP.AG",
  "Serbia_IronGates_Mesolithic.AG",
  "Mbuti.DG"
)

#Store all outgroup sets in a list
outgroup_sets <- list(
  out1 = out1,
  out2 = out2,
  out3 = out3,
  out4 = out4,
  out5 = out5,
  out6 = out6
)


#Load data

#Read population information
ind_file <- paste0(prefix, ".ind")
anno_file <- paste0(prefix, ".anno")

#Load individual and annotation data
individuals <- read_table(ind_file, col_names = c("ID", "Sex", "Population"))
annotations <- read_tsv(anno_file)


#Helper function to get individuals from a population in the .ind file
get_inds_from_pop <- function(ind_file, pop_name) {
  ind_data <- read.table(ind_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(ind_data) <- c("ind", "sex", "pop")
  ind_data$ind[ind_data$pop == pop_name]
}

#Helper function to read .ind file
read_ind_file <- function(prefix) {
  ind_file <- paste0(prefix, ".ind")
  ind_data <- read.table(ind_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(ind_data) <- c("ind", "sex", "pop")
  return(ind_data)
}

#Function to create inds and pops vectors for population merging
create_merged_pop_mapping <- function(prefix, source_groups, target, right_pops) {
  
  ind_data <- read_ind_file(prefix)
  
  inds_list <- c()
  pops_list <- c()
  
  #Add target population (not merged, keep original name)
  target_inds <- ind_data$ind[ind_data$pop == target]
  inds_list <- c(inds_list, target_inds)
  pops_list <- c(pops_list, rep(target, length(target_inds)))
  
  # Add source populations (with potential merging)
  for (group_name in names(source_groups)) {
    orig_pops <- source_groups[[group_name]]
    for (orig_pop in orig_pops) {
      pop_inds <- ind_data$ind[ind_data$pop == orig_pop]
      if (length(pop_inds) > 0) {
        inds_list <- c(inds_list, pop_inds)
        pops_list <- c(pops_list, rep(group_name, length(pop_inds)))
      }
    }
  }
  
  #Add right (outgroup) populations (not merged, keep original names)
  for (right_pop in right_pops) {
    right_inds <- ind_data$ind[ind_data$pop == right_pop]
    if (length(right_inds) > 0) {
      inds_list <- c(inds_list, right_inds)
      pops_list <- c(pops_list, rep(right_pop, length(right_inds)))
    }
  }
  
  return(list(inds = inds_list, pops = pops_list))
}

run_qpadm_analysis <- function(target, source_groups, right_pops, label, f2_directory = "f2_blocks") {
  
  #Get the merged source names (these will be the qpAdm left populations)
  left_pops_merged <- names(source_groups)
  
  tryCatch({
    cat("\nPreparing merged populations for model:", label, "\n")
    cat("Source groups:\n")
    for (g in names(source_groups)) {
      cat("  ", g, ":", paste(source_groups[[g]], collapse = " + "), "\n")
    }
    
    #Create the individual-to-merged-population mapping
    pop_mapping <- create_merged_pop_mapping(prefix, source_groups, target, right_pops)
    
    cat("Total individuals:", length(pop_mapping$inds), "\n")
    cat("Merged populations:", paste(unique(pop_mapping$pops), collapse = ", "), "\n")
    
    #Capture output from extract_f2 to get SNP count
    output_log <- capture.output({
      extract_f2(
        pref = prefix,
        outdir = f2_directory,
        inds = pop_mapping$inds,
        pops = pop_mapping$pops,
        maxmiss = 0.1,
        minmaf = 0.01,
        auto_only = TRUE,
        overwrite = TRUE
      )
    }, type = "output")
    
    #Print the captured output so user can see progress
    cat(paste(output_log, collapse = "\n"), "\n")
    
    #Extract SNP count from the output
    snp_count <- NA
    snp_line <- grep("SNPs remain after filtering", output_log, value = TRUE)
    if (length(snp_line) > 0) {
      snp_match <- regmatches(snp_line[1], regexpr("[0-9]+", snp_line[1]))
      if (length(snp_match) > 0) {
        snp_count <- as.integer(snp_match)
      }
    }
    
    #Read f2 data back from disk
    #The f2 data will now have the merged population names
    f2_data <- f2_from_precomp(f2_directory)
    
    cat("Testing model:", label, "\n")
    cat("Target:", target, "\n")
    cat("Sources (merged):", paste(left_pops_merged, collapse = " + "), "\n")
    cat("Outgroups:", paste(right_pops[1:min(3, length(right_pops))], collapse = ", "), "...\n")
    cat("SNPs used:", ifelse(is.na(snp_count), "NA", snp_count), "\n")

    #Run qpAdm with merged population names
    result <- qpadm(
      data = f2_data,
      left = left_pops_merged,
      right = right_pops,
      target = target,
      boot = FALSE,
      verbose = TRUE
    )
    
    #Extract weights and p-value
    weights <- result$weights
    p_value <- result$popdrop[1, "p"]
    
    #Check model feasibility
    feasible <- all(weights$weight >= 0 & weights$weight <= 1)
    
    #Create summary data frame
    summary_df <- data.frame(
      Model = label,
      Target = target,
      Sources = paste(left_pops_merged, collapse = " + "),
      N_sources = length(left_pops_merged),
      P_value = p_value,
      Feasible = ifelse(feasible, "TRUE", "FALSE"),
      SNPs_used = snp_count,
      stringsAsFactors = FALSE
    )
    
    #Add ancestry proportions with standard errors for each merged source
    for (i in seq_along(left_pops_merged)) {
      source_name <- left_pops_merged[i]
      summary_df[[paste0(source_name, "_prop")]] <- weights$weight[i]
      summary_df[[paste0(source_name, "_se")]] <- weights$se[i]
    }
    
    #Store result before cleanup
    final_result <- list(summary = summary_df, full_result = result, snp_count = snp_count)
    
    #Clean up f2 directory after successful run to save disk space
    if (dir.exists(f2_directory)) {
      unlink(f2_directory, recursive = TRUE, force = TRUE)
      cat("Cleaned up f2 directory:", f2_directory, "\n")
    }
    
    return(final_result)
    
  }, error = function(e) {
    cat("Error in model", label, ":", conditionMessage(e), "\n")
    
    #Clean up f2 directory even on error
    if (dir.exists(f2_directory)) {
      unlink(f2_directory, recursive = TRUE, force = TRUE)
      cat("Cleaned up f2 directory after error:", f2_directory, "\n")
    }
    
    return(list(
      summary = data.frame(
        Model = label,
        Target = target,
        Sources = paste(names(source_groups), collapse = " + "),
        N_sources = length(source_groups),
        P_value = NA,
        Feasible = NA,
        SNPs_used = NA,
        Error = conditionMessage(e),
        stringsAsFactors = FALSE
      ),
      full_result = NULL,
      snp_count = NA
    ))
  })
}


#Create named lists for each region
#Single populations
single_levant <- list(
  "Levant_JordanPPNB" = "Jordan_PPNB.AG",
  "Levant_IsraelPPNB" = "Israel_PPNB.AG",
  "Levant_JordanBaja" = "Jordan_Baja_Late_PPNB.AG"
)

single_anatolia <- list(
  #NW Anatolia (Marmara region)
  "Anatolia_Barcin" = "Turkey_Marmara_Barcin_N.SG",
  #Central Anatolia
  "Anatolia_Catalhoyuk" = "Turkey_Central_Catalhoyuk_N.SG",
  "Anatolia_Boncuklu" = "Turkey_Central_Boncuklu_PPN.AG",
  "Anatolia_BoncukluWGC" = "Turkey_Central_Boncuklu_PPN.WGC.SG"
)

single_iran <- list(
  "Iran_HajjiFiruz" = "Iran_HajjiFiruz_N.AG",
  "Iran_GanjDareh" = "Iran_GanjDareh_N.AG",
  "Iran_TepeAbdulHosein" = "Iran_TepeAbdulHosein_N.SG",
  "Iran_Wezmeh" = "Iran_Wezmeh_N.SG",
  "Iran_SehGabi_C" = "Iran_SehGabi_C.AG"
)

#Single populations - MESOPOTAMIA (Iraq + SE Anatolia per Lazaridis)
single_mesopotamia <- list(
  #Iraq populations
  "Meso_Iraq_Shanidar" = "Iraq_Shanidar.AG",
  "Meso_Iraq_PPNA" = "Iraq_PPNA.AG",
  #SE Anatolia populations (= Upper Mesopotamia per Lazaridis)
  "Meso_SEAnat_Cayonu" = "Turkey_Southeast_Cayonu_PPN.SG",
  "Meso_SEAnat_Mardin" = "Turkey_Southeast_Mardin_PPN.AG"
)

#Merged populations - MESOPOTAMIA
merged_mesopotamia <- list(
  #Iraq only merged
  "Meso_IraqMerged" = iraq_merged,
  #SE Anatolia only merged
  "Meso_SEAnatoliaMerged" = se_anatolia_merged,
  #Cross combinations
  "Meso_PPNA_Mardin" = mesopotamia_ppna_mardin,
  #Full Lazaridis Mesopotamia
  "Meso_Lazaridis_Full" = mesopotamia_lazaridis_full
)

#Merged populations - LEVANT
merged_levant <- list(
  "Levant_AllMerged" = levant_merged
)

merged_anatolia <- list(
  #NW + Central Anatolia combinations (not including SE which is Mesopotamian)
  "Anatolia_BarcinCatalhoyuk" = anatolia_barcin_catalhoyuk,
  "Anatolia_CatalhoyukBoncuklu" = anatolia_catalhoyuk_boncuklu
)

merged_iran <- list(
  "Iran_GanjDarehTepe" = iran_ganjdareh_tepe,
  "Iran_TepeWezmeh" = iran_tepe_wezmeh,
  "Iran_GanjDarehWezmeh" = iran_ganjdareh_wezmeh,
  "Iran_ZagrosMerged" = iran_zagros_merged
)

#Iran-Iraq combined (Zagros + Mesopotamia clades)
merged_iran_iraq <- list(
  "IranIraq_HajjiShanidar" = iran_iraq_hajji_shanidar,
  "IranIraq_GanjShanidar" = iran_iraq_ganj_shanidar,
  "IranIraq_TepeShanidar" = iran_iraq_tepe_shanidar,
  "IranIraq_WezmehShanidar" = iran_iraq_wezmeh_shanidar,
  "IranIraq_ZagrosShanidar" = iran_iraq_zagros_shanidar,
  "IranIraq_FullZagros" = iran_iraq_full_zagros
)

#GENERATE ALL 3-WAY MODEL COMBINATIONS


generate_3way_models <- function() {
  models <- list()
  
  #Combine single and merged populations
  all_levant <- c(single_levant, merged_levant)
  all_anatolia <- c(single_anatolia, merged_anatolia)
  all_iran <- c(single_iran, merged_iran)
  all_mesopotamia <- c(single_mesopotamia, merged_mesopotamia, merged_iran_iraq)
  
  #Helper function to check for population overlap
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #3-way: Levant + Anatolia + Iran
  for (lev_name in names(all_levant)) {
    for (anat_name in names(all_anatolia)) {
      for (iran_name in names(all_iran)) {
        sources <- c(all_levant[[lev_name]], all_anatolia[[anat_name]], all_iran[[iran_name]])
        if (!has_overlap(sources)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", iran_name)
          models[[model_name]] <- list(
            sources = sources,
            label = model_name
          )
        }
      }
    }
  }
  
  #3-way: Levant + Anatolia + Mesopotamia
  for (lev_name in names(all_levant)) {
    for (anat_name in names(all_anatolia)) {
      for (meso_name in names(all_mesopotamia)) {
        sources <- c(all_levant[[lev_name]], all_anatolia[[anat_name]], all_mesopotamia[[meso_name]])
        if (!has_overlap(sources)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", meso_name)
          models[[model_name]] <- list(
            sources = sources,
            label = model_name
          )
        }
      }
    }
  }
  
  #3-way: Levant + Iran + Mesopotamia
  #Merged_iran_iraq contains both Iran and Iraq pops
  #Skip combinations where Iran source overlaps with Mesopotamia source
  for (lev_name in names(all_levant)) {
    for (iran_name in names(all_iran)) {
      for (meso_name in names(all_mesopotamia)) {
        sources <- c(all_levant[[lev_name]], all_iran[[iran_name]], all_mesopotamia[[meso_name]])
        if (!has_overlap(sources)) {
          model_name <- paste0("3way_", lev_name, "_", iran_name, "_", meso_name)
          models[[model_name]] <- list(
            sources = sources,
            label = model_name
          )
        }
      }
    }
  }
  
  #3-way: Anatolia + Iran + Mesopotamia
  for (anat_name in names(all_anatolia)) {
    for (iran_name in names(all_iran)) {
      for (meso_name in names(all_mesopotamia)) {
        sources <- c(all_anatolia[[anat_name]], all_iran[[iran_name]], all_mesopotamia[[meso_name]])
        if (!has_overlap(sources)) {
          model_name <- paste0("3way_", anat_name, "_", iran_name, "_", meso_name)
          models[[model_name]] <- list(
            sources = sources,
            label = model_name
          )
        }
      }
    }
  }
  
  return(models)
}

#GENERATE SIMPLIFIED 3-WAY MODELS (fewer combinations)

generate_3way_models_simplified <- function() {
  models <- list()
  
  #Define key populations with their GROUP NAMES and population vectors

  key_levant <- list(
    "Levant_Israel" = list(Levant = "Israel_PPNB.AG"),
    "Levant_Merged" = list(Levant = levant_merged)
  )
  
  key_anatolia <- list(
    "Anatolia_Barcin" = list(Anatolia = "Turkey_Marmara_Barcin_N.SG"),
    "Anatolia_Catalhoyuk" = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG")
  )
  
  key_iran <- list(
    "Iran_HajjiFiruz" = list(Iran = "Iran_HajjiFiruz_N.AG"),
    "Iran_ZagrosMerged" = list(Iran = iran_zagros_merged),
    "Iran_SehGabi_C" = list(Iran = "Iran_SehGabi_C.AG")
  )
  
  key_mesopotamia <- list(
    "Meso_Shanidar" = list(Mesopotamia = "Iraq_Shanidar.AG"),
    "Meso_PPNA" = list(Mesopotamia = "Iraq_PPNA.AG"),
    "Meso_IraqMerged" = list(Mesopotamia = iraq_merged),
    "Meso_Cayonu" = list(Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"),
    "Meso_SEAnatoliaMerged" = list(Mesopotamia = se_anatolia_merged),
    "Meso_Lazaridis_Full" = list(Mesopotamia = mesopotamia_lazaridis_full)
  )
  
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #3-way: Levant + Anatolia + Iran
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (iran_name in names(key_iran)) {
        source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], key_iran[[iran_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", iran_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  #3-way: Levant + Anatolia + Mesopotamia
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (meso_name in names(key_mesopotamia)) {
        source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], key_mesopotamia[[meso_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", meso_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  #3-way: Levant + Iran + Mesopotamia
  for (lev_name in names(key_levant)) {
    for (iran_name in names(key_iran)) {
      for (meso_name in names(key_mesopotamia)) {
        source_groups <- c(key_levant[[lev_name]], key_iran[[iran_name]], key_mesopotamia[[meso_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", iran_name, "_", meso_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  return(models)
}

#GENERATE SIMPLIFIED 4-WAY MODELS (fewer combinations)

generate_4way_models_simplified <- function() {
  models <- list()
  
  key_levant <- list(
    "Levant_Israel" = list(Levant = "Israel_PPNB.AG"),
    "Levant_Merged" = list(Levant = levant_merged)
  )
  
  key_anatolia <- list(
    "Anatolia_Barcin" = list(Anatolia = "Turkey_Marmara_Barcin_N.SG"),
    "Anatolia_Catalhoyuk" = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG")
  )
  
  key_iran <- list(
    "Iran_HajjiFiruz" = list(Iran = "Iran_HajjiFiruz_N.AG"),
    "Iran_ZagrosMerged" = list(Iran = iran_zagros_merged),
    "Iran_SehGabi_C" = list(Iran = "Iran_SehGabi_C.AG")
  )
  
  key_mesopotamia <- list(
    "Meso_Shanidar" = list(Mesopotamia = "Iraq_Shanidar.AG"),
    "Meso_PPNA" = list(Mesopotamia = "Iraq_PPNA.AG"),
    "Meso_IraqMerged" = list(Mesopotamia = iraq_merged),
    "Meso_Cayonu" = list(Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"),
    "Meso_SEAnatoliaMerged" = list(Mesopotamia = se_anatolia_merged),
    "Meso_Lazaridis_Full" = list(Mesopotamia = mesopotamia_lazaridis_full)
  )
  
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #4-way: Levant + Anatolia + Iran + Mesopotamia
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (iran_name in names(key_iran)) {
        for (meso_name in names(key_mesopotamia)) {
          source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], 
                             key_iran[[iran_name]], key_mesopotamia[[meso_name]])
          all_pops <- unlist(source_groups)
          if (!has_overlap(all_pops)) {
            model_name <- paste0("4way_", lev_name, "_", anat_name, "_", iran_name, "_", meso_name)
            models[[model_name]] <- list(source_groups = source_groups, label = model_name)
          }
        }
      }
    }
  }
  
  return(models)
}


#GENERATE ALL 4-WAY MODEL COMBINATIONS

generate_4way_models <- function() {
  models <- list()
  
  #Combine single and merged populations
  all_levant <- c(single_levant, merged_levant)
  all_anatolia <- c(single_anatolia, merged_anatolia)
  all_iran <- c(single_iran, merged_iran)
  all_mesopotamia <- c(single_mesopotamia, merged_mesopotamia, merged_iran_iraq)
  
  #Helper function to check for population overlap
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #4-way: Levant + Anatolia + Iran + Mesopotamia
  for (lev_name in names(all_levant)) {
    for (anat_name in names(all_anatolia)) {
      for (iran_name in names(all_iran)) {
        for (meso_name in names(all_mesopotamia)) {
          sources <- c(all_levant[[lev_name]], all_anatolia[[anat_name]], 
                       all_iran[[iran_name]], all_mesopotamia[[meso_name]])
          # Skip if any population appears twice
          if (!has_overlap(sources)) {
            model_name <- paste0("4way_", lev_name, "_", anat_name, "_", iran_name, "_", meso_name)
            models[[model_name]] <- list(
              sources = sources,
              label = model_name
            )
          }
        }
      }
    }
  }
  
  return(models)
}


#GENERATE PRIORITY MODELS (subset for testing)
#Now using source_groups format for proper population merging


generate_priority_models <- function() {
  models <- list()
  
  #3-WAY MODELS: Levant + Anatolia + Iran

  
  #Model 1: Levant_merged + Barcin + IranZagros_merged
  models[["3way_LevantMerged_Barcin_IranZagros"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = iran_zagros_merged
    ),
    label = "3way_LevantMerged_Barcin_IranZagros"
  )
  
  #Model 2: Jordan_PPNB + Barcin + HajjiFiruz (single pops)
  models[["3way_JordanPPNB_Barcin_HajjiFiruz"]] <- list(
    source_groups = list(
      Levant = "Jordan_PPNB.AG",
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG"
    ),
    label = "3way_JordanPPNB_Barcin_HajjiFiruz"
  )
  
  #Model 3: Levant_merged + Catalhoyuk + SehGabi (Chalcolithic Iran)
  models[["3way_LevantMerged_Catalhoyuk_SehGabi"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Iran = "Iran_SehGabi_C.AG"
    ),
    label = "3way_LevantMerged_Catalhoyuk_SehGabi"
  )
  

  #3-WAY MODELS: Levant + Anatolia + Mesopotamia

  #Model 4: Levant_merged + Barcin + Iraq_Shanidar
  models[["3way_LevantMerged_Barcin_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "3way_LevantMerged_Barcin_Shanidar"
  )
  
  #Model 5: Levant_merged + Barcin + Iraq_PPNA
  models[["3way_LevantMerged_Barcin_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "3way_LevantMerged_Barcin_PPNA"
  )
  
  # Model 6: Levant_merged + Barcin + Iraq_merged
  models[["3way_LevantMerged_Barcin_IraqMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = iraq_merged
    ),
    label = "3way_LevantMerged_Barcin_IraqMerged"
  )
  
  # Model 7: Levant_merged + Barcin + Cayonu (SE Anatolia = Upper Mesopotamia)
  models[["3way_LevantMerged_Barcin_Cayonu"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"
    ),
    label = "3way_LevantMerged_Barcin_Cayonu"
  )
  
  # Model 8: Levant_merged + Barcin + Mardin
  models[["3way_LevantMerged_Barcin_Mardin"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = "Turkey_Southeast_Mardin_PPN.AG"
    ),
    label = "3way_LevantMerged_Barcin_Mardin"
  )
  
  # Model 9: Levant_merged + Barcin + SE_Anatolia_merged
  models[["3way_LevantMerged_Barcin_SEAnatoliaMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Mesopotamia = se_anatolia_merged
    ),
    label = "3way_LevantMerged_Barcin_SEAnatoliaMerged"
  )
  
  # Model 10: Levant_merged + Catalhoyuk + Shanidar
  models[["3way_LevantMerged_Catalhoyuk_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "3way_LevantMerged_Catalhoyuk_Shanidar"
  )
  
  # 3-WAY MODELS: Levant + Iran + Mesopotamia

  # Model 11: Levant_merged + HajjiFiruz + Shanidar
  models[["3way_LevantMerged_HajjiFiruz_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "3way_LevantMerged_HajjiFiruz_Shanidar"
  )
  
  # Model 12: Levant_merged + HajjiFiruz + PPNA
  models[["3way_LevantMerged_HajjiFiruz_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "3way_LevantMerged_HajjiFiruz_PPNA"
  )
  
  # Model 13: Levant_merged + IranZagros_merged + IraqMerged
  models[["3way_LevantMerged_IranZagros_IraqMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Iran = iran_zagros_merged,
      Mesopotamia = iraq_merged
    ),
    label = "3way_LevantMerged_IranZagros_IraqMerged"
  )
  
  # 3-WAY MODELS: Using Iran+Iraq clades as single source

  # Model 14: Levant_merged + Barcin + Iran-Iraq Zagros clade
  models[["3way_LevantMerged_Barcin_IranIraqZagros"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      ZagrosMeso = iran_iraq_full_zagros
    ),
    label = "3way_LevantMerged_Barcin_IranIraqZagros"
  )
  
  # 4-WAY MODELS: Levant + Anatolia + Iran + Mesopotamia

  # Model 15: Levant_merged + Barcin + HajjiFiruz + Shanidar
  models[["4way_LevantMerged_Barcin_HajjiFiruz_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_Shanidar"
  )
  
  # Model 16: Levant_merged + Barcin + HajjiFiruz + PPNA
  models[["4way_LevantMerged_Barcin_HajjiFiruz_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_PPNA"
  )
  
  # Model 17: Levant_merged + Barcin + HajjiFiruz + Cayonu
  models[["4way_LevantMerged_Barcin_HajjiFiruz_Cayonu"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_Cayonu"
  )
  
  # Model 18: Levant_merged + Barcin + HajjiFiruz + Mardin
  models[["4way_LevantMerged_Barcin_HajjiFiruz_Mardin"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Turkey_Southeast_Mardin_PPN.AG"
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_Mardin"
  )
  
  # Model 19: Levant_merged + Barcin + GanjDareh + Shanidar
  models[["4way_LevantMerged_Barcin_GanjDareh_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_GanjDareh_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_LevantMerged_Barcin_GanjDareh_Shanidar"
  )
  
  # Model 20: Levant_merged + Barcin + GanjDareh + PPNA
  models[["4way_LevantMerged_Barcin_GanjDareh_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_GanjDareh_N.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "4way_LevantMerged_Barcin_GanjDareh_PPNA"
  )
  
  # Model 21: Levant_merged + Barcin + IranZagros_merged + IraqMerged
  models[["4way_LevantMerged_Barcin_IranZagros_IraqMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = iran_zagros_merged,
      Mesopotamia = iraq_merged
    ),
    label = "4way_LevantMerged_Barcin_IranZagros_IraqMerged"
  )
  
  # Model 22: Levant_merged + Catalhoyuk + HajjiFiruz + Shanidar
  models[["4way_LevantMerged_Catalhoyuk_HajjiFiruz_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_LevantMerged_Catalhoyuk_HajjiFiruz_Shanidar"
  )
  
  # Model 23: Levant_merged + Catalhoyuk + HajjiFiruz + PPNA
  models[["4way_LevantMerged_Catalhoyuk_HajjiFiruz_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "4way_LevantMerged_Catalhoyuk_HajjiFiruz_PPNA"
  )
  
  # Model 24: Levant_merged + Catalhoyuk + SehGabi (Chl) + PPNA
  models[["4way_LevantMerged_Catalhoyuk_SehGabi_PPNA"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Iran = "Iran_SehGabi_C.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "4way_LevantMerged_Catalhoyuk_SehGabi_PPNA"
  )
  
  # Model 25: Levant_merged + Catalhoyuk + SehGabi (Chl) + Shanidar
  models[["4way_LevantMerged_Catalhoyuk_SehGabi_Shanidar"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Central_Catalhoyuk_N.SG",
      Iran = "Iran_SehGabi_C.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_LevantMerged_Catalhoyuk_SehGabi_Shanidar"
  )
  
  # Model 26: Jordan_PPNB + Barcin + HajjiFiruz + Shanidar (single pops)
  models[["4way_JordanPPNB_Barcin_HajjiFiruz_Shanidar"]] <- list(
    source_groups = list(
      Levant = "Jordan_PPNB.AG",
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_JordanPPNB_Barcin_HajjiFiruz_Shanidar"
  )
  
  # Model 27: Jordan_PPNB + Barcin + HajjiFiruz + PPNA
  models[["4way_JordanPPNB_Barcin_HajjiFiruz_PPNA"]] <- list(
    source_groups = list(
      Levant = "Jordan_PPNB.AG",
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_PPNA.AG"
    ),
    label = "4way_JordanPPNB_Barcin_HajjiFiruz_PPNA"
  )
  
  # Model 28: Israel_PPNB + Barcin + HajjiFiruz + Shanidar
  models[["4way_IsraelPPNB_Barcin_HajjiFiruz_Shanidar"]] <- list(
    source_groups = list(
      Levant = "Israel_PPNB.AG",
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = "Iraq_Shanidar.AG"
    ),
    label = "4way_IsraelPPNB_Barcin_HajjiFiruz_Shanidar"
  )
  
  # Model 29: Levant_merged + Barcin + HajjiFiruz + IraqMerged
  models[["4way_LevantMerged_Barcin_HajjiFiruz_IraqMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = iraq_merged
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_IraqMerged"
  )
  
  # Model 30: Levant_merged + Barcin + HajjiFiruz + SEAnatoliaMerged
  models[["4way_LevantMerged_Barcin_HajjiFiruz_SEAnatoliaMerged"]] <- list(
    source_groups = list(
      Levant = levant_merged,
      Anatolia = "Turkey_Marmara_Barcin_N.SG",
      Iran = "Iran_HajjiFiruz_N.AG",
      Mesopotamia = se_anatolia_merged
    ),
    label = "4way_LevantMerged_Barcin_HajjiFiruz_SEAnatoliaMerged"
  )
  
  return(models)
}


#MAIN ANALYSIS FUNCTION


run_all_analyses <- function(models, outgroup_sets, output_prefix = "qpAdm_results") {
  
  all_results <- list()
  summary_tables <- list()
  
  for (out_name in names(outgroup_sets)) {
    cat("# Running analysis with outgroup set:", out_name, "\n")

    outgroups <- outgroup_sets[[out_name]]
    results_this_outgroup <- list()
    
    for (model_name in names(models)) {
      model <- models[[model_name]]
      
      #Create unique f2 directory for this combination
      f2_dir_name <- paste0("f2_blocks_", out_name, "_", model_name)
      
      #Run qpAdm with source_groups (properly merged populations)
      result <- run_qpadm_analysis(
        target = target_population,
        source_groups = model$source_groups,
        right_pops = outgroups,
        label = paste0(model$label, "_", out_name),
        f2_directory = f2_dir_name
      )
      
      #Store result
      result$summary$Outgroup_set <- out_name
      results_this_outgroup[[model_name]] <- result
    }
    
    all_results[[out_name]] <- results_this_outgroup
    
    #Create summary table for this outgroup set
    summaries <- lapply(results_this_outgroup, function(x) x$summary)
    summary_table <- bind_rows(summaries)
    summary_tables[[out_name]] <- summary_table
    
    #Save intermediate results
    saveRDS(results_this_outgroup, file = paste0(results_dir, output_prefix, "_", out_name, ".rds"))
    write_xlsx(summary_table, path = paste0(results_dir, output_prefix, "_", out_name, ".xlsx"))
  }
  
  #Combine all summary tables
  combined_summary <- bind_rows(summary_tables)
  
  #Save combined results
  saveRDS(all_results, file = paste0(results_dir, output_prefix, "_all.rds"))
  write_xlsx(combined_summary, path = paste0(results_dir, output_prefix, "_combined.xlsx"))
  
  return(list(
    all_results = all_results,
    summary_tables = summary_tables,
    combined_summary = combined_summary
  ))
}


#Filter and rank the results 

analyze_results <- function(combined_summary) {
  
  #Filter for feasible models with p > 0.05
  #Feasible is now character "TRUE"/"FALSE"
  passing_models <- combined_summary %>%
    filter(Feasible == "TRUE", p > 0.05) %>%
    arrange(desc(p))
  
  # Summary statistics
  cat("# Results summary\n")

  cat("Total models tested:", nrow(combined_summary), "\n")
  cat("Models with p > 0.05 and feasible weights:", nrow(passing_models), "\n\n")
  
  #Best models per outgroup set
  cat("Best models per outgroup set:\n")

  best_by_outgroup <- passing_models %>%
    group_by(Outgroup_set) %>%
    slice_max(p, n = 5) %>%
    select(Model, p, Feasible, SNPs_used, everything())
  
  print(as.data.frame(best_by_outgroup))
  
  #Most consistent models across outgroup sets
  cat("\n\nModels passing in multiple outgroup sets:\n")

  model_counts <- passing_models %>%
    mutate(base_model = gsub("_out[0-9]+$", "", Model)) %>%
    group_by(base_model) %>%
    summarise(
      n_passing = n(),
      mean_p = mean(p),
      min_p = min(p),
      max_p = max(p)
    ) %>%
    filter(n_passing >= 3) %>%
    arrange(desc(n_passing), desc(mean_p))
  
  print(as.data.frame(model_counts))
  
  return(list(
    passing_models = passing_models,
    best_by_outgroup = best_by_outgroup,
    robust_models = model_counts
  ))
}


#Additionnal simplified 3-way models (testing different Levant sources and Iran options)

generate_3way_models_simplified_v2 <- function() {
  models <- list()
  
  #LEVANT: Test individual populations
  key_levant <- list(
    "Levant_Jordan" = list(Levant = "Jordan_PPNB.AG"),
    "Levant_JordanBaja" = list(Levant = "Jordan_Baja_Late_PPNB.AG"),
    "Levant_Israel" = list(Levant = "Israel_PPNB.AG"),
    "Levant_JordanMerged" = list(Levant = c("Jordan_PPNB.AG", "Jordan_Baja_Late_PPNB.AG"))
  )
  
  #ANATOLIA: Keep simple
  key_anatolia <- list(
    "Anatolia_Barcin" = list(Anatolia = "Turkey_Marmara_Barcin_N.SG"),
    "Anatolia_Catalhoyuk" = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG")
  )
  
  #IRAN: All single populations + qpWave-based merges
  key_iran <- list(
    # Single populations
    "Iran_HajjiFiruz" = list(Iran = "Iran_HajjiFiruz_N.AG"),
    "Iran_GanjDareh" = list(Iran = "Iran_GanjDareh_N.AG"),
    "Iran_TepeAbdul" = list(Iran = "Iran_TepeAbdulHosein_N.SG"),
    "Iran_Wezmeh" = list(Iran = "Iran_Wezmeh_N.SG"),
    "Iran_SehGabi" = list(Iran = "Iran_SehGabi_C.AG"),
    # qpWave-based merges
    "Iran_GanjDareh_Tepe" = list(Iran = c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG")),
    "Iran_Tepe_Wezmeh" = list(Iran = c("Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG")),
    "Iran_ZagrosMerged" = list(Iran = iran_zagros_merged)
  )
  
  #MESOPOTAMIA: Core options
  key_mesopotamia <- list(
    "Meso_Shanidar" = list(Mesopotamia = "Iraq_Shanidar.AG"),
    "Meso_PPNA" = list(Mesopotamia = "Iraq_PPNA.AG"),
    "Meso_IraqMerged" = list(Mesopotamia = iraq_merged),
    "Meso_Cayonu" = list(Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"),
    "Meso_Mardin" = list(Mesopotamia = "Turkey_Southeast_Mardin_PPN.AG")
  )
  
  #Iran-Iraq merged clades (qpWave-based)
  key_iran_iraq <- list(
    "ZagrosMeso_GanjDareh_Shanidar" = list(ZagrosMeso = c("Iran_GanjDareh_N.AG", "Iraq_Shanidar.AG")),
    "ZagrosMeso_Tepe_Shanidar" = list(ZagrosMeso = c("Iran_TepeAbdulHosein_N.SG", "Iraq_Shanidar.AG")),
    "ZagrosMeso_IranZagros_IraqMerged" = list(ZagrosMeso = c(iran_zagros_merged, iraq_merged))
  )
  
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #3-way: Levant + Anatolia + Iran (all Iran single/merged options)
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (iran_name in names(key_iran)) {
        source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], key_iran[[iran_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", iran_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  #3-way: Levant + Anatolia + Mesopotamia
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (meso_name in names(key_mesopotamia)) {
        source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], key_mesopotamia[[meso_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", meso_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  #3-way: Levant + Anatolia + Iran-Iraq merged clade
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (zirc_name in names(key_iran_iraq)) {
        source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], key_iran_iraq[[zirc_name]])
        all_pops <- unlist(source_groups)
        if (!has_overlap(all_pops)) {
          model_name <- paste0("3way_", lev_name, "_", anat_name, "_", zirc_name)
          models[[model_name]] <- list(source_groups = source_groups, label = model_name)
        }
      }
    }
  }
  
  return(models)
}


#Additionnal simplified 4-way models (testing different Levant and Iran options)


generate_4way_models_simplified_v2 <- function() {
  models <- list()
  
  #LEVANT: Test individual populations
  key_levant <- list(
    "Levant_Jordan" = list(Levant = "Jordan_PPNB.AG"),
    "Levant_JordanBaja" = list(Levant = "Jordan_Baja_Late_PPNB.AG"),
    "Levant_Israel" = list(Levant = "Israel_PPNB.AG"),
    "Levant_JordanMerged" = list(Levant = c("Jordan_PPNB.AG", "Jordan_Baja_Late_PPNB.AG"))
  )
  
  #ANATOLIA: Keep simple
  key_anatolia <- list(
    "Anatolia_Barcin" = list(Anatolia = "Turkey_Marmara_Barcin_N.SG"),
    "Anatolia_Catalhoyuk" = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG")
  )
  
  #IRAN: All single populations + qpWave-based merges
  key_iran <- list(
    # Single populations
    "Iran_HajjiFiruz" = list(Iran = "Iran_HajjiFiruz_N.AG"),
    "Iran_GanjDareh" = list(Iran = "Iran_GanjDareh_N.AG"),
    "Iran_TepeAbdul" = list(Iran = "Iran_TepeAbdulHosein_N.SG"),
    "Iran_Wezmeh" = list(Iran = "Iran_Wezmeh_N.SG"),
    "Iran_SehGabi" = list(Iran = "Iran_SehGabi_C.AG"),
    # qpWave-based merges
    "Iran_GanjDareh_Tepe" = list(Iran = c("Iran_GanjDareh_N.AG", "Iran_TepeAbdulHosein_N.SG")),
    "Iran_Tepe_Wezmeh" = list(Iran = c("Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG")),
    "Iran_ZagrosMerged" = list(Iran = iran_zagros_merged)
  )
  
  #MESOPOTAMIA: Core options (excluding Iran-Iraq merged since we have separate Iran)
  key_mesopotamia <- list(
    "Meso_Shanidar" = list(Mesopotamia = "Iraq_Shanidar.AG"),
    "Meso_PPNA" = list(Mesopotamia = "Iraq_PPNA.AG"),
    "Meso_IraqMerged" = list(Mesopotamia = iraq_merged),
    "Meso_Cayonu" = list(Mesopotamia = "Turkey_Southeast_Cayonu_PPN.SG"),
    "Meso_Mardin" = list(Mesopotamia = "Turkey_Southeast_Mardin_PPN.AG"),
    "Meso_SEAnatoliaMerged" = list(Mesopotamia = se_anatolia_merged)
  )
  
  has_overlap <- function(...) {
    pops <- unlist(list(...))
    length(pops) != length(unique(pops))
  }
  
  #4-way: Levant + Anatolia + Iran + Mesopotamia
  for (lev_name in names(key_levant)) {
    for (anat_name in names(key_anatolia)) {
      for (iran_name in names(key_iran)) {
        for (meso_name in names(key_mesopotamia)) {
          source_groups <- c(key_levant[[lev_name]], key_anatolia[[anat_name]], 
                             key_iran[[iran_name]], key_mesopotamia[[meso_name]])
          all_pops <- unlist(source_groups)
          if (!has_overlap(all_pops)) {
            model_name <- paste0("4way_", lev_name, "_", anat_name, "_", iran_name, "_", meso_name)
            models[[model_name]] <- list(source_groups = source_groups, label = model_name)
          }
        }
      }
    }
  }
  
  return(models)
}


#Run the analysis 


#Generate models (choose which set to run) 
cat("Generating model combinations\n\n")

#Option 1: Priority models (hand-selected)
#priority_models <- generate_priority_models()

#Option 2: Simplified 3-way and 4-way (key populations only)
simplified_3way <- generate_3way_models_simplified_v2()
simplified_4way <- generate_4way_models_simplified_v2()
simplified_models <- c(simplified_3way, simplified_4way)

#Option 3: Full analysis (all combinations (VERY LARGE), use with caution)
# all_3way_models <- generate_3way_models()
# all_4way_models <- generate_4way_models()
# all_models <- c(all_3way_models, all_4way_models)

#Show model counts
cat("Model counts:\n")
#cat("  Priority models:", length(priority_models), "\n")
cat("  Simplified 3-way:", length(simplified_3way), "\n")
cat("  Simplified 4-way:", length(simplified_4way), "\n")
cat("  Simplified total:", length(simplified_models), "\n")
# cat("  Full 3-way:", length(all_3way_models), "\n")
# cat("  Full 4-way:", length(all_4way_models), "\n")
# cat("  Full total:", length(all_models), "\n")
cat("\n")
cat("Number of outgroup sets:", length(outgroup_sets), "\n")
#cat("Total analyses with priority models:", length(priority_models) * length(outgroup_sets), "\n")
cat("Total analyses with simplified models:", length(simplified_models) * length(outgroup_sets), "\n")
cat("\n")

#Choose which model to run 
#Change 'models_to_run' to: priority_models, simplified_models, or all_models
models_to_run <- simplified_models  #Start with this

cat("Running", length(models_to_run), "models x", length(outgroup_sets), "outgroup sets =", 
    length(models_to_run) * length(outgroup_sets), "total analyses\n")
cat("f2 directories will be cleaned up after each analysis to save disk space.\n\n")

#Run analyses
cat("Starting qpAdm analyses...\n\n")

results <- run_all_analyses(
  models = models_to_run,
  outgroup_sets = outgroup_sets,
  output_prefix = "qpAdm_LevantChL_simplified_v2"
)

analyze_results2 <- function(combined_summary) {

  check_significant_proportions <- function(row) {
    # Find all columns ending in "_prop" and "_se"
    prop_cols <- grep("_prop$", names(row), value = TRUE)
    se_cols <- grep("_se$", names(row), value = TRUE)
    
    # For each source, check if |proportion| > SE
    all_significant <- TRUE
    for (prop_col in prop_cols) {
      # Get corresponding SE column
      source_name <- sub("_prop$", "", prop_col)
      se_col <- paste0(source_name, "_se")
      
      if (se_col %in% names(row)) {
        prop_val <- as.numeric(row[[prop_col]])
        se_val <- as.numeric(row[[se_col]])
        
        #Check if proportion is significantly different from zero
        #Using absolute value since negative proportions can occur
        if (!is.na(prop_val) && !is.na(se_val)) {
          if (abs(prop_val) <= se_val) {
            all_significant <- FALSE
            break
          }
        }
      }
    }
    return(all_significant)
  }
  
  #Apply the significance check to each row
  combined_summary$Significant_props <- sapply(1:nrow(combined_summary), function(i) {
    check_significant_proportions(combined_summary[i, ])
  })
  
  #Filter for feasible models with p > 0.05 AND all proportions significant
  passing_models <- combined_summary %>%
    filter(Feasible == "TRUE", 
           p > 0.05,
           Significant_props == TRUE) %>%
    arrange(desc(p))
  
  #Also create a secondary list: feasible + p>0.05 but some non-significant proportions
  marginal_models <- combined_summary %>%
    filter(Feasible == "TRUE", 
           p > 0.05,
           Significant_props == FALSE) %>%
    arrange(desc(p))
  
  #Summary statistics
  cat("# Results summary\n")

  cat("Total models tested:", nrow(combined_summary), "\n")
  cat("Models with p > 0.05 and feasible weights:", 
      nrow(filter(combined_summary, Feasible == "TRUE", p > 0.05)), "\n")
  cat("Models with ALL proportions significant (prop > SE):", nrow(passing_models), "\n")
  cat("Models with SOME non-significant proportions:", nrow(marginal_models), "\n\n")
  
  # Best models per outgroup set
  cat("Best passing models (all proportions significant)\n")

  if (nrow(passing_models) > 0) {
    best_by_outgroup <- passing_models %>%
      group_by(Outgroup_set) %>%
      slice_max(p, n = 5) %>%
      select(Model, p, Feasible, SNPs_used, everything())
    
    print(as.data.frame(best_by_outgroup))
  } else {
    cat("No models passed all criteria.\n")
  }
  
  cat("\n\n Marginal models (some proportions have SE > prop)\n")

  if (nrow(marginal_models) > 0) {
    marginal_by_outgroup <- marginal_models %>%
      group_by(Outgroup_set) %>%
      slice_max(p, n = 3) %>%
      select(Model, p, Feasible, SNPs_used, everything())
    
    print(as.data.frame(marginal_by_outgroup))
  } else {
    cat("No marginal models.\n")
  }
  
  #Save results
  write_xlsx(passing_models, path = paste0(results_dir, "passing_models_significant.xlsx"))
  write_xlsx(marginal_models, path = paste0(results_dir, "marginal_models_nonsignificant.xlsx"))
  
  cat("\n\nResults saved to:\n")
  cat("  - passing_models_significant.xlsx (all proportions > SE)\n")
  cat("  - marginal_models_nonsignificant.xlsx (some proportions <= SE)\n")
  
  return(list(
    passing = passing_models,
    marginal = marginal_models,
    all = combined_summary
  ))
}

results_import <- read_excel("results_aDNA/qpAdm_LevantChL_simplified_v2_combined.xlsx")
#Analyze results
analysis <- analyze_results(results$combined_summary)

analysis2 <- analyze_results2(results_import)

#Save final analysis
saveRDS(analysis2, file = paste0(results_dir, "2qpAdm_analysis_simplified_summary_v2.rds"))
write_xlsx(analysis2$passing, path = paste0(results_dir, "2qpAdm_passing_simplified_models_v2.xlsx"))
write_xlsx(analysis2$marginal, path = paste0(results_dir, "2qpAdm_marginal_simplified_models_v2.xlsx"))

cat("\n\nAnalysis complete! Results saved to:", results_dir, "\n")





#Bidirectional gene flow analysis 

#Use outgroup set 3 (most SNPs for these populations)
bidir_outgroup <- out3

#MESOPOTAMIA AS TARGET (without merging)
cat("Levant contribution to Mesopotamian populations\n\n")

bidir_meso_models <- list(
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_GanjDareh_N.AG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_GanjDareh_N.AG", Levant = "Jordan_PPNB.AG")),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = "Jordan_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = "Jordan_PPNB.AG")),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = "Jordan_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = "Jordan_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_GanjDareh_N.AG", Levant = "Israel_PPNB.AG")),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Levant = "Israel_PPNB.AG"))
)

bidir_meso_results <- list()
for (m in bidir_meso_models) {
  label <- paste0("Bidir_", m$target, "_from_zagros_", m$sources$Zagros, "_from_levant_", m$sources$Levant)
  result <- run_qpadm_analysis(m$target, m$sources, bidir_outgroup, label)
  bidir_meso_results[[label]] <- result$summary
  cat(sprintf("%s: p=%.3f, %s=%.1f%%, %s=%.1f%%\n", 
              m$target, result$summary$P_value,
              names(m$sources)[1], result$summary[[paste0(names(m$sources)[1], "_prop")]]*100,
              names(m$sources)[2], result$summary[[paste0(names(m$sources)[2], "_prop")]]*100))
}

#MESOPOTAMIA AS TARGET (with merging)
cat("\nWith merged Levant source\n\n")

bidir_meso_merged <- list(
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_GanjDareh_N.AG", Levant = levant_merged)),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_GanjDareh_N.AG", Levant = levant_merged)),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = levant_merged)),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_HajjiFiruz_N.AG", Levant = levant_merged)),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = levant_merged)),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = "Iran_TepeAbdulHosein_N.SG", Levant = levant_merged)),
  list(target = "Iraq_Shanidar.AG", 
       sources = list(Zagros = iran_ganjdareh_tepe, Levant = levant_merged)),
  list(target = "Iraq_PPNA.AG", 
       sources = list(Zagros = iran_ganjdareh_tepe, Levant = levant_merged))
)

for (m in bidir_meso_merged) {
  label <- paste0("Bidir_merged_", m$target, "from_zagros_", m$sources$Zagros)
  result <- run_qpadm_analysis(m$target, m$sources, bidir_outgroup, label)
  bidir_meso_results[[label]] <- result$summary
  cat(sprintf("%s: p=%.3f\n", m$target, result$summary$P_value))
}

#LEVANT PPNB AS TARGET (Mesopotamian into Levant)
cat("\n Mesopotamian contribution to Levant PPNB\n\n")

bidir_levant_models <- list(
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = "Iraq_Shanidar.AG")),
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = "Iraq_Shanidar.AG")),
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = "Iraq_PPNA.AG")),
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = "Iraq_PPNA.AG")),
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = iraq_merged)),
  list(target = "Israel_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = iraq_merged)),
  
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = "Iraq_Shanidar.AG")),
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = "Iraq_Shanidar.AG")),
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = "Iraq_PPNA.AG")),
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = "Iraq_PPNA.AG")),
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Central_Catalhoyuk_N.SG", Mesopotamia = iraq_merged)),
  list(target = "Jordan_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = iraq_merged)),

  list(target = "Jordan_Baja_Late_PPNB.AG", 
       sources = list(Anatolia = "Turkey_Marmara_Barcin_N.SG", Mesopotamia = "Iraq_PPNA.AG"))
)

bidir_levant_results <- list()
for (m in bidir_levant_models) {
  label <- paste0("Bidir_", m$target, "_from_Anatolia_", m$sources$Anatolia, "_from_Meso_", m$sources$Mesopotamia)
  result <- run_qpadm_analysis(m$target, m$sources, bidir_outgroup, label)
  bidir_levant_results[[label]] <- result$summary
  cat(sprintf("%s: p=%.3f\n", m$target, result$summary$P_value))
}

#SAVE BIDIRECTIONAL RESULTS

#Function to safely extract single values from results
safe_extract <- function(x) {
  if (is.null(x)) return(NA)
  if (length(x) > 1) return(x[1])  # Take first element if vector
  return(x)
}

#Function to convert a single result to a clean dataframe row
clean_result <- function(res, model_id) {
  # Extract all values safely
  df <- data.frame(
    Model_ID = model_id,
    Model = safe_extract(res$Model),
    Target = safe_extract(res$Target),
    Sources = safe_extract(res$Sources),
    N_sources = safe_extract(res$N_sources),
    p = safe_extract(res$p),
    Feasible = safe_extract(res$Feasible),
    SNPs_used = safe_extract(res$SNPs_used),
    stringsAsFactors = FALSE
  )
  
  # Add proportion columns if they exist
  prop_cols <- grep("_prop$|_se$", names(res), value = TRUE)
  for (col in prop_cols) {
    df[[col]] <- safe_extract(res[[col]])
  }
  
  return(df)
}

#Process Mesopotamia target results
meso_clean <- lapply(names(bidir_meso_results), function(nm) {
  clean_result(bidir_meso_results[[nm]], nm)
})

#Process Levant target results  
levant_clean <- lapply(names(bidir_levant_results), function(nm) {
  clean_result(bidir_levant_results[[nm]], nm)
})

#Combine using plyr::rbind.fill to handle different columns
bidir_meso_df <- plyr::rbind.fill(meso_clean)
bidir_meso_df$Direction <- "Levant_to_Mesopotamia"

bidir_levant_df <- plyr::rbind.fill(levant_clean)
bidir_levant_df$Direction <- "Mesopotamia_to_Levant"

bidir_all <- plyr::rbind.fill(bidir_meso_df, bidir_levant_df)

#Reorder columns for clarity
col_order <- c("Direction", "Model_ID", "Target", "Sources", "N_sources", "p", "Feasible", "SNPs_used")
other_cols <- setdiff(names(bidir_all), col_order)
bidir_all <- bidir_all[, c(col_order, other_cols)]

write_xlsx(bidir_all, file.path(results_dir, "Bidirectional_GeneFlow_Results.xlsx"))
cat("\nBidirectional results saved.\n")
cat("Total models:", nrow(bidir_all), "\n")
cat("Passing models (p > 0.05 & Feasible):", sum(bidir_all$p > 0.05 & bidir_all$Feasible == TRUE, na.rm = TRUE), "\n")


