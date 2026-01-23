library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#Load and preprocess AADR data
aadr_path <- "../AADR_v62_HO/"
prefix <- file.path(aadr_path, "v62.0_HO_public")
results_dir <- "../results_aDNA/"
#dir.create(file.path(results_dir, "pca_analysis"), showWarnings = FALSE)

#Read population information
ind_file <- paste0(prefix, ".ind")
anno_file <- paste0(prefix, ".anno")

#Load individual and annotation data
individuals <- read_table(ind_file, col_names = c("ID", "Sex", "Population"))

prepare_population_lists <- function(ind_data) {
  
  cat("Reading .ind file and categorizing populations...\n")
  
  #Get unique populations
  all_pops <- unique(ind_data$Population)
  
  #Define MODERN populations for computing PCs
  #These should be present-day populations with good coverage
  modern_pops <- all_pops[!grepl("\\.AG|\\.SG|\\.DG|_", all_pops)]
  print(modern_pops)
  
  #Common modern populations in HO dataset
  modern_pops_expected <- c(
    #Middle East
    "Lebanese.HO", "Druze.HO", "Palestinian.HO", "Jordanian.HO", "Syrian.HO",
    "BedouinA.HO", "BedouinB.HO", "Saudi.HO", "Yemenite.HO", "Egyptian.HO",
    
    #Anatolia/Caucasus
    "Turkish.HO", "Armenian.HO", "Georgian.HO", "Abkhasian.HO", "Chechen.HO",
    "Lezgin.HO", "Adygei.HO", "Balkar.HO", "North_Ossetian.HO",
    
    #Iran/Central Asia
    "Iranian.HO", "Iranian_Bandari.HO", "Balochi.HO", "Brahui.HO", "Makrani.HO",
    "Pathan.HO", "Kalash.HO", "Burusho.HO", "Hazara.HO", "Uygur.HO",
    
    #Europe (for context)
    "Sardinian.HO", "Basque.HO", "French.HO", "Italian_North.HO", "Spanish.HO",
    "Orcadian.HO", "Russian.HO", "Hungarian.HO", "Czech.HO", "Greek.HO",
    
    #Africa (outgroups)
    "Yoruba.HO", "Mbuti.HO", "San.HO", "Ethiopian_Jew.HO", "Somali.HO",
    
    #East Asia (for ANE detection)
    "Han.HO", "She.HO", "Japanese.HO", "Yakut.HO", "Mongola.HO"
  )
  
  #Get actual modern populations available
  modern_pops_available <- intersect(modern_pops_expected, modern_pops)
  
  #Define ANCIENT populations to project
  ancient_pops <- all_pops[grepl("\\.AG|\\.SG|\\.DG|_UP|_HG|_N|_C|_BA|_EBA|_MBA|_LBA|_IA", all_pops)]
  #print(ancient_pops)
  
  #Focus on Near Eastern ancient populations relevant to your study
  ancient_near_east <- c(
    "Israel_C.AG",
    "Israel_Natufian.AG", "Israel_PPNB.AG", "Jordan_PPNB.AG", "Jordan_Baja_Late_PPNB.AG",
    "Turkey_Marmara_Barcin_N.SG", "Turkey_Central_Catalhoyuk_N.SG", "Turkey_Southeast_Cayonu_PPN.SG",
    "Turkey_Southeast_Mardin_PPN.AG", "Turkey_Central_Boncuklu_PPN.WGC.SG",
    "Iran_GanjDareh_N.AG", "Iran_HajjiFiruz_N.AG", "Iran_SehGabi_C.AG", 
    "Iran_TepeAbdulHosein_N.SG", "Iran_Wezmeh_N.SG",
    "Iraq_PPNA.AG", "Iraq_Shanidar.AG",
    "Georgia_Kotias_Mesolithic.SG", "Armenia_Aknashen_N.AG"
  )
  
  ancient_near_east_available <- intersect(ancient_near_east, ancient_pops)
  
  return(list(
    modern = modern_pops_available,
    ancient = ancient_near_east_available,
    n_modern = length(modern_pops_available),
    n_ancient = length(ancient_near_east_available)
  ))
}

#Get population lists
pop_lists <- prepare_population_lists(individuals)

cat(sprintf("Found %d modern populations for PC calculation\n", pop_lists$n_modern))
cat(sprintf("Found %d ancient populations to project\n", pop_lists$n_ancient))

create_smartpca_params <- function(prefix, modern_pops, ancient_pops, results_dir) {
  
  modern_file <- file.path(results_dir, "pca_analysis/modern_populations.txt")
  writeLines(modern_pops, modern_file)

  ancient_file <- file.path(results_dir, "pca_analysis/ancient_populations.txt")
  writeLines(ancient_pops, ancient_file)
  
  param_content <- sprintf("
genotypename:    %s.geno
snpname:         %s.snp
indivname:       %s.ind
evecoutname:     %s/pca_analysis/HO_pca.evec
evaloutname:     %s/pca_analysis/HO_pca.eval
poplistname:     %s
lsqproject:      YES
numoutevec:      10
numoutlieriter:  0
numthreads:      8
", prefix, prefix, prefix, results_dir, results_dir, modern_file)
  
  param_file <- "smartpca_HO.par"
  writeLines(param_content, param_file)
  
  return(param_file)
}

#Create parameter file
param_file <- create_smartpca_params(prefix, pop_lists$modern, 
                                     pop_lists$ancient, results_dir)

parse_smartpca_output <- function(results_dir, ancient_pops) {
  
  evec_file <- file.path(results_dir, "pca_analysis/HO_pca.evec")
  eval_file <- file.path(results_dir, "pca_analysis/HO_pca.eval")
  
  evec <- read.table(evec_file, header = FALSE, stringsAsFactors = FALSE)
  eval <- scan(eval_file, quiet = TRUE)
  
  #PCs are numeric columns
  is_num <- sapply(evec, function(x) all(!is.na(suppressWarnings(as.numeric(head(x, 20))))))
  pc_mat <- evec[, is_num]
  colnames(pc_mat) <- paste0("PC", seq_len(ncol(pc_mat)))
  
  #Non-numeric columns -> IDs and POPs
  non_num <- evec[, !is_num, drop = FALSE]
  if (ncol(non_num) == 2) {
    IDs  <- non_num[,1]
    POPs <- non_num[,2]
  } else {
    IDs  <- non_num[,1]
    POPs <- IDs
  }
  
  #Mark Ancient vs Modern
  Type <- ifelse(POPs %in% ancient_pops, "Ancient", "Modern")
  
  #Assemble dataframe
  pca_data <- data.frame(ID = IDs, pc_mat, Population = POPs, Type = Type)
  
  # Variance explained
  var_explained <- (eval / sum(eval)) * 100
  
  return(list(pca_data = pca_data,
              eigenvalues = eval,
              var_explained = var_explained))
}

#Parse output
pca_results <- parse_smartpca_output(results_dir, pop_lists$ancient)


create_pca_plot_zoom <- function(pca_results) {
  
  pca_data <- pca_results$pca_data
  var_explained <- pca_results$var_explained
  
  group_map <- c(
    "Israel_C.AG" = "Israel_C.AG",
    "Israel_Natufian.AG" = "Levant", "Israel_PPNB.AG" = "Levant", 
    "Jordan_PPNB.AG" = "Levant", "Jordan_Baja_Late_PPNB.AG" = "Levant",
    "Turkey_Marmara_Barcin_N.SG" = "Anatolia", "Turkey_Central_Catalhoyuk_N.SG" = "Anatolia",
    "Turkey_Central_Boncuklu_PPN.WGC.SG" = "Anatolia",
    "Turkey_Southeast_Cayonu_PPN.SG" = "Mesopotamia", "Turkey_Southeast_Mardin_PPN.AG" = "Mesopotamia",
    "Iraq_PPNA.AG" = "Mesopotamia", "Iraq_Shanidar.AG" = "Mesopotamia",
    "Iran_GanjDareh_N.AG" = "Zagros", "Iran_HajjiFiruz_N.AG" = "Zagros",
    "Iran_SehGabi_C.AG" = "Zagros", "Iran_TepeAbdulHosein_N.SG" = "Zagros", "Iran_Wezmeh_N.SG" = "Zagros",
    "Georgia_Kotias_Mesolithic.SG" = "Caucasus", "Armenia_Aknashen_N.AG" = "Caucasus"
  )
  
  # Add group column
  pca_data$Group <- ifelse(pca_data$Population %in% names(group_map),
                           group_map[pca_data$Population],
                           "Other")
  
  # Define periods
  pca_data$Period <- case_when(
    grepl("_C.AG|ChL", pca_data$Population) ~ "Chalcolithic",
    grepl("_N|PPN|Shanidar", pca_data$Population) ~ "Neolithic",
    grepl("Israel_Natufian.AG", pca_data$Population) ~ "Natufian",
    grepl("MBA", pca_data$Population) ~ "MBA",
    grepl("Jordan_EBA.AG", pca_data$Population) ~ "EBA",
    grepl("Georgia_Kotias_Mesolithic.SG", pca_data$Population) ~ "Mesolithic"
  )
  
  # Separate modern vs ancient
  modern_df <- dplyr::filter(pca_data, Type == "Modern")
  ancient_df <- dplyr::filter(pca_data, Type == "Ancient")
  
  
  euro_highlight <- c("Sardinian.HO", "Basque.HO", "French.HO", "Italian_North.HO", "Spanish.HO")
  euros <- modern_df %>%
    dplyr::filter(Population %in% euro_highlight) %>%
    dplyr::mutate(Group = "Modern European",   # color legend key
                  Period = "Modern")   
  
  # Define colors per group
  group_colors <- c(
    "Israel_C.AG" = "red",
    "Levant" = "darkgreen",
    "Anatolia" = "blue",
    "Zagros" = "orange",
    "Mesopotamia" = "purple",
    "Caucasus" = "brown",
    "Modern European" = "pink",
    "Other" = "gray60"
  )
  
  # Plot with zoom
  p <- ggplot() +
    
    geom_point(data = modern_df,
               aes(x = PC1, y = PC2),
               color = "gray90", size = 1, alpha = 0.2) +
    
    # Ancient grouped
    geom_point(data = ancient_df,
               aes(x = PC1, y = PC2, color = Group, shape = Period),
               size = 3, alpha = 0.9, stroke = 0.3) +
    
    geom_point(data = dplyr::filter(ancient_df, Group == "Israel_C.AG"),
               aes(PC1, PC2, color = Group, shape = Period), size = 4.2, show.legend = TRUE) +
    
    geom_point(data = euros,
               aes(PC1, PC2, color = Group, shape = Period),
               size = 2.6, alpha = 0.95, stroke = 0.9) +
    scale_color_manual(
      values = group_colors,
      breaks = c("Israel_C.AG","Levant","Anatolia","Zagros","Mesopotamia","Caucasus","Modern European"),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("Neolithic" = 16, "Chalcolithic" = 17, "Mesolithic" = 15, "MBA" = 14, "EBA" = 13, "Modern" = 4),
      breaks = c("Neolithic","Chalcolithic","Mesolithic","EBA","MBA","Modern"),
      drop = TRUE
    )  +
    labs(
      title = "PCA (Zoomed) of Ancient Near Eastern Populations",
      subtitle = "Ancient samples colored by region; moderns in gray",
      x = paste0("PC1"),
      y = paste0("PC2"),
      color = "Ancient Region",
      shape = "Ancient Period"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(fill = NA, color = "black")) + 
    scale_x_continuous(labels = function(x) gsub("-", "\u2212", scales::number_format()(x))) +
    scale_y_continuous(labels = function(y) gsub("-", "\u2212", scales::number_format()(y))) 
    
    #Zoom to focus on ancient variation
    #coord_cartesian(xlim = c(-0.025, 0.01), ylim = c(0.00, 0.02)) +
    #coord_cartesian(xlim = c(-0.025, 0.025), ylim = c(-0.025, 0.025))
  
  return(p)
}

pca_zoom <- create_pca_plot_zoom(pca_results)
ggsave(file.path(results_dir, "pca_analysis/Figure_PCA_HO_projected.png"),
       pca_zoom, width = 10, height = 7, dpi = 600)


