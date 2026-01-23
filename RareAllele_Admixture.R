library(tidyverse)
library(data.table)
library(admixtools)
library(ggplot2)
library(ggtext)
library(pheatmap)
library(writexl)

aadr_path <- "AADR_v62/"
prefix <- file.path(aadr_path, "v62.0_1240k_public")
prefix_merged <- paste0(prefix, "_merged")
results_dir <- "results_aDNA/"
dir.create(file.path(results_dir, "rare_allele_analysis"), showWarnings = FALSE)

target_pop <- "Israel_C.AG"
source_populations_merged <- c(
  "Levant_N_merged",
  "Iran_N_merged",
  "Iran_ChL_merged",
  "Anatolia_N_merged",
  "Mesopotamia_N_merged"
)

outgroup_populations <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Chukchi.DG",
  "Karitiana.DG",
  "Mbuti.DG"
)

#First create a parameter file for convertf
convertf_params <- "
genotypename:    v62.0_1240k_public_merged.geno
snpname:         v62.0_1240k_public_merged.snp
indivname:       v62.0_1240k_public_merged.ind
outputformat:    PACKEDPED
genotypeoutname: ancient_merged.bed
snpoutname:      ancient_merged.bim
indivoutname:    ancient_merged.fam
familynames:     NO
"

writeLines(convertf_params, "convertf_eigenstrat_to_plink.par")

#Run convertf (requires EIGENSOFT installed)
system("convertf -p convertf_eigenstrat_to_plink.par")

system("plink --bfile ancient_merged \\
        --maf 0.01 \\
        --max-maf 0.05 \\
        --make-bed \\
        --out ancient_rare_variants")

#Read the filtered rare variant data
fam_data <- fread("AADR_v62/ancient_rare_variants.fam", header = FALSE)
bim_data <- fread("AADR_v62/ancient_rare_variants.bim", header = FALSE)
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
colnames(bim_data) <- c("Chr", "SNP", "cM", "BP", "A1", "A2")

#Read .ind file to get population assignments
ind_data <- fread(paste0(prefix, ".ind"), header = FALSE)
colnames(ind_data) <- c("IID", "Sex", "Population")

fam_data <- fam_data %>%
  left_join(ind_data[, c("IID", "Population")], by = "IID")

calculate_rare_allele_sharing <- function(bed_file, fam_data, populations_of_interest) {
  
  library(BEDMatrix)
  
  bed_matrix <- BEDMatrix(bed_file)
  
  pop_indices <- which(fam_data$Population %in% populations_of_interest)
  
  if(length(pop_indices) == 0) {
    stop("No individuals found for specified populations")
  }
  
  geno_subset <- bed_matrix[pop_indices, ]
  fam_subset <- fam_data[pop_indices, ]
  
  n_ind <- nrow(geno_subset)
  n_snp <- ncol(geno_subset)
  
  cat(sprintf("Analyzing %d individuals and %d rare SNPs\n", n_ind, n_snp))
  
  pairwise_sharing <- matrix(0, n_ind, n_ind)
  rownames(pairwise_sharing) <- fam_subset$IID
  colnames(pairwise_sharing) <- fam_subset$IID
  
  unique_pops <- unique(fam_subset$Population)
  pop_sharing <- matrix(0, length(unique_pops), length(unique_pops))
  rownames(pop_sharing) <- unique_pops
  colnames(pop_sharing) <- unique_pops
  
  pb <- txtProgressBar(min = 0, max = n_snp, style = 3)
  
  for (snp_idx in 1:n_snp) {
    setTxtProgressBar(pb, snp_idx)
    
    geno <- as.numeric(geno_subset[, snp_idx])
    
    if (sum(!is.na(geno)) < 2 || length(unique(na.omit(geno))) <= 1) next
    
    rare_carriers <- which(geno > 0 & !is.na(geno))
    
    if (length(rare_carriers) > 1) {
      for (i in 1:(length(rare_carriers)-1)) {
        for (j in (i+1):length(rare_carriers)) {
          idx1 <- rare_carriers[i]
          idx2 <- rare_carriers[j]
          
          weight <- (geno[idx1] * geno[idx2]) / 4
          
          pairwise_sharing[idx1, idx2] <- pairwise_sharing[idx1, idx2] + weight
          pairwise_sharing[idx2, idx1] <- pairwise_sharing[idx2, idx1] + weight
        }
      }
    }
  }
  
  close(pb)
  
  pairwise_sharing <- pairwise_sharing / n_snp
  
  for (i in 1:length(unique_pops)) {
    for (j in 1:length(unique_pops)) {
      pop1_idx <- which(fam_subset$Population == unique_pops[i])
      pop2_idx <- which(fam_subset$Population == unique_pops[j])
      
      if (i == j) {
        if (length(pop1_idx) > 1) {
          within_pop <- pairwise_sharing[pop1_idx, pop1_idx]
          diag(within_pop) <- NA
          pop_sharing[i, j] <- mean(within_pop, na.rm = TRUE)
        }
      } else {
        pop_sharing[i, j] <- mean(pairwise_sharing[pop1_idx, pop2_idx])
      }
    }
  }
  
  return(list(
    pairwise = pairwise_sharing,
    population = pop_sharing,
    fam_data = fam_subset
  ))
}


target_population <- "Israel_C.AG" 

source_populations <- c(
  "Jordan_PPNB.AG", "Israel_PPNB.AG", "Jordan_Baja_Late_PPNB.AG",
  "Turkey_Marmara_Barcin_N.SG", "Turkey_Central_Catalhoyuk_N.SG", "Turkey_Southeast_Cayonu_PPN.SG",
  "Iran_SehGabi_C.AG", "Iran_HajjiFiruz_N.AG", "Iran_TepeAbdulHosein_N.SG", 
  "Iran_GanjDareh_N.AG", "Iran_Wezmeh_N.SG",
  "Iraq_Shanidar.AG", "Iraq_PPNA.AG", "Turkey_Southeast_Mardin_PPN.AG"
)

outgroup_populations <- c(
  "Russia_UstIshim_IUP.DG",
  "Russia_Kostenki14_UP.SG",
  "Turkey_Central_Pinarbasi_Epipaleolithic.AG",
  "Han.DG",
  "Papuan.DG",
  "Chukchi.DG",
  "Karitiana.DG",
  "Mbuti.DG"
)

all_pops <- c(target_population, source_populations)

rare_sharing_results <- calculate_rare_allele_sharing(
  bed_file = "AADR_v62/ancient_rare_variants",
  fam_data = fam_data,
  populations_of_interest = all_pops
)

pdf(file.path(results_dir, "out_population_rare_sharing_heatmap.pdf"), 
    width = 10, height = 8)

pheatmap(rare_sharing_results$population,
         main = "Population-level Rare Allele Sharing",
         color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
         display_numbers = FALSE,
         number_format = "%.4f",
         fontsize_number = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2")

dev.off()

#Focus on target population relationships
target_sharing <- rare_sharing_results$population[target_population, ]
target_sharing <- target_sharing[names(target_sharing) != target_population]
target_sharing <- sort(target_sharing, decreasing = TRUE)

# Bar plot of rare allele sharing with Israel_C
pdf(file.path(results_dir, "target_population_sharing.pdf"), 
    width = 12, height = 6)

ggplot(data.frame(Population = names(target_sharing), Sharing = target_sharing),
       aes(x = reorder(Population, Sharing), y = Sharing)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = paste("Rare Allele Sharing with", target_population),
       x = "Population",
       y = "Rare Allele Sharing Coefficient") +
  theme_minimal() +
  theme(text = element_text(size = 12))

dev.off()

perform_sharing_test <- function(sharing_matrix, fam_data, target_pop, source_pops) {
  
  results <- data.frame()
  
  target_idx <- which(fam_data$Population == target_pop)
  
  for (source in source_pops) {
    source_idx <- which(fam_data$Population == source)
    
    if (length(target_idx) > 0 && length(source_idx) > 0) {
      # Get sharing values
      sharing_vals <- as.vector(sharing_matrix[target_idx, source_idx])
      
      # Get background sharing (all other populations)
      other_idx <- which(!fam_data$Population %in% c(target_pop, source))
      background_vals <- as.vector(sharing_matrix[target_idx, other_idx])
      
      # Perform Wilcoxon test
      if (length(sharing_vals) > 0 && length(background_vals) > 0) {
        test_result <- wilcox.test(sharing_vals, background_vals, alternative = "greater")
        
        results <- rbind(results, data.frame(
          Target = target_pop,
          Source = source,
          Mean_sharing = mean(sharing_vals),
          Background_mean = mean(background_vals),
          Excess_sharing = mean(sharing_vals) - mean(background_vals),
          P_value = test_result$p.value,
          Significant = test_result$p.value < 0.05
        ))
      }
    }
  }
  
  results <- results %>% arrange(desc(Excess_sharing))
  return(results)
}

#Test for excess sharing
sharing_test_results <- perform_sharing_test(
  sharing_matrix = rare_sharing_results$pairwise,
  fam_data = rare_sharing_results$fam_data,
  target_pop = target_population,
  source_pops = source_populations
)

print(sharing_test_results)

rare_allele_excel <- list(
  "Summary" = data.frame(
    Analysis = "Rare Allele Sharing Analysis",
    Date = Sys.Date(),
    Target_population = target_population,
    Total_rare_SNPs = ncol(BEDMatrix("AADR_v62/ancient_rare_variants")),
    MAF_range = "0.01-0.05",
    N_populations = length(unique(rare_sharing_results$fam_data$Population)),
    N_individuals = nrow(rare_sharing_results$fam_data)
  ),
  
  "Population_Sharing_Matrix" = as.data.frame(rare_sharing_results$population),
  
  "Statistical_Tests" = sharing_test_results,
  
  "Target_Population_Sharing" = data.frame(
    Population = names(target_sharing),
    Sharing_coefficient = target_sharing
  )
)

write_xlsx(rare_allele_excel, 
           file.path(results_dir, "rare_allele_sharing_results.xlsx"))

cat("\nRare allele sharing analysis complete!\n")
cat("Results saved to:", file.path(results_dir), "\n\n")

#Read and plot admixture 

cv <- read.table("results_aDNA/cv_summary.tsv", header=TRUE, sep="\t")

pdf(file.path(results_dir, "Admixture_CV.pdf"), 
    width = 12, height = 6)
plot(cv$K, cv$CV_Error, type="b", pch=19, col="blue",
     xlab="K (number of clusters)", ylab="CV error",
     main="ADMIXTURE Cross-Validation Error")
grid()

dev.off()

library(tidyverse)
library(RColorBrewer)

cv_errors <- read.table("results_aDNA/cv_errors.txt", col.names = c("K", "CV_error"))

#Find optimal K
optimal_k <- cv_errors$K[which.min(cv_errors$CV_error)]
cat("Optimal K based on CV error:", optimal_k, "\n")

pdf(file.path(results_dir,"admixture_cv_plot.pdf"), width = 8, height = 6)
ggplot(cv_errors, aes(x = K, y = CV_error)) +
  geom_line() +
  geom_point(size = 3) +
  geom_point(data = cv_errors[which.min(cv_errors$CV_error),], 
             color = "red", size = 4) +
  labs(title = "ADMIXTURE Cross-Validation Error",
       x = "K (number of ancestral populations)",
       y = "Cross-validation error") +
  theme_minimal()
dev.off()


prefix <- "results_aDNA/ancient_pruned_final"        
ind_file <- "AADR_v62/v62.0_1240k_public.ind"        
out_pdf  <- "results_aDNA/admixture_barplots_by_population2.pdf"

fam <- read.table(paste0(prefix, ".fam"), header = FALSE, stringsAsFactors = FALSE)
iid <- fam$V2  

ind <- read.table(ind_file, header = FALSE, stringsAsFactors = FALSE)
colnames(ind) <- c("sample", "sex", "pop")

#Create a vector of populations ordered like fam/IID
iid2pop <- setNames(ind$pop, ind$sample)
pop_vec <- iid2pop[iid]

#Sanity check
if (any(is.na(pop_vec))) {
  missing <- iid[is.na(pop_vec)]
  stop(sprintf("These %d samples in .fam are missing in .ind, e.g.: %s",
               length(missing), paste(head(missing, 5), collapse = ", ")))
}

population_order_df <- tibble(
  Population = sort(unique(pop_vec))
) %>%
  mutate(
    group = case_when(
      Population == "Israel_C.AG" ~ "Israel_C",
      stringr::str_detect(Population, "Levant|Israel|Jordan|Natufian") ~ "Levant",
      stringr::str_detect(Population, "Anatolia|Turkey|Marmara|Barcin") ~ "Anatolia",
      stringr::str_detect(Population, "^Iran|Zagros") ~ "Zagros",
      stringr::str_detect(Population, "Iraq|Mesopotamia") ~ "Mesopotamia",
      TRUE ~ "Others"
    ),
    group_rank = match(group, c("Israel_C", "Levant", "Anatolia", "Zagros", "Mesopotamia", "Others")),
    group_rank = replace_na(group_rank, length(c("Israel_C", "Levant", "Anatolia", "Zagros", "Mesopotamia", "Others")))
  ) %>%
  arrange(group_rank, Population)

desired_population_order <- population_order_df$Population
display_label_lookup <- population_order_df %>%
  mutate(Display = stringr::str_replace(Population, "^Iran", "Zagros")) %>%
  select(Population, Display)

get_cols <- function(k) {
  if (k <= 9) {
    brewer.pal(max(k, 3), "Set1")[1:k]
  } else {
    colorRampPalette(brewer.pal(9, "Set1"))(k)
  }
}

cols_for_k <- function(k) {
  if (k <= 9) {
    brewer.pal( max(3, k), "Set1" )[1:k]
  } else {
    colorRampPalette(brewer.pal(9, "Set1"))(k)
  }
}

order_populations <- function(pop_means_matrix, method = c("by_firstK_desc","alphabetical")) {
  method <- match.arg(method)
  if (method == "by_firstK_desc") {
    ord <- order(pop_means_matrix[,1], decreasing = TRUE)
  } else {
    ord <- order(rownames(pop_means_matrix))
  }
  pop_means_matrix[ord, , drop = FALSE]
}

pdf(out_pdf, width = 11, height = 7)  # landscape-ish page
for (k in 2:12) {
  qfile <- sprintf("%s.%d.Q", prefix, k)
  if (!file.exists(qfile)) {
    # blank page placeholder if missing
    plot.new(); title(main = sprintf("ADMIXTURE by population (K = %d) — missing %s", k, basename(qfile)))
    next
  }
  
  Q <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(Q) != length(pop_vec)) {
    stop(sprintf("Row mismatch at K=%d: Q rows=%d, fam IIDs=%d", k, nrow(Q), length(pop_vec)))
  }
  
  # mean ancestry per population (rows = pops, cols = K components)
  pops <- unique(pop_vec)
  pop_means <- t(sapply(pops, function(p) {
    rows <- which(pop_vec == p)
    colMeans(as.matrix(Q[rows, , drop = FALSE]))
  }))
  rownames(pop_means) <- pops
  colnames(pop_means) <- paste0("K", seq_len(k))
  
  # order pops (change method to "alphabetical" if preferred)
  pop_means <- order_populations(pop_means, method = "by_firstK_desc")
  
  pop_levels <- desired_population_order[desired_population_order %in% rownames(pop_means)]
  pop_means <- pop_means[pop_levels, , drop = FALSE]
  
  display_levels <- display_label_lookup$Display[match(pop_levels, display_label_lookup$Population)]
  display_levels[is.na(display_levels)] <- pop_levels[is.na(display_levels)]
  
  # long format for ggplot
  df <- as.data.frame(pop_means) %>%
    mutate(
      Population = factor(pop_levels, levels = pop_levels),
      DisplayPopulation = factor(display_levels, levels = rev(display_levels))
    )
  
  df_long <- df |>
    pivot_longer(cols = starts_with("K"), names_to = "Cluster", values_to = "Ancestry")
  
  axis_levels <- levels(df$DisplayPopulation)
  axis_labels <- axis_levels
  israel_label <- display_label_lookup$Display[display_label_lookup$Population == "Israel_C.AG"]
  israel_label <- if (length(israel_label)) israel_label else "Israel_C.AG"
  if (israel_label %in% axis_levels) {
    axis_labels[axis_levels == israel_label] <- paste0("<b>", israel_label, "</b>")
  }
  names(axis_labels) <- axis_levels
  
  p <- ggplot(df_long, aes(x = DisplayPopulation, y = Ancestry, fill = Cluster)) +
    geom_col(width = 0.8) +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(limits = c(0, 1.001), expand = expansion(mult = c(0, 0))) +
    coord_flip(ylim = c(0, 1)) +
    scale_fill_manual(values = cols_for_k(k), guide = guide_legend(nrow = 1)) +
    labs(title = sprintf("ADMIXTURE by population (K = %d)", k),
         x = NULL, y = "Mean ancestry proportion") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      axis.text.y = element_markdown(size = 9),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank()
    )
  
  print(p)
}
dev.off()
