#both samples from iraq with anitolia
library(admixtools)

# CONFIGURATION 
aadr_path <- "path/to/aadr/folder"
prefix <- file.path(aadr_path, "v62.0_1240k_public")

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

iraq_pops <- c(
  "Iraq_Shanidar.AG",
  "Iraq_PPNA.AG"
)

anitolia_pops <- c(
  "Turkey_Marmara_Barcin_N.SG",   
  "Turkey_Central_AsikliHoyuk_PPN_o.SG", 
  "Turkey_Central_Catalhoyuk_N.SG",
  "Turkey_Southeast_Cayonu_PPN.SG",
  "Turkey_Southeast_Mardin_PPN.AG",
  "Turkey_Central_Boncuklu_PPN.WGC.SG",
  "Turkey_Central_Boncuklu_PPN.AG"
)

# GENERATE VALID COMBINATIONS (only sizes 2,3,4)
valid_combos <- function(all_pops, sizes) {
  combos <- list()
  idx <- 1
  
  for(k in sizes){
    cmb <- combn(all_pops, k, simplify = FALSE)
    
    for(c in cmb){
      # keep only mixed Iran-Iraq sets
      if(any(c %in% iraq_pops) && any(c %in% anitolia_pops)){
        combos[[idx]] <- c
        idx <- idx + 1
      }
    }
  }
  return(combos)
}
all_pops <- c(iraq_pops, anitolia_pops)
combo_list <- valid_combos(all_pops, sizes = c(2,3,4))

length(combo_list)   # number of valid sets
combo_list            # inspect sets

#Loop

results <- list()

for(i in seq_along(combo_list)){
  
  left_set <- combo_list[[i]]
  name <- paste(left_set, collapse="__")
  
  cat("Running combination:", paste(left_set, collapse=", "), "\n")
  
  pops_for_f2 <- c(left_set, outgroups)
  
  # Capture the printed output from f2_from_geno 
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
  formatted_text <- paste0(
    "For \"", paste(left_set, collapse="\" + \""), "\"\n\n",
    "F2 blocks\n",
    snps_total, "\n",
    snps_filtered, "\n\n",
    "QpWave\n",
    capture.output(print(qp$rankdrop)),
    collapse = "\n"
  )
  
  results[[name]] <- list(
    snps_total    = snps_total,
    snps_filtered = snps_filtered,
    rankdrop      = qp$rankdrop,
    formatted     = formatted_text
  )
  
  
  print(qp$rankdrop)
}


results

print_result <- function(results_list) {
  for(name in names(results_list)) {
    
    entry <- results_list[[name]]
    
    cat("For:", gsub("__", " + ", name), "\n")
    
    cat("F2 blocks\n")
    if(!is.null(entry$snps_total))    cat(entry$snps_total, "\n")
    if(!is.null(entry$snps_filtered)) cat(entry$snps_filtered, "\n\n")
    
    cat("qpWave\n")
    print(entry$rankdrop)
    cat("\n")
  }
}

print_result(results)



# writing the results in docx file 

library(flextable)
library(officer)

write_results_to_docx <- function(results_list, file_name = "qpwave_results.docx") {
  
  library(officer)
  library(flextable)
  
  # bold formatting
  bold_style <- fp_text(bold = TRUE)
  
  doc <- read_docx()
  
  for(name in names(results_list)) {
    entry <- results_list[[name]]
    
    # H2 TITLE
    h2_title <- paste("For:", gsub("__", " + ", name))
    doc <- body_add_par(doc, h2_title, style = "heading 2")
    doc <- body_add_par(doc, "")   # spacing
    
    # F2 BLOCKS (BOLD HEADING)
    doc <- body_add_fpar(doc, fpar(ftext("F2 blocks", prop = bold_style)))
    
    if(!is.null(entry$snps_total))
      doc <- body_add_par(doc, entry$snps_total)
    
    if(!is.null(entry$snps_filtered))
      doc <- body_add_par(doc, entry$snps_filtered)
    
    doc <- body_add_par(doc, "")
    
    # QPWAVE (BOLD HEADING)
    doc <- body_add_fpar(doc, fpar(ftext("QpWave", prop = bold_style)))
    
    # Convert tibble to flextable (makes a nice proper Word table)
    ft <- flextable(entry$rankdrop)
    ft <- autofit(ft)
    doc <- body_add_flextable(doc, ft)
    
    doc <- body_add_par(doc, "")
  }
  
  print(doc, target = file_name)
  cat("Saved:", file_name, "\n")
}


write_results_to_docx(results, "qpwave_results.docx")




