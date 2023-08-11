stop()
rm(list=ls())

source("common.R")
source("common-RNAseq.R")
source("common-gsea.R")

load("data/de_bay.RData")
load("data/de_mn.RData")
load("data/de_rupture.RData")

load("data/de_crasta.RData")
load("data/de_santaguida_cycling.RData")
load("data/de_santaguida_arrested.RData")

# Hypergeometric test: 
# If we have a bucket with some proportion p of red balls and 1-p green balls
# and we draw n balls, did we draw more red balls than expected by chance alone?

# The Mps1i DEGs are the drawn balls (or He, et al or Santaguida, et al), and the bucket contain either MN DEGs and non-DEGs
# rupture DEGs and non-DEGs.
hg_test <- function(drawn, categories, total, fc=0) {
  filtered_drawn <- list() # this is the Mps1i DEG list?
  filtered_categories <- list()
  
  # This loops through the drawn list and pulls out the gene names
  # filtering by DEG and abs(FC) >= the fc parameter above
  
  # Filtering out genes with no gene name
  # Using gene name rather than ensemble ID because multiple IDs 
  # can reference the same gene
  for(name in names(drawn)) {
    filtered_drawn[[name]] <- drawn[[name]]$res_shrunk_df_reduced %>%
      filter(padj <= 0.05 & ext_gene != "" & abs(log2FoldChange) >= fc) %>%
      distinct(ext_gene, .keep_all=T) %>%
      pull(ext_gene)
    print(sprintf("%d %s DEGs", length(filtered_drawn[[name]]), name))
  }
  
  # This loops through the categories list and pulls out the gene names
  # filtering by DEG and abs(FC) >= the fc parameter above
  for(name in names(categories)) {
    filtered_categories[[name]] <- categories[[name]]$res_shrunk_df_reduced %>%
      filter(padj <= 0.05 & ext_gene != "" & abs(log2FoldChange) >= fc) %>%
      distinct(ext_gene, .keep_all=T) %>%
      pull(ext_gene)
    print(sprintf("%d %s DEGs", length(filtered_categories[[name]]), name))
  }
  
  er <- hyperEnrichment(filtered_drawn, filtered_categories, ntotal=total)
  er
}

# From CBMRtools
if (!require(CBMRtools, quietly = TRUE)) {
  library(devtools)
  install_github("tetomonti/CBMRtools/CBMRtools")
  require(CBMRtools)
}

drawn <- list(
  "Mps1i"=de_bay,
  "He, et al."=de_crasta,
  "Santaguida, et al."=de_santaguida_arrested
)
categories <- list(
  "MN+/-"=de_mn,
  "Rupture+/-"=de_rupture
)

# The total number of annotated genes
back_pop <- de_bay$res_shrunk_df_reduced %>% 
  filter(ext_gene != "") %>% 
  distinct(ext_gene, .keep_all=T) %>% nrow()

# No FC filter
result <- hg_test(drawn, categories, back_pop)
tibble(result)

# FC >= 1.5 filter
result <- hg_test(drawn, categories, back_pop, fc=log2(1.5))
tibble(result)

# FC >= 2.0 filter
# There are simply an insufficient number of DEGs for this to be meaningful, I think
result <- hg_test(drawn, categories, back_pop, fc=log2(2.0))
tibble(result)



gsea_test <- function(drawn, categories) {
  filtered_categories <- list()
  
  results <- list()
  
  for(name in names(categories)) {
    filtered_categories[[name]] <- categories[[name]]$res_shrunk_df_reduced %>%
      filter(padj <= 0.05 & ext_gene != "") %>%
      distinct(ext_gene, .keep_all=T) %>%
      pull(ext_gene)
    print(sprintf("%d %s DEGs", length(filtered_categories[[name]]), name))  
  }
  
  for(name in names(drawn)) {
    degs <- drawn[[name]]$res_shrunk_df_reduced %>%
      filter(padj <= 0.05 & ext_gene != "") %>%
      distinct(ext_gene, .keep_all=T) %>%
      select(ext_gene, log2FoldChange) %>%
      arrange(desc(log2FoldChange)) %>%
      deframe()
    print(sprintf("%d %s DEGs", length(degs), name)) 
    
    results[[name]] <- list(
      "table"=fgsea(
        pathways=filtered_categories,
        minSize=5,
        maxSize=2500,
        stats=degs,
        nproc=1,
        nPermSimple=10000
      ) %>%
        as_tibble() %>%
        arrange(padj),
      "plots"=list()
    )
    
    for(i in names(filtered_categories)) {
      results[[name]]$plots[[i]] <- plotEnrichment(
        filtered_categories[[i]], 
        degs
      ) + labs(subtitle=i, caption=format_p_value(results[[name]]$table %>% filter(pathway == i) %>% pull(padj), type=""))
    }
    
    results[[name]]$plot <- wrap_plots(results[[name]]$plots, nrow=2) + plot_annotation(title=name)
  }
  
  results
}

results <- gsea_test(drawn, categories)
results$Mps1i$table
results$Mps1i$plot

results$`He, et al.`$table
results$`He, et al.`$plot

results$`Santaguida, et al.`$table
results$`Santaguida, et al.`$plot
