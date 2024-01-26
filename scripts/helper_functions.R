# Made by: Dennis Scheper
# Contact: d.j.scheper@umcg.nl / d.j.scheper@st.hanze.nl / djscheper@gmail.com

# Libraries used
library(dplyr)
library(cluster)
library(ReactomePA)
library(circlize)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggupset)
library(purrr)
library(ComplexHeatmap)
library(enrichR)
library(glue)

#' Reorder meta and count datasets
#'
#' This function will reorder a counts dataset based on the order of samples in a meta dataset
#' @param meta_data - a meta dataset with sample identifiers
#' @param count_data - a counts dataset which will be reordered
#' @return returns reordered count dataset
#' @examples
#' meta <- read.csv('path/to/meta.csv') 
#' counts <- read.csv('path/to/counts.csv')
#' ordered_counts <- reordered(meta, counts)
reordered <- function(meta_data, count_data) {
  all.in <- all(meta_data$sample_identifier %in% colnames(count_data))
  ord <- all(meta_data$sample_identifier == colnames(count_data))
  if (ord == F) {
    reord <- match(meta_data$sample_identifier, colnames(count_data))
    new.counts <- count_data[, reord]
    return(new.counts)
  } 
  
  return(count_data)
}

#' Perform PCA
#'
#' This function performs PCA on any counts dataset and returns the results and eigenvalues
#' @param df - a normalized/standardized (counts) dataset
#' @return returns a results list with PCA components (PCs) and eigenvalues
#' @examples
#' counts <- read.csv('path/to/counts.csv')
#' pca_res <- perform_pca(counts)
perform_pca <- function(df, norm=F, center=T, scale=T) {
  
  # Transform and assign to matrix
  matrix_pca <- t(df) %>% as.matrix()
  
  if (norm) {
    matrix_pca <- log2(matrix_pca + 0.01)
  }
  
  # Perform PCA
  res_pca <- prcomp(matrix_pca, center = center, scale. = scale)
  
  # How many components do we need (depending on the number of columns and rows)
  n <- number_pca_components(df)
  
  eigvals <- res_pca$sdev^2 %>%
    data.frame(PC = factor(1:n), var = .) %>%
    mutate(pva = var/sum(var) * 100) %>% # variance per PC
    mutate(p_cum = cumsum(pva))          # cumulative percentage per additional PC
  
  result <- list(pca = res_pca, eigenvalues = eigvals)
  
  return(result)
}

#' Determines the amount of PCs needed
#'
#' Determines the amount of PCs needed for eigenvalues dataframe (no need to use this function yourself)
#' @param df - a normalized/standardized (counts) dataset
#' @param count_data - a counts dataset which will be reordered
#' @return returns the number of PCs
#' @examples
#' Don't use this function outside of 'perform_pca()'

number_pca_components <- function(df) {
  
  # Calculate row and column length
  row_length <- nrow(df)
  col_length <- ncol(df)
  
  # The number of meaningful PCs is the smaller of the two
  n_comp <- min(row_length, col_length)
  
  return(n_comp)
}

#' Function to decide the amount of PCs
#'
#' This function determines how many PCs we need to explain a percentage of variance
#' @param df - a dataframe containing eigenvalues and variance per principal component
#' @param pva - a number in percentage representing the amount of explained variance
#' @return returns the sufficient number of components explaining the number represented in pva
#' @examples
#' pca_res <- perform_pca(counts_dataset)
#' number_components <- explain_variance(pca_res$eigenvalues$pva, 80)
explain_variance <- function(df, pva) {
  df <- as.data.frame(df)
  num_components <- which(df >= pva)[1]
  return(num_components)
}

#' Function to perform DE analysis with
#'
#' This function preprocesses count data, 
#' performs DE analysis via DESeq2 for all levels in a design variable,
#' and handles the processing of DEGs per comparison.
#' @param countData - RNA-Seq count dataframe
#' @param metaData - meta dataframe
#' @param design.element - which variable needs to be assessed
#' @param low.exp - threshold for discarding low-expressed genes (standard = 10)
#' @param smallest.group.size - how many samples do need to have at least less than low.exp counts for a gene to be discarded
#' @param versus_healthy - DE analysis against healthy controls (standard = FALSE)
#' @param hc_name -  starting prefix of healthy controls (standard = 'hc')
#' @param pvalue_thres - threshold p-value (standard = 0.05)
#' @param fold_thres - threshold fold change (standard = 1.2)
#' @param batch - perform batch correction (standard = TRUE)
#' @return a list with DESeq2's results ($dds), DE analysis results per comparison ($results), and levels that were compared ($lvls)
#' @examples
#' meta_data <- read.csv('path/to/meta.csv')
#' count_data <- read.csv('path/to/counts.csv)
#' de_analysis <- de(count_data, meta_data, "sepsis_severity")
de_analysis <- function(countData, metaData, design.element, low.exp = 10, 
                        smallest.group.size = 10, versus_healthy = F, hc_name = "hc", 
                        pvalue_thres = 0.05, fold_thres = 1.2, batch = T) {
  
  # Check if design parameters are in metaData
  if (!all(design.element %in% colnames(metaData)) && design.element != '1') {
    stop(glue::glue('variable not found in colData. \n Your input: {design.element}'))
  }
  
  # Insert batch correct variable if true
  if (batch) {
    design <- as.formula(paste("~ sequencing_month_year + ", paste(design.element), collapse = " + "))
  } else {
    design <- as.formula(paste("~", paste(design.element), collapse = " + "))
  }
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = design)
  
  # Remove low-expressed genes
  if (low.exp > 0) {
    keep <- rowSums(counts(dds) >= low.exp) >= smallest.group.size
    dds <- dds[keep, ]
  }
  
  dds <- DESeq(dds)
  
  lvls <- levels(colData(dds)[[design.element]]) # Get all levels
  comparisons <- combn(lvls, 2, simplify = F) # Make all combinations
  
  if (versus_healthy == TRUE) {
    # Only keep comparisons against healthy controls
    comparisons <- keep(comparisons, ~any(.x == hc_name))
  }
  
  # Retrieve results for all comparisons
  all_results <- list()
  for (i in 1:length(comparisons)) {
    res <- DESeq2::results(dds, contrast = 
                             c(design.element, comparisons[[i]][2], comparisons[[i]][1]),
                           independentFiltering=FALSE) %>%
      as.data.frame() %>%
      mutate(absolute_change = abs(log2FoldChange)) %>%
      # Remove the log2 scale
      mutate(fold_change = case_when(
        log2FoldChange > 0 ~ 2^absolute_change,
        log2FoldChange < 0 ~ -2^absolute_change,
        TRUE ~ 0)) %>%
      mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
      # A DE is when: padj is lower than pvalue_thres and 
      # fold_change is lower than -fold_thres or higher than fold_thres
      mutate(de = ifelse(padj <= pvalue_thres & 
                           (fold_change >= fold_thres | 
                              fold_change <= -fold_thres),
                         "DE", "NO")) %>%
      # Assign which direction a DE goes
      mutate(direction = case_when(
        fold_change >= 0 & de == "DE" ~ "UP",
        fold_change <= 0 & de == "DE" ~ "DOWN",
        padj > 0.05 ~ 'NOT_SIG'
      ))
    
    # Outline which comparison was made  
    name <- paste0(comparisons[[i]][2], "_against_", comparisons[[i]][1])
    all_results[[name]] <- res
  }
  
  return(list(dds = dds, results = all_results, levels = lvls))
}

#' Summary function for DE analysis
#'
#' This function summarizes findings for DE analysis by extracting DEGs per comparison
#' @param df - a dataframe with results from DE analysis
#' @return a dataframe with the amount of DEGs in every direction
#' @examples
#' meta_data <- read.csv('path/to/meta.csv')
#' count_data <- read.csv('path/to/counts.csv)
#' de_analysis <- de(count_data, meta_data, "sepsis_severity")
#' print(summary_de(de_analysis$results$level1_against_level2))
summary_de <- function(df) {
  res <- data.frame(
    all_degs = df %>% filter(de == "DE") %>% nrow(),
    up_degs = df %>% filter(de == "DE" & direction == "UP") %>% nrow(),
    down_degs = df %>% filter(de == "DE" & direction == "DOWN") %>% nrow()
  )
  
  return(res)
}

#' Function normalize count data
#'
#' This function normalizes count data the same as DESeq2 internal functionality does. This is intended for downstream analysis.
#' @param counts - RNA-Seq count data
#' @param pva - a dataframe with meta data
#' @return a normalized count dataframe
#' @examples
#' meta_data <- read.csv('path/to/meta.csv')
#' count_data <- read.csv('path/to/counts.csv)
#' normalized_counts <- normalize(count_data, meta_data)
normalize <- function(counts, meta) {
  normalized_counts <- counts %>% as.matrix() %>%
    DESeq2::varianceStabilizingTransformation() %>%
    sva::ComBat(batch = meta$sequencing_month_year)
  return(normalized_counts)
}

#' Pathway analysis through MSigdb database
#'
#' Given a list of genes as input, the function will search the MSigdb databases in search of significant pathways
#' @param genes - a list of gene names
#' @param pvalue - the p-value cutoff (standard = 0.05)
#' @param dbs - a vector of databases to search in
#' @return returns the sufficient number of components explaining the number represented in pva
#' @examples
#' pca_res <- perform_pca(counts_dataset)
#' number_components <- explain_variance(pca_res$eigenvalues$pva, 80)
go_pathway_analysis <- function(genes, pvalue = 0.05, 
                                dbs = c("MSigDB_Hallmark_2020",
                                        "GO_Molecular_Function_2018", 
                                        "GO_Cellular_Component_2018")) {
  # Remove any NAs
  genes <- na.omit(genes) %>% as.character()
  
  res <- enrichr(genes = genes, databases = dbs)
  
  # Filter on pvalue for every database used
  res <- lapply(res, function(x) {
    filter(x, Adjusted.P.value <= pvalue)
  })
  
  return(res)
}

#' Pathway analysis through Reactome database
#'
#' Given a list of genes as input, the function will search the Reactome database in search of significant pathways
#' @param genes - a list of gene names
#' @param pvalue - the p-value cutoff (standard = 0.05)
#' @param universe - a list of background gene names
#' @return returns a list of up- and downregulated pathways and information about them
#' @examples
#' gene_list <- c('ARG1', 'ARG2', 'THEM4')
#' reactome_res <- reactome_pathway_analysis(gene_list)
reactome_pathway_analysis <- function(genes, pvalue = 0.05, universe = NULL) {
  
  # Remove any NAs if necessary
  genes <- na.omit(genes)
  universe <- na.omit(universe)
  
  # Arguments for enrichPathway function
  arguments <- list(gene = genes, pvalueCutoff = pvalue, readable = TRUE, minGSSize = 2)
  
  # Add universe to arguments list if not null
  if(!is.null(universe)) {
    arguments$universe <- as.character(universe$entrez_id)
  }
  
  # Call the enrichPathway function like this so the universe parameter is optional
  # without the need for an if-else statement; only take the results list via '@results'
  res <- do.call("enrichPathway", arguments)@result %>%
    filter(p.adjust <= 0.05)
  
  return(res)
}

#' Function to perform pathway analysis with
#'
#' This function performs pathway analysis for all comparisons and directions (if necessary) for Reactome and MgsigDB
#' @param comparisons - a list of comparisons of DEGs (for example, levelone_against_leveltwo)
#' @param split - if true, split between up and downregulated pathways (standard = FALSE)
#' @param universe - background gene names
#' @return returns a list of pathway analysis results from Reactome and MgSigDB databases
#' @examples
#' res <- pathway_analysis(results_from_de_analysis, split = TRUE)
pathway_analysis <- function(comparisons, split = F, universe = NULL) {
  
  # Retrieve Entrez IDs for background genes
  if (!is.null(universe)) {
    entrez_ids_universe <- AnnotationDbi::select(org.Hs.eg.db, keys=universe$gene, 
                                                 column=c("ENTREZID", "SYMBOL"), keytype="SYMBOL") %>%
      rename(entrez_id = ENTREZID) %>%
      dplyr::select(-SYMBOL)
  }
  
  # For every comparison, retrieve gene names and extract Entrez IDs
  degs_list <- comparisons %>%
    map(~ .x %>% 
          mutate(entrez_id = AnnotationDbi::select(org.Hs.eg.db, keys = gene, "ENTREZID", keytype="SYMBOL")) %>%
          tidyr::unnest(c(entrez_id)) %>% 
          dplyr::select(-SYMBOL) %>%
          rename(entrez_id = ENTREZID))
  
  if (split) {
    # Split between down- and upregulated for every comparison (for a list with sublists)
    degs_list <- degs_list %>%
      map(~ .x %>% split(.x$direction))
    if (!is.null(universe)) {
      # Reactome for every sublist in the sublists in the deg list (for example, deg_list$level1_against_level2$DOWN) with universe
      reactome_res <- degs_list %>%
        map(~ map(., ~ reactome_pathway_analysis(.x$entrez_id, universe = entrez_ids_universe)))
    } else {
      # Reactome for every sublist in the sublists in the deg list (for example, deg_list$level1_against_level2$DOWN) without universe
      reactome_res <- degs_list %>%
        map(~ map(., ~ reactome_pathway_analysis(.x$entrez_id)))
    }
    go_res <- degs_list %>%
      # MgSigDB for every sublist in the sublists in the deg list (for example, deg_list$level1_against_level2$DOWN)
      map(~ map(., ~ go_pathway_analysis(.x$gene)))
    
  } 
  # The same as above but without the split
  else {
    if (!is.null(universe)) {
      reactome_res <- degs_list %>%
        map(~ reactome_pathway_analysis(.x$entrez_id, universe = entrez_ids_universe))
    } else {
      reactome_res <- degs_list %>%
        map(~ reactome_pathway_analysis(.x$entrez_id))
    }
    
    go_res <- degs_list %>%
      map(~ go_pathway_analysis(.x$gene))
  }
  
  return(list(reactome_res = reactome_res, go_res = go_res))
}

#' Function to plot dotplots based on MgSigDB results
#'
#' This function plots a dotplot regarding MgSigDB results
#' @param comparison_list - a dataframe containing comparisons of the DESeq2 DE analysis
#' @return returns the result of the plotting
#' @examples
#' pathway_analysis(comparisons_list)
plot_go <- function(comparison_list) {
  
  comparison_list <- comparison_list %>%
    map(~ map(., ~  map(., ~ .x %>%
                          separate_wider_delim(Overlap, delim = '/', 
                                               names = c('background_gene_count', 'total')) %>%
                          mutate(Ratio = as.numeric(background_gene_count) / as.numeric(total))))) %>%
    map(~ map(., ~ bind_rows(.x, .id = "database"))) %>%
    map(~ bind_rows(.x, .id = 'direction')) %>%
    bind_rows(.id = 'comparison') %>%
    mutate(Term = ifelse(database != 'MSigDB_Hallmark_2020', 
                         sapply(Term, extract_words), Term)) %>%
    
    mutate(comparison = str_replace_all(comparison, 
                                        "([^\\s_]+)_against_([^\\s_]+)", # works for any character of any length
                                        "\\1 vs. \\2"))
  
  res <- lapply(unique(comparison_list$database), function(x) {
    comparison_list %>% filter(x == database) %>%
      ggplot(aes(y = Term, 
                 x = as.factor(direction), 
                 fill = Ratio, 
                 size = -log10(Adjusted.P.value), 
                 color = direction)) +
      geom_point(pch = 16) +
      facet_wrap(~ comparison) + 
      labs(size = "-Log10PAdjusted", color = "Direction", 
           y = "Pathway", x = "Direction") +
      scale_color_manual(values = c("UP" = "orange", "DOWN" = "magenta"),
                         labels = c("UP" = "Upregulated", "DOWN" = "Downregulated")) +
      scale_x_discrete(labels = c("UP" = "Upregulated", "DOWN" = "Downregulated")) +
      theme(axis.text.y=element_text(size=7)) + 
      guides(fill = FALSE, # Remove RATIO legend
             size = guide_legend(override.aes = list(pch = 1))) +
      ggtitle(glue("Comparison between clusters ({x})"))
  })
  
  return(res)
}

extract_words <- function(description) {
  number_words <- str_count(description, "\\S+")
  words_extract <- min(number_words, 3)
  word(description, 1, words_extract)
}

#' Function to plot dotplots based on Reactome results
#'
#' This function plots a dotplot regarding Reactome results
#' @param comparison_list - a dataframe containing comparisons of the DESeq2 DE analysis
#' @return returns the result of the plotting
#' @examples
#' pathway_analysis(comparisons_list)
plot_reactome <- function(comparison_list) {
  
  comparison_list <- comparison_list %>%
    map(~ map(., ~ .x %>% 
                separate_wider_delim(BgRatio, delim = '/', 
                                     names = c('background_gene_count', 'total')) %>% # split BgRatio into two columns
                mutate(Ratio = Count / as.numeric(background_gene_count)) %>%
                mutate(Description = sapply(Description, extract_words)))) %>%
    map(~ bind_rows(.x, .id = 'direction')) %>%
    bind_rows(.id = 'comparison') %>%
    mutate(comparison = str_replace_all(comparison, 
                                        "([^\\s_]+)_against_([^\\s_]+)", # works for any character of any length
                                        "\\1 vs. \\2"))
  
  res <- ggplot(comparison_list, 
                aes(y = Description, 
                    x = as.factor(direction), 
                    fill = Ratio, 
                    size = -log10(p.adjust), 
                    color = direction)) +
    facet_wrap(~ comparison) +
    geom_point(pch = 16) +
    labs(size = "-Log10PAdjusted", color = "Direction", 
         y = "Pathway", x = "Direction") +
    scale_color_manual(values = c("UP" = "orange", "DOWN" = "magenta"),
                       labels = c("UP" = "Upregulated", "DOWN" = "Downregulated")) +
    scale_x_discrete(labels = c("UP" = "Upregulated", "DOWN" = "Downregulated")) +
    theme(axis.text.y=element_text(size=7, vjust = 0.8)) + 
    guides(fill = FALSE, # Remove RATIO legend
           size = guide_legend(override.aes = list(pch = 1))) +
    ggtitle("Comparison between clusters")
  
  return(res)
}

#' Function to generate a heatmap
#'
#' This function generates a heatmap based on gene expression per cluster/endotype
#' @param meta - meta dataset
#' @param counts - RNA-Seq counts dataset
#' @param type - indicating whether a cohort belongs to ER or ICU (due to different SOFA score variables)
#' @return returns the sufficient number of components explaining the number represented in pva
#' @examples
#' meta_data <- read.csv('path/to/meta.csv')
#' count_data <- read.csv('path/to/counts.csv')
#' generate_heatmap(meta_data, count_data, 'er')

generate_heatmap <- function(meta, count, type) {
  
  # Check whether names are already adjusted
  if ("pam.clusters" %in% names(meta)) {
    meta <- meta %>%
      dplyr::rename(cluster = pam.clusters)
  }
  
  # Amount of clusters
  num_cluster <- unique(meta$cluster)
  
  # Scale gene count
  count <- scale(count, center = TRUE, scale = TRUE)
  
  # Scale color usage
  breaks <- c(min(count), 0, max(count))
  colors <- c("blue", "black", "yellow")
  used_colors <- colorRamp2(breaks, colors)
  
  # Generate different column annotation based on type parameter
  if (type == "er") {
    sofa <- meta$worst_within_72_sofa
    column_annotation <- HeatmapAnnotation(
      Severity = meta$sepsis_severity, 
      Mortality = meta$mortality, 
      SOFA = anno_barplot(sofa, which = "column", border = T),
      Cluster = anno_block(gp = gpar(fill = c('purple', 'green'),
                                     col="white"),
                           labels = num_cluster),
      col = list(Severity = c("High" = "red",
                              "Intermediate" = "blue",
                              "Low" = "green"),
                 Mortality = c("deceased" = "black",
                               "survived" = "yellow",
                               "unknown" = "grey"))
    )
    
  } else if (type == "icu") {
    sofa <- meta$icu_sofa
    column_annotation <- HeatmapAnnotation(
      Severity = meta$sepsis_severity,
      SOFA = anno_barplot(sofa, which = "column", border = T),
      Cluster = anno_block(gp = gpar(fill = c('orange', 'blue'),
                                     col="white"),
                           labels = num_cluster),
      col = list(Severity = c("High" = "red",
                              "Intermediate" = "blue",
                              "Low" = "green"))
    )
  } else {
    stop(glue::glue('
      Please make sure you use either "icu" or "er" as your type parameter 
      \n Your input: {type}'))
  }
  
  # Generate the heatmap
  heatmaps_cluster <- count %>% t() %>%
    Heatmap(column_split = meta$cluster,
            show_column_names = FALSE, 
            show_row_names = FALSE, 
            name = "Z-score",
            top_annotation = column_annotation,
            column_title = NULL, col = used_colors,
            show_row_dend = FALSE, show_column_dend = FALSE)
  
  draw(heatmaps_cluster)
}

