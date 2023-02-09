# ----- D
# This code performs differential gene expression analysis using the limma
# package in R. The first step is to create the design matrix (`int_design`) 
# which specifies the factors and interactions that will be used in the model.
# The `model.matrix` function takes the `sample_lookup` data frame and modifies
# the "protocol" and "passage" factors to be re-leveled with the `relevel` 
# function such that the reference_passage and "p" (protocol) are the reference
# levels, respectively.
#
# Next, the code filters the wide expression matrix (`expr_wide`) to only keep 
# the rows that correspond to the samples in the design matrix (`int_keep`). 
# This filtered expression matrix is then transformed into a voom-normalized 
# object (`int_voomed`) using the `voom` function. The linear regression model 
# is fit to the voom-normalized expression data (`int_fit`), and then an empirical
# Bayesian adjustment is applied to the model (`int_fit_ebayes`) using the 
# `eBayes` function.
# 
# Finally, the top differentially expressed genes are selected using the 
# `topTable` function, which returns a data frame of the genes with adjusted 
# P-values below a specified threshold (in this case, not specified). The 
# resulting table is converted to a tibble data structure and joined with gene 
# annotation information from `gene_lookup` using the `left_join` function. The 
# resulting data frame is then modified to add a column indicating the row 
# number (`hit_num`).jh
do_medium_analysis <- function(filename, reference_passage) {
  int_design <- model.matrix(~ passage * protocol,
                             data = sample_lookup |>
                               mutate(
                                 protocol = relevel(factor(protocol), "p"),
                                 passage = relevel(factor(passage), reference_passage)
                               )
  )
  
  int_keep <- filterByExpr(expr_wide |> select(-gene_id), int_design)
  
  
  int_voomed <- voom(
    expr_wide[int_keep, sample_lookup$sample],
    int_design,
    plot = TRUE
  )
  
  int_fit <- lmFit(int_voomed, int_design)
  
  
  int_fit_ebayes <- eBayes(int_fit)
  
  int_limma_res <- topTable(int_fit_ebayes,
                            coef = "protocols",
                            genelist = expr_wide |> filter(int_keep) |> select(gene_id),
                            number = 20000000,
                            confint = TRUE
  ) |>
    as_tibble() |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(hit_num = row_number())
  
  
  # Gene ontologies ---------------
  
  int_enricher_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_gene
  )
  
  int_enricher_res <- int_enricher_obj@result |> as_tibble()
  
  go_lookup <- go2term(int_enricher_res$ID) |>
    as_tibble() |>
    left_join(go2ont(int_enricher_res$ID), by = "go_id") |>
    mutate(Ontology = case_when(
      Ontology == "BP" ~ "Biological Process",
      Ontology == "MF" ~ "Molecular Function",
      Ontology == "CC" ~ "Cellular Component"
    ))
  
  int_enricher_res_tbl <- int_enricher_res |>
    left_join(go_lookup, by = c("ID" = "go_id")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(Term, Ontology, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  
  data_long <- expr_long_cpm |>
    left_join(sample_lookup, by = "sample") |>
    left_join(gene_lookup, by = "gene_id") |>
    mutate(gene_name = gene_id)
  
  
  ## KEGG
  kegg_obj <- enricher(
    gene = int_limma_res |> filter(adj.P.Val < 0.05) |> pull(gene_id),
    universe = int_limma_res |> pull(gene_id),
    TERM2GENE = term_to_brite
  )
  # enricher: A universal enrichment analyzer
  # genes are signficant genes, universe is all genes. Later used for Fisher's exact test?
  
  kegg_res_tbl <- kegg_obj@result |> 
    as_tibble() |>
    inner_join(kegg_categories, by = c("ID" = "brite")) |>
    rowwise() |>
    mutate(ratio = eval(parse(text = GeneRatio)) / eval(parse(text = BgRatio))) |>
    select(-Description) |>
    select(category, ratio, pvalue, p.adjust, qvalue, everything()) |>
    relocate(geneID, .after = last_col())
  
  write.csv(int_limma_res, paste("output/dataframes/df",filename,"genes.csv",sep = "_"))
  write.csv(int_enricher_res_tbl, paste("output/dataframes/df",filename,"pathways.csv",sep = "_"))
  write.csv(kegg_res_tbl, paste("output/dataframes/df",filename,"kegg.csv",sep = "_"))
  
}

medium_p27 <- do_medium_analysis(filename = "media_p27", reference_passage = "p27")
medium_p77 <- do_medium_analysis(filename = "media_p77", reference_passage = "p77")
