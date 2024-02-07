# MotifScan <- function(
#   seurat_obj,
#   species_genome, # hg38, mm10, etc...
#   pfm, # matrix set from JASPAR2020 for example
#   EnsDb, # Ensembl database such as EnsDb.Mmusculus.v79
#   wgcna_name=NULL
# ){

#   if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}

#   # get a dataframe of just the motif name and the motif ID:
#   motif_df <- data.frame(
#     motif_name = purrr::map(1:length(pfm), function(i){pfm[[i]]@name}) %>% unlist,
#     motif_ID = purrr::map(1:length(pfm), function(i){pfm[[i]]@ID}) %>% unlist
#   )

#   # get promoters of protein coding genes from the given Ensdb
#   # note: everything breaks if I try to use X & Y chromosomes.
#   # gene.promoters <- ensembldb::promoters(EnsDb, filter = ~ gene_biotype == "protein_coding") %>%
#   #   subset(seqnames %in% c(1:100))
#   # gene.coords <- ensembldb::genes(EnsDb, filter = ~ gene_biotype == "protein_coding") %>%
#   #   subset(seqnames %in% c(1:100))

#   gene.promoters <- ensembldb::promoters(EnsDb) %>%
#     subset(seqnames %in% c(1:100))
#   gene.coords <- ensembldb::genes(EnsDb) %>%
#     subset(seqnames %in% c(1:100))


#   # add the gene name to the promoter object
#   gene.promoters$symbol <- gene.coords$symbol[match(gene.promoters$gene_id, names(gene.coords))]

#   # drop unnecessary chromosomes
#   gene.promoters <- GenomeInfoDb::keepSeqlevels(gene.promoters, value= levels(droplevels(seqnames(gene.promoters))))

#   # rename seqlevels to add 'chr', remove X&Y chromosomes because they break everything
#   old_levels <- levels(GenomeInfoDb::seqnames(gene.promoters))
#   new_levels <- ifelse(old_levels %in% c('X', 'Y'), old_levels, paste0('chr', old_levels))
#   gene.promoters <- GenomeInfoDb::renameSeqlevels(gene.promoters, new_levels)

#   # set the genome (not sure if we NEED to do this...)
#   genome(seqinfo(gene.promoters)) <- species_genome

#   # set up promoters object that only has the necessary info for motifmatchr
#   my_promoters <- GenomicRanges::GRanges(
#     seqnames =  droplevels(seqnames(gene.promoters)),
#     IRanges(
#       start = start(gene.promoters),
#       end = end(gene.promoters)
#     ),
#     symbol = gene.promoters$symbol,
#     genome=species_genome
#   )

#   # scan these promoters for motifs:
#   print('Matching motifs...')
#   motif_ix <- motifmatchr::matchMotifs(pfm, my_promoters, genome=species_genome)

#   # get the matches
#   tf_match <- motifmatchr::motifMatches(motif_ix)
#   rownames(tf_match) <- my_promoters$symbol

#   # use motif names as the column names:
#   colnames(tf_match) <- motif_df$motif_name

#   # only keep genes that are in the Seurat object and in the given EnsDb:
#   gene_list <- rownames(seurat_obj)
#   gene_list <- gene_list[gene_list %in% rownames(tf_match)]
#   tf_match <- tf_match[gene_list,]

#   # get list of target genes for each TF:
#   print('Getting TF target genes...')
#   tfs <- motif_df$motif_name
#   tf_targets <- list()
#   n_targets <- list()
#   for(cur_tf in tfs){
#     tf_targets[[cur_tf]] <- names(tf_match[,cur_tf][tf_match[,cur_tf]])
#     n_targets[[cur_tf]] <- length(tf_targets[[cur_tf]] )
#   }
#   n_targets <- unlist(n_targets)

#   # add number of target genes to motif df
#   motif_df$n_targets <- n_targets

#   # add info to seurat object
#   seurat_obj <- SetMotifMatrix(seurat_obj, tf_match)
#   seurat_obj <- SetMotifs(seurat_obj, motif_df)
#   seurat_obj <- SetMotifTargets(seurat_obj, tf_targets)
#   seurat_obj <- SetPFMList(seurat_obj, pfm)

#   seurat_obj
# }

ConstructTFNetwork <- function(
    seurat_obj,
    model_params,
    nfold=5,
    wgcna_name=NULL
){

    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # TODO: add checks? 
   
    # get the motif information from the Seurat object:
    motif_matrix <- GetMotifMatrix(seurat_obj)
    motif_df <- GetMotifs(seurat_obj)

    # check that gene names column is in the motif names
    if(! 'gene_name' %in% colnames(motif_df)){
        stop('gene_name column missing in motif table (GetMotifs(seurat_obj)). Please add a column indicating the gene_name in the seurat_obj for each motif.' )
    }

    # subset the motif_df by genes that are in the seurat obj:
    motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

    # get the expression matrix:
    datExpr <- as.matrix(GetDatExpr(seurat_obj, wgcna_name=wgcna_name))
    genes_use <- colnames(datExpr)

    # set up output dataframes
    importance_df <- data.frame()
    eval_df <- data.frame()

    # set up the progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(genes_use), style = 3, width = 50, char = "=")
    counter <- 1

    ts <- tic()
    for(cur_gene in genes_use){

      print(cur_gene)

        setTxtProgressBar(pb, counter)

        # check if this gene is in the motif matrix:
        if(! cur_gene %in% rownames(motif_matrix)){
            print(paste0('Gene not found in the motif_matrix, skipping ', cur_gene))
            next
        }

        # get the list of TFs that regulate this gene:
        cur_tfs <- names(which(motif_matrix[cur_gene,]))
        cur_tfs <- subset(motif_df, motif_name %in% cur_tfs) %>% .$gene_name %>% unique
        cur_tfs <- cur_tfs[cur_tfs %in% genes_use]

        # set up the expression matrices
        if(cur_gene %in% cur_tfs){
            cur_tfs <- cur_tfs[cur_tfs != cur_gene]
            x_vars <- datExpr[,cur_tfs]
        } 
        x_vars <- datExpr[,cur_tfs]
        y_var <- as.numeric(datExpr[,cur_gene])

        if(length(cur_tfs) < 2){
            print(paste0('Not enough putative TFs, skipping ', cur_gene))
            next
        }

        # correlation:
        tf_cor <- as.numeric(cor(
            x=as.matrix(x_vars),
            y=y_var
        ))
        names(tf_cor) <- cur_tfs

        # run xgboost model
        if(all(y_var == 0)){
            print(paste0('skipping ', cur_gene))
            next
        }
        xgb <- xgboost::xgb.cv(
            params = model_params,
            data = x_vars,
            label = y_var,
            nrounds = 100,
            showsd = FALSE,
            nfold = nfold,
            callbacks = list(cb.cv.predict(save_models=TRUE)),
            verbose=FALSE
        )

        # get the CV evaluation info
        xgb_eval <- as.data.frame(xgb$evaluation_log)
        xgb_eval$variable <- cur_gene

        # average the importance score from each fold
        importance <- Reduce('+', lapply(1:nfold, function(i){
            cur_imp <- xgb.importance(feature_names = colnames(x_vars), model = xgb$models[[i]])
            ix <- match(colnames(x_vars),  as.character(cur_imp$Feature))
            cur_imp <- as.matrix(cur_imp[ix,-1])
            cur_imp[is.na(cur_imp)] <- 0
            cur_imp
        })) / nfold
        importance <- as.data.frame(importance)

        # add tf and source info
        importance$tf <- colnames(x_vars)
        importance$gene<- cur_gene

        # add the tf correlation information
        importance$Cor <- as.numeric(tf_cor)

        # re-order columns, and re-order rows by gain:
        importance <- importance %>% dplyr::select(c(tf, gene, Gain, Cover, Frequency, Cor))
        importance <- arrange(importance, -Gain)

        # append
        importance_df <- rbind(importance_df, importance)
        eval_df <- rbind(eval_df, xgb_eval)

        # update progress bar
        counter <- counter+1

    } 

    # close the progress bar
    te <- toc()
    (te$toc - te$tic) / 60 
    close(pb)

    # add the results to the Seurat object
    seurat_obj <- SetTFNetwork(seurat_obj, importance_df, wgcna_name=wgcna_name)
    seurat_obj <- SetTFEval(seurat_obj, eval_df, wgcna_name=wgcna_name)

    #return(list(importance=importance_df, eval=eval_df))
    seurat_obj

}


AssignTFRegulons <- function(
    seurat_obj,
    strategy = "A", # A, B, or C
    reg_thresh = 0.01,
    n_genes = 50,
    n_tfs = 10,
    wgcna_name=NULL
){

    # get data from active assay if wgcna_name is not given
    if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
    CheckWGCNAName(seurat_obj, wgcna_name)

    # TODO: Checks 

    # get the tf_net from the seurat_obj 
    tf_net <- GetTFNetwork(seurat_obj, wgcna_name)

    if(strategy == 'A'){

        # Take the top TFs for each gene
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            group_by(gene) %>%
            slice_max(order_by=Gain, n=n_tfs) %>% 
            ungroup()

    } else if(strategy == 'B'){

        # Take the top target genes for each TF
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) %>% 
            group_by(tf) %>%
            slice_max(order_by=Gain, n=n_genes) %>% 
            ungroup()

    } else if(strategy == 'C'){

        # Take all interactions above a certain score
        tf_regulons <- tf_net %>% 
            subset(Gain >= reg_thresh) 

    } else {
        stop('Invalid choice for strategy. Valid choices are A, B, or C.')
    }

    # add the regulons to the seurat object
    seurat_obj <- SetTFRegulons(seurat_obj, tf_regulons, wgcna_name)

}


############################
# getters and setters
###########################

#' SetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param tf_net dataframe storing the TF network info in ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFNetwork <- function(seurat_obj, tf_net, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_net <- tf_net
  seurat_obj
}

#' GetTFNetwork
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFNetwork <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_net 
}


#' SetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param tf_eval dataframe storing the TF network evaluation info from ConstructTFNetwork
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFEval <- function(seurat_obj, tf_eval, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_eval <- tf_eval
  seurat_obj
}

#' GetTFEval
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFEval <- function(seurat_obj, wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_eval
}

#' SetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param tf_regulons dataframe storing the TF regulon info from AssignTFRegulons
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
SetTFRegulons <- function(seurat_obj, tf_regulons, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)

  seurat_obj@misc[[wgcna_name]]$tf_regulons <- tf_regulons
  seurat_obj
}

#' GetTFRegulons
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetTFRegulons <- function(seurat_obj, wgcna_name=NULL){

  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  CheckWGCNAName(seurat_obj, wgcna_name)
  seurat_obj@misc[[wgcna_name]]$tf_regulons

}




#' GetAvgModuleExpr
#'
#' @param seurat_obj A Seurat object
#' @param wgcna_name The name of the hdWGCNA experiment in the seurat_obj@misc slot
#' @keywords scRNA-seq
#' @export
GetAvgModuleExpr <- function(seurat_obj,  wgcna_name=NULL){
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  seurat_obj@misc[[wgcna_name]]$avg_modules
}
