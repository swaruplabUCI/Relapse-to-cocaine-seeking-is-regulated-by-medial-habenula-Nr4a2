
custom_vln <- function(
  seurat_obj,
  features,
  group.by,
  split.by=NULL,
  groups = NULL,
  selected_split = NULL,
  plot.margin = margin(0,0,0,0, "cm"),
  slot = 'data',
  assay = NULL,
  raster_dpi=200,
  add_boxplot = TRUE,
  pt.size = 0,
  add.noise = TRUE,
  line.size=NULL,
  adjust=1,
  quantiles=c(0.5),
  stat_method = 'wilcox.test',
  comparisons = NULL,
  pval_y_adjust=0.7,
  ref_group = '.all.',
  group_color_df = NULL,
  split_colors = NULL,
  add_colorbar = TRUE,
  plot_ymin = 0
){

  if(is.null(assay)){assay <- seurat_obj@active.assay}

  seurat_meta <- seurat_obj@meta.data

  if(class(seurat_meta[[group.by]]) != 'factor'){
    seurat_meta[[group.by]] <- as.factor(seurat_meta[[group.by]])
  }

  # get list of groups from seurat object
  if(is.null(groups)){
      groups <- as.character(unique(seurat_obj@meta.data[[group.by]]))
      groups <- groups[order(groups)]
  } else{
    if(!all(groups %in% seurat_obj@meta.data[[group.by]])){
      stop('Some groups not present in seurat_obj@meta.data[[group.by]]')
    }
  }

  if(!is.null(selected_split)){
    if(!all(selected_split %in% seurat_obj@meta.data[[split.by]])){
      stop('Some selected_split not present in seurat_obj@meta.data[[split.by]] ')
    }
  }

  # number of individual plots to make
  n_plots <- length(features) * length(groups)


  # colors for groups
  if(is.null(group_color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(seurat_meta[[group.by]])),
      group_name = levels(seurat_meta[[group.by]])
    )

    p <- seurat_meta %>%
      ggplot(aes_string(x=1, y=1, color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    group_color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    group_color_df$group <- factor_df$group_name

  } else if(is.character(group_color_df)){
    # check if it's a named character vec:
    if(is.null(names(group_color_df))){
      stop("color_df must be a dataframe with a group column and a colour column, or color_df can be a named character vector where the values are the colors and the names are the corresponding groups.")
    }
    group_color_df <- data.frame(
      colour = as.character(group_color_df),
      group = names(group_color_df)
    )
  }

  # only keep groups that we are using
  group_color_df <- subset(group_color_df, group %in% groups)
  group_color_df$var <- 1

    # color scheme
  group_colors <- group_color_df$colour
  names(group_colors) <- as.character(group_color_df$group)

  # makr color bar in ggplot
  if(add_colorbar){
    colorbar <- group_color_df %>%
      ggplot(aes(x=group, y=var, fill=group)) +
      geom_tile() +
      scale_fill_manual(values=group_color_df$colour) +
      NoLegend() +
      theme(
        plot.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.margin=margin(0,0,0,0)
      ) + RotatedAxis()
  }

  # initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_plots,style = 3, width = 50, char = "=")

  patch_list <- list()
  i <- 1
  for(feature in features){

    plot_df <- seurat_meta

    if(feature %in% colnames(plot_df)){
      plot_df$PlotFeature <- plot_df[[feature]]
    }
    else{
      plot_df$PlotFeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    }

    plot_range <- range(plot_df$PlotFeature)

    if(add.noise){
      noise <- rnorm(n = length(x = plot_df$PlotFeature)) / 100000
      plot_df$PlotFeature <- plot_df$PlotFeature + noise
    }

    # subset selected groups to compare only:
    if(!is.null(selected_split)){
      plot_df <- plot_df[plot_df[[split.by]] %in% selected_split,]
    }

    plot_list <- list()

    for(cur_group in groups){
      cur_df <- plot_df[as.character(plot_df[[group.by]]) == cur_group,]

      if(is.null(split.by)){
        p <- cur_df %>% ggplot(aes_string(x=group.by, y='PlotFeature', fill=group.by)) +
        scale_fill_manual(values=group_colors, na.value = 'grey90')
      } else{
        p <- cur_df %>% ggplot(aes_string(x=group.by, y='PlotFeature', fill=split.by))
        if(!is.null(split_colors)){
          p <- p + scale_fill_manual(values=split_colors)
        }
      }

      # add violin
      p <- p + geom_violin(
          trim=FALSE, adjust=adjust,
          scale='width', draw_quantiles = quantiles,
          color='black', lwd=0.5
        )

      if(add_boxplot){
        p <- p + geom_boxplot(width=0.1, fill='white', outlier.shape=NA)
      }

      if(pt.size > 0){
        p <- p +  ggrastr::rasterise(ggbeeswarm::geom_beeswarm(size=0.5), dpi=raster_dpi)
      }

      # code for stacked vln plot theme from CellChat:
      p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
        theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
      p <- p + theme(plot.title= element_blank(), # legend.position = "none",
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_text(size = rel(1), angle = 0),
                     axis.text.y = element_text(size = rel(1)),
                     plot.margin = plot.margin ) +
        theme(axis.text.y = element_text(size = 8))
      p <- p + theme(element_line(size=line.size))

      if(length(plot_list) > 0){
        p <- p + theme(
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
        )
      }


      # y limits:
      ymax <- ceiling(plot_range[2]); ymin <- floor(plot_range[1])
      p <- p + expand_limits(y=c(ymin, ymax)) +
        scale_y_continuous(expand=c(0,0), breaks = c(ymax), limits=c(plot_ymin, NA))

      # stats:
      if(!is.null(split.by)){
        if(is.null(comparisons)){
          p <- p + stat_compare_means(method=stat_method, aes(label = ..p.signif..), label.y=ymax*pval_y_adjust)
        } else{
          p <- p + stat_compare_means(
            method=stat_method,
            label='p.signif',
            comparisons=comparisons,
            label.y=ymax-1.5,
            ref.group=ref_group
          )
        }
      }

      p <- p + ylab(feature)

      # different settings for the last gene:
      if(!add_colorbar & feature == features[length(features)]){
        p <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
          RotatedAxis()
      }

      plot_list[[cur_group]] <- p

      # update progress bar
      setTxtProgressBar(pb, i)
      i <- i+1

    }
    patch_list[[feature]] <- wrap_plots(plot_list, ncol=length(groups))

  }

  # close progress bar
  close(pb)

  if(length(patch_list) == 1){
    return(patch_list[[1]])
  }

  plot_heights <- rep(50 / length(features), length(features))

  if(add_colorbar){
    plot_heights <- c(plot_heights, 1)
    patch_list <- c(patch_list, list(colorbar))
  }

  wrap_plots(patch_list, ncol=1) + plot_layout(heights=plot_heights, guides='collect')
}



FeatureEmbedding <- function(
  seurat_obj,
  features,
  raster = TRUE,
  slot = 'data',
  assay = NULL,
  reduction = 'umap',
  dpi = 200,
  dpi_scale=0.5,
  ncol = 4,
  combine=TRUE,
  order_points=TRUE, # 'shuffle'
  point_size = 0.5,
  plot_max = 'q100',
  plot_min = 'q0',
  same_range = FALSE,
  colfunc = viridis::inferno,
  rev_colors = TRUE
){

  input_plot_max <- plot_max
  input_plot_min <- plot_min

  plot_list <- list()

  if(same_range){
    tmp <- GetAssayData(seurat_obj, slot=slot, assay=assay)[features,]
    plot_range <- range(tmp)
    print(plot_range)

    if(is.null(plot_max)){plot_max <- plot_range[2]}
    if(is.null(plot_min)){plot_min <- plot_range[1]}
    print(plot_max)
    print(plot_min)
  }


  for(feature in features){

    plot_df <- seurat_obj@meta.data
    plot_df$plot_x_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$plot_y_coord <-  seurat_obj@reductions[[reduction]]@cell.embeddings[,2]


    # check if the feature is in the meta-data
    if(feature %in% rownames(seurat_obj)){
      if(is.null(assay)){assay <- seurat_obj@active.assay}
      plot_df$plotfeature <- GetAssayData(seurat_obj, slot=slot, assay=assay)[feature,]
    } else if(feature %in% colnames(plot_df)){
      if(!is.numeric(plot_df[[feature]])){
        stop("Specified feature is not numeric. Try plotting with VisDimPlot?")
      }
      plot_df$plotfeature <- plot_df[[feature]]
    } else{
      stop("feature not found in rownames(seurat_obj) or in colnames(seurat_obj@meta.data).")
    }

    if(!same_range){
      plot_max <- input_plot_max
      plot_min <- input_plot_min

      plot_range <- range(plot_df$plotfeature)
      if(!is.null(plot_max)){
        if(is.character(plot_max)){
          quant <- as.numeric(gsub('q', '', plot_max)) / 100
          plot_max <- as.numeric(quantile(plot_df$plotfeature, quant))
        }
        plot_range[2] <- plot_max
        print(plot_max)
        plot_df$plotfeature <- ifelse(
          plot_df$plotfeature > plot_max,
          plot_max,
          plot_df$plotfeature
        )
      }

      if(is.character(plot_min)){
        quant <- as.numeric(gsub('q', '', plot_min)) / 100
        plot_min <- as.numeric(quantile(plot_df$plotfeature, quant))
      }
      plot_range[1] <- plot_min
    }
    #
    # shuffle points:
    if(order_points == TRUE){
      plot_df <- plot_df %>% dplyr::arrange(plotfeature)
    } else if(order_points == "shuffle"){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }


    # initialize gpgplot
    p <- plot_df %>% subset(plotfeature > plot_min) %>%
      ggplot(aes_(x=~plot_x_coord, y=~plot_y_coord,color=~plotfeature))

    # add all the grey dots with low/zero expression
    if(raster){
      p <- p +
        ggrastr::rasterise(
          geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size),
          dpi=dpi, scale=dpi_scale
        )  +
        ggrastr::rasterise(
          geom_point(size=point_size),
          dpi=dpi, scale=dpi_scale
        )
    } else{
      p <- p +
        geom_point(inherit.aes=FALSE, data = subset(plot_df, plotfeature <= plot_min), aes(x=plot_x_coord, y=plot_y_coord),color='lightgrey', size=point_size) +
        geom_point(size=point_size)
    }

    # add extras to plot:
    colors <- colfunc(256)
    if(rev_colors){colors <- rev(colors)}
    p <- p +
      labs(color = feature) +
      scale_color_gradientn(colors=colors, limits = plot_range) +
      coord_equal() +
      theme(
        plot.title = element_text(hjust=0.5, face='plain'),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
      )

    plot_list[[feature]] <- p

  }

  if(length(plot_list) == 1){
    return(plot_list[[1]])
  }

  if(combine){
    patch <- wrap_plots(plot_list, ncol=ncol)
    return(patch)
  } else{
    return(plot_list)
  }


}




PlotEmbedding <- function(
  seurat_obj,
  group.by,
  reduction = 'umap',
  split.by = NULL,
  plot_under = FALSE,
  label = TRUE,
  point_size = 1,
  legend_point_size = 5,
  text_size = 3,
  raster = TRUE,
  order_points = "shuffle",
  x = NULL,
  y = NULL,
  raster_dpi = 400,
  raster_scale = 1,
  selected = NULL,
  plot_theme = NULL,
  color_df = NULL,
  plot_ratio=TRUE
){

  plot_df <- seurat_obj@meta.data
  if(!is.null(x) & !is.null(y)){
    plot_df$x <- plot_df[[x]]
    plot_df$y <- plot_df[[y]]
  } else if(reduction %in% names(seurat_obj@reductions)){
    plot_df$x <- seurat_obj@reductions[[reduction]]@cell.embeddings[,1]
    plot_df$y <- seurat_obj@reductions[[reduction]]@cell.embeddings[,2]
  }


  # convert to a factor:
  if(!is.factor(plot_df[[group.by]])){
    cur_groups <- unique(plot_df[[group.by]])
    cur_groups <- cur_groups[order(cur_groups)]
    plot_df[[group.by]] <- factor(
      as.character(plot_df[[group.by]]),
      levels = cur_groups
    )
  }

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  if(label){
    for(cur_cluster in unique(plot_df[[group.by]])){
      cur_meta <- plot_df[plot_df[[group.by]] == cur_cluster,]
      df <- data.frame(
        cluster = cur_cluster,
        x = mean(cur_meta$x),
        y = mean(cur_meta$y)
      )
      centroid_df <- rbind(centroid_df, df)
    }
  }

  # make a dummy ggplot to extract color scheme
  if(is.null(color_df)){

    factor_df <- data.frame(
      level = 1:length(levels(plot_df[[group.by]])),
      group_name = levels(plot_df[[group.by]])
    )

    p <- plot_df %>%
      ggplot(aes_string(x='x', y='y', color=group.by)) +
      geom_point()
    g <- ggplot_build(p)
    g_df <- g$data[[1]]
    color_df <- dplyr::select(g_df, c(colour, group)) %>% distinct() %>% arrange(group)
    color_df$group <- factor_df$group_name

  }

  # only show selected groups
  if(!is.null(selected)){
    plot_df[[group.by]][!(plot_df[[group.by]] %in% selected)] <- NA

    if(is.factor(plot_df[[group.by]])){
      plot_df[[group.by]] <- droplevels(plot_df[[group.by]])
    }

    if(label){
      centroid_df <- subset(centroid_df, cluster %in% selected)
    }
  }

  print(levels(plot_df[[group.by]]))

  # shuffle points:
  if(order_points == "shuffle"){
    plot_df <- plot_df[sample(nrow(plot_df)),]
  }

  # plot a single embedding
  if(is.null(split.by)){
    p <- .PlotSingleEmbedding(plot_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, plot_theme=plot_theme, plot_ratio=plot_ratio)
  } else{

    split_groups <- unique(plot_df[[split.by]])
    plot_list <- lapply(split_groups, function(cur_split){
      print(cur_split)
      cur_df <- plot_df[plot_df[[split.by]] == cur_split,]
      split_df <- plot_df[plot_df[[split.by]] != cur_split,]
      .PlotSingleEmbedding(cur_df, group.by, label, raster, raster_dpi, raster_scale, point_size, legend_point_size, text_size, color_df, centroid_df, split_df, cur_split, plot_theme, plot_under, plot_ratio=plot_ratio)
    })

    return(plot_list)

  }

  p

}

.PlotSingleEmbedding <- function(
  plot_df,
  group.by,
  label,
  raster,
  raster_dpi,
  raster_scale,
  point_size,
  legend_point_size,
  text_size,
  color_df,
  centroid_df,
  split_df = NULL,
  cur_split = NULL,
  plot_theme = NULL,
  plot_under=FALSE,
  plot_ratio = TRUE
){

  p <- plot_df %>%inc
    ggplot(aes_string(x='x', y='y', color=group.by))

  # add the under layer for the split plots
  if(!is.null(split_df) & plot_under){
    # add points
    if(raster){
      p <-  p + ggrastr::rasterise(geom_point(data = split_df, size=point_size/2, color='lightgrey'), dpi=raster_dpi, scale=raster_scale)
    } else{
      p <- p + geom_point(data = split_df, size=point_size/2, color='lightgrey')
    }
  }

  # add points
  if(raster){
    p <-  p + ggrastr::rasterise(geom_point(size=point_size), dpi=raster_dpi,  scale=raster_scale)
  } else{
    p <- p + geom_point(size=point_size)
  }

  # add labels
  if(label){
    p <- p + ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=text_size)
  }

  # color scheme
  plot_colors <- color_df$colour
  names(plot_colors) <- as.character(color_df$group)
  p <- p + scale_color_manual(values=plot_colors, na.value = 'grey90')

  # add title:
  if(is.null(cur_split)){
    p <- p + ggtitle(group.by) + theme(plot.title=element_text(hjust=0.5))
  } else{
    p <- p + ggtitle(cur_split) + theme(plot.title=element_text(hjust=0.5))
  }

  # add theme:
  if(!is.null(plot_theme)){
    p <- p + plot_theme
  }

  # fixed coords
  if(plot_ratio){
    p <- p + coord_equal()
  }

  # adjust legend point size
  p <- p + guides( color = guide_legend(override.aes = list(size=legend_point_size)))

  p

}
