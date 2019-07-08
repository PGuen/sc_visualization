# Visualization of signature enrichment scores  
sc_signature_enrichment <- 
  function(seurat_object,
           signature_upregulated,signature_downregulated,
           plot_type = "dim.reduction_plot", # {'dim.reduction_plot','box_plot','ridges_plot'}
           reduction_use="tsne", # must be slot of seurat_object@dr
           cell_ids = NULL, # must be colnames of seurat_object@data
           anno_name = "ident", # either ident or name of a character column from seurat_object@meta.data
           do.scale=T,pt.size=1, group.cols=NULL,plot.legends=F,
           scales_x_axis=NULL, # numeric vector of lower limit and upper limit for x axis
           scales_gradient=NULL, # numeric vector of lower limit and upper limit for color
           do.raster=F,
           seurat.assay=NULL
  ){
    # check some requirements
    packages_to_load=c("RColorBrewer","ggplot2","scales","ggridges","dplyr","ggrastr")
    if(any(!packages_to_load %in% installed.packages())) stop(paste("Please install the",packages_to_load[!packages_to_load %in% installed.packages()]))
    if(!as.character(class(seurat_object))=="Seurat") stop("Please provide an object of class seurat")
    if(!plot_type %in% c('dim.reduction_plot','box_plot','violin_plot','ridges_plot')) stop("Please use one of the following plot_types: dim.reduction_plot, box_plot, violin_plot ,ridges_plot")
    if(plot_type == 'dim.reduction_plot' & !reduction_use %in% names(seurat_object@reductions)) stop(paste("If dim.reduction_plot is specified,",reduction_use,"has to be computed"))
    if(!anno_name %in% c("ident",colnames(seurat_object@meta.data))) stop("anno_name has to be ident or column name of seurat_object@meta.data")
    if(is.null(seurat.assay)) seurat.assay <- DefaultAssay(seurat_object)
    if(is.null(cell_ids)) cell_ids <- colnames(seurat_object@assays[[seurat.assay]]@data)
    
    # load required packages
    invisible(lapply(packages_to_load, require, character.only = TRUE))  
    
    # fetch raw data
    count_data <- Matrix::as.matrix(seurat_object@assays[[seurat.assay]]@counts)[,cell_ids] 
    norm_data <- Matrix::as.matrix(seurat_object@assays[[seurat.assay]]@data)[,cell_ids] 
    
    # remove genes from signatures that are not present in the filtered dataset and report
    size_up = length(signature_upregulated);size_down= length(signature_downregulated)
    signature_upregulated <- signature_upregulated[signature_upregulated %in%
                                                     rownames(norm_data)]
    print(paste("Signature.upregulated:",length(signature_upregulated %in% 
                                                  rownames(norm_data)),"out of", size_up,
                " genes are contained in the dataset"))
    signature_downregulated <- signature_downregulated[signature_downregulated %in%
                                                         rownames(norm_data)]
    print(paste("Signature.downregulated:",length(signature_downregulated %in%
                                                    rownames(norm_data)),"out of", size_down,
                " genes are contained in the dataset"))
    # calculate mean expression for the signature_upregulated
    expr_df = data.frame(
      mean_upregulated=colMeans(count_data[as.character(signature_upregulated),,drop=F],na.rm = T),
      mean_downregulated=colMeans(count_data[as.character(signature_downregulated),,drop=F],na.rm = T),
      mean_total=colMeans(count_data,na.rm = T)
    ) 
    
    
    if (anno_name == "ident") {
      expr_df$group=as.character(seurat_object@active.ident[rownames(expr_df)])
    } else if (is.element(anno_nam,colnames(seurat_object@meta.data))) {
      expr_df$group=as.character(seurat_object$orig.ident)
    }

        # calculate fractions
    # since we use the raw data the counts are scaled by total counts
    expr_df$fraction_upregulated <- expr_df$mean_upregulated/expr_df$mean_total*100
    expr_df$fraction_downregulated <- expr_df$mean_downregulated/expr_df$mean_total*100
    expr_df$fraction_diff <- expr_df$fraction_upregulated-expr_df$fraction_downregulated
    
    if (is.null(signature_downregulated)) expr_df$fraction_diff <- expr_df$fraction_upregulated
    
    if (do.scale==T){
      expr_df$fraction_diff <- as.vector(scale(expr_df$fraction_diff))
    }
    
    if (is.null(group.cols)){
      group.cols <- RColorBrewer::brewer.pal(n = length(unique(expr_df$group)),name = "Paired")
    }
    
    
    ########### plotting 
    plot= ggplot(data = expr_df)
    
    # boxplot
    if (plot_type== "box_plot") {
      plot= plot +
        geom_boxplot(aes(x=group,y=fraction_diff,fill=group))+
        scale_fill_manual(values = group.cols)+
        xlab("Cell identity") + ylab("Scaled signature enrichment")
    } 
    
    # violin plot
    if (plot_type== "violin_plot") {
      plot= plot +
        geom_violin(aes(x=group,y=fraction_diff,fill=group))+
        geom_boxplot(aes(x=group,y=fraction_diff),width=.1) +
        scale_fill_manual(values = group.cols) +
        xlab("Cell identity") + ylab("Scaled signature enrichment")
    }
    
    # ridges plot
    if (plot_type== "ridges_plot") {
      plot= plot +
        geom_density_ridges_gradient(aes(x = fraction_diff, y= group ,
                                         fill=fraction_diff,fill = ..x..),
                                     scale = 2,size=.1) +  
        scale_y_discrete(expand = c(.01,0)) +
        scale_fill_distiller(name = paste("Scaled\nsignature\nenrichment"),
                             palette = "RdBu",limits=c(-2,2),oob=squish) +
        xlab("Scaled signature enrichment") +ylab("Cell identity")
      if(!is.null(scales_x_axis)){
        plot=plot+ coord_cartesian(xlim = c(scales_x_axis[1],scales_x_axis[2]))
      }
    }
    
    # dim.reduxtion plot 
    if (plot_type == "dim.reduction_plot"){
      dim.reduction_df = as.data.frame(seurat_object@reductions[[reduction_use]]@cell.embeddings[rownames(expr_df),])
      axes= colnames(dim.reduction_df)[1:2]
      expr_df= merge(expr_df, dim.reduction_df, by='row.names', all=TRUE) 
      expr_df <- expr_df[order(expr_df$fraction_diff,decreasing = F),]
      plot= ggplot(data = expr_df) 
      if(do.raster==T) {
        plot= plot + geom_point_rast(aes(x=expr_df[[axes[1]]],y=expr_df[[axes[2]]],
                                         color=fraction_diff),size = pt.size) 
      }
      if(do.raster==F) {
        plot= plot + geom_point(aes(x=expr_df[[axes[1]]],y=expr_df[[axes[2]]],
                                    color=fraction_diff),size = pt.size) 
      }
      plot= plot + 
        scale_colour_distiller(name = paste("Scaled\nsignature\nenrichment"),
                               palette = "RdBu",limits=c(-2,2),oob=squish) +
        xlab(axes[1])+ ylab(axes[2])
      
      if (!is.null(scales_gradient)) {
        plot=plot+scale_colour_distiller(name = paste("Scaled\nsignature\nenrichment"),
                                         palette = "RdBu",
                                         limits=c(scales_gradient[1],scales_gradient[2]),oob=squish)
      }
    }
    
    plot = plot +
      theme_bw() + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         text = element_text(size = 10)) 
    if (!plot.legends) plot = plot + theme(legend.position = "none")
    return(plot)
  }

# not run
# sc_signature_enrichment(seurat_object = combined_bal_minf50,
#                         signature_upregulated = c("MARCO","LYZ") ,
#                         signature_downregulated = NULL,
#                         plot_type = "dim.reduction_plot",
#                         reduction_use = "umap")
# 
# sc_signature_enrichment(seurat_object = combined_bal_minf50,
#                         signature_upregulated = c("MARCO","LYZ") ,
#                         signature_downregulated = NULL,
#                         plot_type = "box_plot")
# 
# sc_signature_enrichment(seurat_object = combined_bal_minf50,
#                         signature_upregulated = c("MARCO","LYZ") ,
#                         signature_downregulated = NULL,
#                         plot_type = "violin_plot")
# 
# sc_signature_enrichment(seurat_object = combined_bal_minf50,
#                         signature_upregulated = c("MARCO","LYZ") ,
#                         signature_downregulated = NULL,
#                         plot_type = "ridges_plot")



