#' plotEmbedding
#' 
#' Plot gene expression in embedded space
#'
#' @param sce 
#' @param dim 
#' @param genes 
#' @param q 
#' @param average_expr 
#' @param only_average 
#' @param assay.type
#' @param theme.args
#' @param return_plotlist
#'
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotEmbedding <- function(sce, genes, dim = "UMAP", q = 0.95, average_expr = FALSE, only_average = FALSE, assay.type = 'logcounts', theme.args = NULL, return_plotlist = FALSE) {
    
    .checkGenes(sce, genes)

    df <- data.frame(
        colData(sce), 
        PCA_X = reducedDim(sce, 'PCA')[, 1], 
        PCA_Y = reducedDim(sce, 'PCA')[, 2], 
        TSNE_X = reducedDim(sce, 'TSNE')[, 1], 
        TSNE_Y = reducedDim(sce, 'TSNE')[, 2], 
        UMAP_X = reducedDim(sce, 'UMAP')[, 1], 
        UMAP_Y = reducedDim(sce, 'UMAP')[, 2], 
        force_X = reducedDim(sce, 'force')[, 1], 
        force_Y = reducedDim(sce, 'force')[, 2]
    )
    df[, "Dim_1"] <- df[, paste0(dim, "_X")]
    df[, "Dim_2"] <- df[, paste0(dim, "_Y")]
    
    # ---- Only 1 gene
    if (all(length(genes) == 1)) {
        df$gene <- bindByQuantiles(assay(sce, assay.type)[genes, ], q_low = 1 - q, q_high = q)
        p <- ggplot(df, aes_string(x = "Dim_1", y = "Dim_2", fill = "gene")) + 
            geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
            theme_bw() + 
            # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
            scale_fill_distiller(palette = 'YlOrBr', direction = 1) + 
            coord_fixed((max(df[, "Dim_1"])-min(df[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df[, "Dim_2"]))) + 
            labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = genes) + 
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.ticks = element_blank()
            )
        if (!is.null(theme.args)) p <- p + theme.args
    }
    
    # ---- Multiple genes
    if (all(length(genes) > 1)) {
        for (gene in genes) {
            expr <- bindByQuantiles(assay(sce, assay.type)[gene, ], q_low = 1 - q, q_high = q)
            df[, paste0(gene, "_expr")] <- expr
        }
        df_bound <- gather(df, "gene", "expr", -grep("_expr$", colnames(df), value = TRUE, invert = TRUE))
        plotFUN <- function(df, gene, theme.args) {
            p <- ggplot(df, aes_string(x = "Dim_1", y = "Dim_2", fill = "expr")) + 
                geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
                theme_bw() + 
                # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
                scale_fill_distiller(palette = 'YlOrBr', direction = 1) + 
                coord_fixed((max(df[, "Dim_1"])-min(df[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df[, "Dim_2"]))) + 
                labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = '', title = gsub('_expr', '', gene)) + 
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.ticks = element_blank()
                )
            if (!is.null(theme.args)) p <- p + theme.args
            return(p)
        }
        pl <- df_bound %>% 
            group_by(gene) %>% 
            nest() %>% 
            mutate(plots = purrr::map2(data, gene, plotFUN, theme.args)) %>% 
            # mutate(plots = purrr::map2(data, gene, plotFUN)) %>% 
            pull(plots)
        
        if (average_expr == TRUE) {
            df$gene <- bindByQuantiles(Matrix::colMeans(assay(sce, assay.type)[genes, ]), q_low = 1 - q, q_high = q)
            p <- ggplot(df, aes_string(x = "Dim_1", y = "Dim_2", fill = "expr")) + 
                    geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
                    theme_bw() + 
                    # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
                    scale_fill_distiller(palette = 'YlOrBr', direction = 1) + 
                    coord_fixed((max(df[, "Dim_1"])-min(df[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df[, "Dim_2"]))) + 
                    labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = "Average expr.") + 
                    theme(
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), 
                        axis.ticks = element_blank()
                    )
            if (!is.null(theme.args)) p <- p + theme.args
            pl[[length(pl) + 1]] <- p
        }
        
        if (!only_average) {
            if (return_plotlist) {
                p <- pl
            } else {
                p <- cowplot::plot_grid(plotlist = pl, align = 'hv')
            }
        } else {
            p <- pl[[length(pl)]]
        }

    }
    return(p)
}

#' plotAnimatedEmbedding
#' 
#' Plot gene expression in embedded space
#'
#' @param sce 
#' @param dim 
#' @param genes 
#' @param q 
#' @param assay.type
#' @param theme.args
#'
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotAnimatedEmbedding <- function(sce, genes, dim = "UMAP", q = 0.95, assay.type = 'logcounts', theme.args = NULL) {
    
    .checkGenes(sce, genes)

    df <- data.frame(
        colData(sce), 
        PCA_X = reducedDim(sce, 'PCA')[, 1], 
        PCA_Y = reducedDim(sce, 'PCA')[, 2], 
        TSNE_X = reducedDim(sce, 'TSNE')[, 1], 
        TSNE_Y = reducedDim(sce, 'TSNE')[, 2], 
        UMAP_X = reducedDim(sce, 'UMAP')[, 1], 
        UMAP_Y = reducedDim(sce, 'UMAP')[, 2], 
        force_X = reducedDim(sce, 'force')[, 1], 
        force_Y = reducedDim(sce, 'force')[, 2]
    )
    df[, "Dim_1"] <- df[, paste0(dim, "_X")]
    df[, "Dim_2"] <- df[, paste0(dim, "_Y")]

    for (gene in genes) {
        expr <- bindByQuantiles(assay(sce, assay.type)[gene, ], q_low = 1 - q, q_high = q)
        df[, paste0(gene, "_expr")] <- expr
    }
    df_bound <- gather(df, "gene", "expr", -grep("_expr$", colnames(df), value = TRUE, invert = TRUE))
    df_bound$gene <- factor(df_bound$gene, levels = paste0(genes, "_expr"))

    p <- ggplot(df_bound, aes_string(x = "Dim_1", y = "Dim_2", fill = "expr")) + 
        geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
        theme_bw() + 
        # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
        scale_fill_distiller(palette = 'YlOrBr', direction = 1) + 
        coord_fixed((max(df_bound[, "Dim_1"])-min(df_bound[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df_bound[, "Dim_2"]))) + 
        labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = '') + 
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.ticks = element_blank()
        )
    if (!is.null(theme.args)) p <- p + theme.args
    p <- p + gganimate::transition_states(
        gene,
        transition_length = 1,
        state_length = 1
    )
    return(p)
}

#' plotGE
#' 
#' Quickly plot gene expression of genes as boxplots
#'
#' @param sce 
#' @param genes 
#' @param by 
#' @param assay.type
#'
#' @return plot
#' 
#' @import SingleCellExperiment
#' @import ggplot2
#'
#' @export

plotBoxplots <- function(sce, genes, by = 'annotation', assay.type = 'logcounts') {
    
    .checkGenes(sce, genes)

    df <- data.frame(
        by = colData(sce)[[by]],
        t(as.matrix(assay(sce, assay.type)[genes,]))
    ) %>% 
        gather('gene', 'expr', -by) %>% 
        filter(!is.na(by)) %>% 
        mutate(gene = factor(gene, genes))
    p <- ggplot(df, aes(x = by, y = expr, fill = by)) +  
        geom_boxplot(outlier.shape = 19, outlier.fill = NA) + 
        theme_bw() + 
        labs(y = "log counts", x = "", fill = "Sets") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    if (length(genes) > 1) p <- p + facet_wrap("gene") 
    return(p)
}

#' plotBoth
#' 
#' Plot gene expression in embedded space and a boxplot below
#'
#' @param sce 
#' @param genes 
#' @param by 
#' @param dim 
#' @param q 
#' @param assay.type
#'
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotBoth <- function(sce, genes, by, dim = "UMAP", q = 0.95, assay.type = 'logcounts') {
    p1 <- plotEmbedding(sce, genes[[1]], dim, average_expr = FALSE, only_average = FALSE)
    p2 <- plotBoxplots(sce, genes[[1]], by)
    p <- cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
    return(p)
}
