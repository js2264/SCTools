#' plotGeneExpression
#' 
#' Plot gene expression in embedded space
#'
#' @param sce 
#' @param dim 
#' @param genes 
#' @param q 
#' @param average_expr 
#' @param only_average 
#'
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotEmbedding <- function(sce, genes, dim = "UMAP", q = 0.95, average_expr = FALSE, only_average = FALSE) {
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
    if (all(length(genes) == 1 & genes %in% row.names(sce))) {
        df$gene <- bindByQuantiles(logcounts(sce)[genes, ], q_low = 1 - q, q_high = q)
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
    }
    
    # ---- Multiple genes
    if (all(length(genes) > 1 & genes %in% row.names(sce))) {
        for (gene in genes) {
            expr <- bindByQuantiles(logcounts(sce)[gene, ], q_low = 1 - q, q_high = q)
            df[, paste0(gene, "_expr")] <- expr
        }
        df_bound <- gather(df, "gene", "expr", -grep("_expr$", colnames(df), value = TRUE, invert = TRUE))
        plotFUN <- function(df, gene) {
            ggplot(df, aes_string(x = "Dim_1", y = "Dim_2", fill = "expr")) + 
                geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
                theme_bw() + 
                # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
                scale_fill_distiller(palette = 'YlOrBr', direction = 1) + 
                coord_fixed((max(df[, "Dim_1"])-min(df[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df[, "Dim_2"]))) + 
                labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = gene) + 
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.ticks = element_blank()
                )
        }
        pl <- df_bound %>% 
            group_by(gene) %>% 
            nest() %>% 
            mutate(plots = purrr::map2(data, gene, plotFUN)) %>% 
            '[['("plots")
        
        if (average_expr == TRUE) {
            df$gene <- bindByQuantiles(Matrix::colMeans(logcounts(sce)[genes, ]), q_low = 1 - q, q_high = q)
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
            pl[[length(pl) + 1]] <- p
        }
        
        if (!only_average) {
            p <- cowplot::plot_grid(plotlist = pl, align = 'hv')
        } else {
            p <- pl[[length(pl)]]
        }

    }
    return(p)
}

#' plotGE
#' 
#' Quickly plot gene expression of genes as boxplots
#'
#' @param sce 
#' @param genes 
#' @param by 
#'
#' @return plot
#' 
#' @import SingleCellExperiment
#' @import ggplot2
#'
#' @export

plotBoxplots <- function(sce, genes, by) {
    df <- data.frame(
        by = colData(sce)[[by]],
        t(as.matrix(logcounts(sce[genes,])))
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
#'
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotBoth <- function(sce, genes, by, dim = "UMAP", q = 0.95) {
    p1 <- plotEmbedding(sce, genes[[1]], dim, average_expr = FALSE, only_average = FALSE)
    p2 <- plotBoxplots(sce, genes[[1]], by)
    p <- cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
    return(p)
}
