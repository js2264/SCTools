#' plotEmbedding
#' 
#' Plot gene expression in embedded space
#'
#' @param sce 
#' @param dim 
#' @param by 
#' @param q 
#' @param average_expr 
#' @param only_average 
#' @param assay.type 
#' @param theme.args 
#' @param return_plotlist 
#'
#' @import SingleCellExperiment
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotEmbedding <- function(sce, by, dim = "UMAP", q = 0.95, average_expr = FALSE, only_average = FALSE, assay.type = 'logcounts', theme.args = NULL, return_plotlist = FALSE) {
    
    if (any(by %in% rownames(sce))) {
        bygene <- TRUE 
        .checkGenes(sce, by)
    }
    else if (any(by %in% colnames(colData(sce)))) {
        bygene <- FALSE 
        .checkColData(sce, by)
    }
    else {
        stop("'by' argument not found in genes nor colData. Aborting. ")
    }

    .checkEmbedding(sce, dim)

    df <- data.frame(colData(sce))
    df[, "Dim_1"] <- reducedDim(sce, dim)[, 1]
    df[, "Dim_2"] <- reducedDim(sce, dim)[, 2]
    
    # ---- Only 1 gene
    if (length(by) == 1) {
        if (bygene) {
            df$by <- bindByQuantiles(assay(sce, assay.type)[by, ], q_low = 1 - q, q_high = q)
            scale <- scale_fill_distiller(palette = 'YlOrBr', direction = 1)
        } 
        else {
            if (length(unique(colData(sce)[[by]])) >= 12) {
                df$by <- colData(sce)[[by]]
                scale <- scale_fill_distiller(palette = 'YlOrBr', direction = 1)
            }
            else {
                df$by <- factor(colData(sce)[[by]])
                scale <- scale_fill_brewer(palette = 'Spectral', direction = 1)
            }
       }
        p <- ggplot(df, aes(x = Dim_1, y = Dim_2, fill = by)) + 
            geom_point(pch = 21, alpha = 0.5, col = '#bcbcbc', stroke = 0.2) + 
            theme_bw() + 
            # scale_fill_gradient(low = 'white', high = '#8b4e36') + 
            scale + 
            coord_fixed((max(df[, "Dim_1"])-min(df[, "Dim_1"]))/(max(df[, "Dim_2"])-min(df[, "Dim_2"]))) + 
            labs(y = paste(dim, " 2"), x = paste(dim, " 1"), fill = by) + 
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.ticks = element_blank()
            )
        if (!is.null(theme.args)) p <- p + theme.args
    }
    
    # ---- Multiple by
    if (length(by) > 1) {
        for (gene in by) {
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
            df$gene <- bindByQuantiles(Matrix::colMeans(assay(sce, assay.type)[by, ]), q_low = 1 - q, q_high = q)
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
#' @import SingleCellExperiment
#' @import dplyr
#' @import tidyr
#' 
#' @return plot
#'
#' @export

plotAnimatedEmbedding <- function(sce, genes, dim = "UMAP", q = 0.95, assay.type = 'logcounts', theme.args = NULL) {
    
    .checkGenes(sce, genes)

    df <- data.frame(colData(sce))
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
    p <- p + 
    gganimate::transition_states(
        gene,
        transition_length = 1,
        state_length = 1
    ) + 
    ggtitle('Now showing {closest_state}') + 
    gganimate::enter_fade() + 
    gganimate::exit_fade() 
    return(p)
}

#' plotBoxplots
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
    .checkColData(sce, by)

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

#' plotBarplots
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

plotBarplots <- function(sce, genes, color_by = 'annotation', order_by = NULL, bins = NULL, assay.type = 'logcounts') {
    
    .checkGenes(sce, genes)
    .checkColData(sce, color_by)
    .checkColData(sce, order_by)
    .checkColData(sce, 'Barcode')

    `%<>%` <- magrittr::`%<>%`

    df <- data.frame(
        color_by = colData(sce)[[color_by]],
        order_by = colData(sce)[[order_by]],
        cell = colData(sce)[['Barcode']], 
        t(as.matrix(assay(sce, assay.type)[genes,]))
    ) %>% 
        arrange(color_by) %>% 
        gather('gene', 'expr', -color_by, -order_by, -cell) %>% 
        filter(!is.na(color_by)) %>% 
        mutate(gene = factor(gene, genes))

    if (!is.null(order_by)) {
        df %<>% arrange(order_by)
    }
    df %<>% mutate(cell = factor(cell, unique(cell)))

    if (!is.null(bins)) {
        df <- df %>% 
            mutate(bins = cut(1:nrow(.), breaks = seq(1, nrow(.), length.out = bins), include.lowest = TRUE)) %>% 
            group_by(gene, bins) %>% 
            summarize(expr = mean(expr)) %>% 
            mutate(cell = bins) %>% 
            mutate(color_by = 1)
    }

    p <- ggplot(df, aes(x = cell, y = expr, fill = color_by, alpha = expr)) +  
        geom_col(width=1) + 
        theme_bw() + 
        labs(y = "log counts", x = "", fill = "Sets") + 
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()
        ) + 
        facet_grid(rows = "gene", scales = 'free')
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
