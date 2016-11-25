#' Calculate distances of Tn5 insertions from motif.
#' 
#' @param df A data.frame returned from funciton calc_v.
#' @param ... Arguments for geom_point in ggplot.
#' @param size Numeric. Point size. 
#' @param alpha Numeric value between 0 and 1. Point transparency. 
#' @return An object of class 'ggplot'.
#' @import ggplot2
#' @export
plot_v <- function(df, size=1, alpha=0.01, ...){
        # plot the data
        gg <- ggplot2::ggplot(df, aes(x = df[ ,1], y = df[ ,2])) +        
                geom_point(alpha=alpha, size=size, ...) +
                scale_x_continuous(expand = c(0, 0)) + 
                scale_y_continuous(expand = c(0, 0)) +
                xlab("Distance from motif") +
                ylab("Fragment length") +
                theme_bw() +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank())
        
        return(gg)
}