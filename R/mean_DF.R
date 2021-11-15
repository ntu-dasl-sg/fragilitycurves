#' Compute the mean damage factors and plot mean damage factor against peak ground acceleration (PGA).
#'
#' @param model An ordinal model fitted for \code{data} using \code{polr} from the R package \code{MASS}.
#' @param data A dataframe with the columns \code{CDF} (ordered factor), \code{logPGA} (numeric), \code{PGA} (numeric), \code{Easting} (numeric) and \code{Northing} (numeric).
#' @param upper.bin A vector containing the unique upper bounds of the damage bins. This should have length one less than the number of damage states.
#' @param bin.length A vector containing the length of each damage bin. This should have length equal to the number of damage states, and start and end with 1 since we assume point masses at the damage bin according to zero and complete damage.
#' @param ex.prob A dataframe with the columns \code{PGA}, \code{CDF} and \code{ex.prob} (the exceedance probability of the \code{CDF} level at the \code{PGA} value) such as the output of \code{frag_curve}. This is only required if we wish to plot the mean damage factor curve, i.e. set \code{plot = TRUE}.
#' @param plot A logical value indicating if the example probability density and mean damage factor curve should be plotted.
#' @return A vector of estimated mean damage factors corresponding to the rows in \code{data}.
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' library(MASS)
#'
#' data(damage_simulation)
#'
#' data.subset.1 <- damage_simulation[damage_simulation$building_cat == 1, ]
#'
#' frag.model.1 <- polr(CDF ~ logPGA, data = data.subset.1,
#'                      method = "probit", Hess = TRUE)
#' ex.prob.1 <- frag_curve(frag.model.1, data = data.subset.1)
#'
#' upper.bin <- c(0, 1, 10, 30, 60, 100)
#' bin.length <- c(1, 1, 10, 20, 30, 40, 1)
#' mDF_vector <- mean_DF(frag.model.1, data = data.subset.1, upper.bin,
#'                       bin.length, ex.prob = ex.prob.1, plot = TRUE)

mean_DF <- function(model, data, upper.bin, bin.length, ex.prob = NULL, plot = FALSE){

  # Function to compute the expected DF from a vector of damage state probabilities:

  mDF <- function(plot.bin.prob, upper.bin, bin.length){
    plot.y <- plot.bin.prob/bin.length
    weighted_ave <-  upper.bin[length(upper.bin)]*plot.bin.prob[length(plot.bin.prob)]
    for (i in 2:length(upper.bin)){
      x1 <- upper.bin[i-1]; x2 <- upper.bin[i]
      p_x1x2 <- plot.y[i]
      ave_contr <- (p_x1x2/2)*(x2^2 - x1^1)
      weighted_ave <- weighted_ave + ave_contr
    }
    return(weighted_ave)
  }

  # Compute the mean DF for dataset:
  damage.prob <- model$fitted.values

  mean.DF <- rep(NA, nrow(data))

  PGA.list.2 <- unique(data$PGA)
  for (i in 1:length(PGA.list.2)){ # Can do by data row if accounting for spatial corr in IM later.
    row.id <- which(data$PGA == PGA.list.2[i])
    if(length(row.id)>=1){
      mean.DF[row.id] <- mDF(damage.prob[row.id[1], ], upper.bin, bin.length)
    }
  }

  if(plot){

    PGA.list <- sort(unique(ex.prob$PGA), decreasing = FALSE)

    # Illustration of weighted average for CDF (mean CDF):
    mean.CDF.df <- data.frame("PGA" = PGA.list, "mean.CDF" = rep(NA, length(PGA.list)))
    plot.j <- ceiling(length(PGA.list)/2)

    CDF_val <- as.character(unique(data$CDF))
    plot.midpts <- as.numeric(CDF_val)

    for (j in 1:length(PGA.list)){
      plot.eg.data <- ex.prob[ex.prob$PGA == PGA.list[j], ]
      plot.bin.prob <- c(1, plot.eg.data$ex.prob) - c(plot.eg.data$ex.prob, 0)
      mean.CDF.df$mean.CDF[j] <- mDF(plot.bin.prob, upper.bin, bin.length)

      if (j == plot.j){
        plot.y <-plot.bin.prob/bin.length
      }

    }

    d1 <- data.frame(x=c(0, upper.bin), y=plot.y)
    d2 <- d1[c(1, nrow(d1)), ]

    p1 <- ggplot() + geom_step(data=d1, mapping=aes(x=x, y=y)) + geom_point(data=d2, mapping=aes(x=x, y=y), color="red") + ylab("Probability density") + xlab(paste("Damage factor (PGA = ", round(PGA.list[plot.j], 3), ")", sep = "")) + ggtitle("(a)") + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

    p2 <- ggplot() + geom_line(data = mean.CDF.df, aes(x=PGA, y=mean.CDF)) + ylab("Mean damage factor") + xlab("PGA") + ggtitle("(b)") + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

    grid.arrange(p1, p2, ncol = 2)

  }

  return(mean.DF)
}
