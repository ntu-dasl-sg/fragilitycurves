#' Compute fitted exceedance probabilities and plot fragility curves.
#'
#' @param model An ordinal model fitted for \code{data} using \code{polr} from the R package \code{MASS} or a list of probit regression models fitted for the different damage states in \code{data} using \code{glm}. E.g. Grade 0 against Grades 1-5, Grades 0-1 against Grades 2-5 etc.
#' @param type "ordinal" (default) or "probit" if supplying a list of probit regressions for \code{model}.
#' @param data A dataframe with the columns \code{CDF} (ordered factor), \code{logPGA} (numeric), \code{PGA} (numeric), \code{Easting} (numeric) and \code{Northing} (numeric).
#' @param plot A logical value indicating if the fragility curves should be plotted.
#' @return A dataframe with the columns \code{PGA}, \code{CDF} and \code{ex.prob} (the exceedance probability of the \code{CDF} level at the \code{PGA} value).
#' @import ggplot2
#' @importFrom stats coef pnorm
#' @importFrom gridExtra grid.arrange
#' @export
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

frag_curve <- function(model, type = "ordinal", data, plot = TRUE){

  master.PGA.list <- sort(unique(data$PGA), decreasing = FALSE)

  damage.df<- data.frame(PGA = master.PGA.list)

  damage.states <- levels(data$CDF)

  if(type == "probit" & length(model)!=(length(damage.states)-1)){print("Error: Number of models in list does not match with the number of damage states.")}

  temp.matrix <- rep(NA, length(damage.states)*3)
  for (i in 1:length(master.PGA.list)){
    PGA.states <- data$CDF[data$PGA == master.PGA.list[i]]
    temp.n <- rep(NA, length(damage.states))
    temp.x <- rep(NA, length(damage.states))
    for (j in 1:length(damage.states)){
      temp.n[j] <- length(PGA.states)
      temp.x[j] <-sum(PGA.states >= damage.states[j])
    }
    temp.p <- temp.x/temp.n
    temp.matrix <- rbind(temp.matrix, c(temp.n, temp.x, temp.p))
  }

  temp.matrix <- temp.matrix[-1, ] # Remove dummy first row.
  damage.df <- cbind(damage.df, temp.matrix)
  colnames(damage.df) <- c("PGA", paste(rep(c("n", "x", "p"), each = length(damage.states)), damage.states, sep = "."))


  p.names <- paste("p", damage.states, sep = ".")
  n.names <- paste("n", damage.states, sep = ".")
  if(type == "ordinal"){
    res.table <- coef(summary(model))
    res.names <- row.names(res.table)
  }

  # Plot empirical exceedance probabilities of damage states (to fit ordinal fragility curves to):

  p.list <- list()
  p.list[[1]] <- ggplot(damage.df, aes(x=PGA, y=p.0)) + geom_point(shape = 1, aes(size = n.0)) + labs(size = "No. of buildings", y = "P(CDF >= 0)") + ylim(c(0, 1))

  ex.prob <- numeric()
  for (i in 2:length(damage.states)){
    temp.df <- damage.df[, c("PGA", p.names[i], n.names[i])]; colnames(temp.df) <- c("PGA", "p", "n")
    p.list[[i]] <- ggplot(temp.df, aes(x=PGA, y=p)) + geom_point(shape = 1, aes(size = n)) + labs(size = "No. of buildings", y = paste("P(CDF >= ", damage.states[i], ")", sep = "")) + ylim(c(0, 1)) + xlim(c(0, max(damage.df$PGA)))
    if(type == "ordinal"){
      curve.fn <- function(x){pnorm(res.table[res.names[i], "Value"] - res.table["logPGA", "Value"]*log(x + 0.001), lower.tail = FALSE)}
    }else{
      curve.fn <- function(x){pnorm(-model[[i-1]]$coefficients[1] - model[[i-1]]$coefficients[2]*log(x + 0.001), lower.tail = FALSE)}
    }
    ex.prob <- c(ex.prob, curve.fn(seq(0, max(master.PGA.list), length.out = 200)))
  }

  ex.prob.df <- data.frame("PGA" = rep(seq(0, max(master.PGA.list), length.out = 200), length(damage.states[-1])), "CDF" = rep(damage.states[-1], each = 200))
  ex.prob.df <- cbind(ex.prob.df, ex.prob)
  # Change CDF to ordinal damage scale:
  ex.prob.df$CDF <- ordered(ex.prob.df$CDF, levels = damage.states)

  p.list[[length(damage.states)+1]] <- ggplot(data = ex.prob.df, aes(x = PGA, y = ex.prob)) + geom_line(aes(lty = CDF)) + labs(y = "P(CDF >= DS)", lty = "Damage State (DS)") + ylim(c(0, 1)) + xlim(c(0, max(damage.df$PGA)))

  for (i in 2:length(damage.states)){
    p.list[[i]]  <- p.list[[i]] + geom_line(data = ex.prob.df[ex.prob.df$CDF == damage.states[i], ],aes(x = PGA, y = ex.prob))
  }

  if (plot){
    do.call(grid.arrange, c(p.list, ncol = 2))
  }

  return(ex.prob.df)
}
