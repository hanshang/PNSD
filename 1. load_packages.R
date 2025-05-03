################
# name packages
################

packages <- c("ftsa", "dplyr", "doMC", "tidyr", "ggplot2") 

## Now load or install and load all 

package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#################
# interval score
#################

# holdout: holdout data
# lb: lower bound
# ub: upper bound
# alpha: level of significance

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(cover, cpd, mean(score)))
}
