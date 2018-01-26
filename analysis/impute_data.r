# ------------------------------------------------------------------------------ #
#                   This function controls the data imputation.                  #
# ------------------------------------------------------------------------------ #
impute_data <- function(data, n_imputed = 20, max_iters = 500) {

    # Multiple Imputation by Chained Equations (MICE)
    library(Rcpp)
    library(mice)

    to_impute <- c('sj1a', 'sj1b', 'sj1c', 'tci_gen_imp', 'tci_fin_imp',
                   'surps_imp', 'log10.k.', 'binge', data$names$Age, data$names$CANTAB,
                   'espad_6_life_nic', 'espad_8a_alc_life', data$names$genes,
                   data$names$sst, data$names$mid)

    imp <- mice(data$raw[to_impute][data$train_inds,], print = FALSE)
    imp <- mice(data$raw[to_impute][data$train_inds,],
                pred = imp$predictorMatrix, seed = 23109, m = n_imputed,
                maxit = max_iters, print = FALSE)

    data$imp <- imp
    data$imp_names <- to_impute

    return(data)
}
# ------------------------------------------------------------------------------ #