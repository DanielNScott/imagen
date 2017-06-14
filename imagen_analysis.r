
data <- pre_process(re_read = TRUE, shuffle = FALSE, n_imputed = 30, max_iters = 1000)

analysis <- get_analysis()

plot_all(data, analysis)