gen_addnl_flds <- function(data) {

  # TCI impulsivity really looks like it falls into categories 'general' and 'financial' to me...
  # Also, it looks like taking rowSums is fine -- the summary vars get NAs if there are NAs in row.
  tci_gen_imp <- c('tci010', 'tci014', 'tci047', 'tci071', 'tci102', 'tci123', 'tci193', 'tci210')
  tci_fin_imp <- c('tci024', 'tci059', 'tci105', 'tci215')

  data$raw['tci_gen_imp'] <- rowSums(data$raw[, tci_gen_imp])
  data$raw['tci_fin_imp'] <- rowSums(data$raw[, tci_fin_imp])
  data$names$TCI <- c(data$names$TCI, 'tci_gen_imp', 'tci_fin_imp')

  #SURPS impulsivity just needs to be summed.
  # Question 22 about being manipulative is not included here.
  surps_imp   <- c('surps5', 'surps2', 'surps11', 'surps15')

  data$raw['surps_imp']   <- rowSums(data$raw[, surps_imp])
  data$names$SURPS <- c(data$names$SURPS, 'surps_imp')

  # Need some ESPAD summaries and other stuff...
  # Individual fields:
  # C.18i    - likelihood of regret
  # C.prev31 - scaled num drinks on typical day in which drinking
  # C.21     - scaled num drinks needed to get drunk
  # C.19b    - scaled num times drunk in last year
  # C.6      - scaled num times smoked in life -- should match 'espad_6_life_nic'
  #
  data$raw['binge'] <- data$raw[,'C.prev31'] / data$raw[,'bmi']
  data$names$ESPAD  <- c(data$names$ESPAD, 'binge')


  return(data)
}

