# List of potential ancestral graphs for consideration
define_ag_models <- function() {
  # classic hyper-indirect model (hindi)
  ag0 <- makeMG(dg = DAG(
          rCaudate  ~ rPreSMA + rIFG,
          rSTN      ~ rPreSMA + rIFG,
          rGPe      ~ rCaudate,
          rGPi      ~ rGPe  +  rSTN,
          rThalamus ~ rGPi
          ),ug = UG(  ~ rPreSMA*rIFG))

  # hyperdirect (hyp)
  ag1 <- makeMG(dg = DAG(
         #rCaudate  ~ rPreSMA + rIFG,
          rSTN      ~ rPreSMA + rIFG,
         #rGPe      ~ rCaudate,
          rGPi      ~ rSTN ,
          rThalamus ~ rGPi
          ),ug = UG( ~ rPreSMA*rIFG + rCaudate + rGPe ))

  # indirect (indir)
  ag2 <- makeMG(dg = DAG(
          rCaudate  ~ rPreSMA + rIFG,
         #rSTN      ~ rPreSMA + rIFG,
          rGPe      ~ rCaudate,
          rGPi      ~ rGPe ,
          rThalamus ~ rGPi
          ),ug = UG(  ~ rPreSMA*rIFG + rSTN ))

  # direct (dir)
  ag3 <- makeMG(dg = DAG(
          rCaudate  ~ rPreSMA  +  rIFG,
          #rSTN     ~ rPreSMA + rIFG,
          #rGPe     ~ rCaudate,
          rGPi      ~ rCaudate,
          rThalamus ~ rGPi
          ),ug = UG( ~ rPreSMA*rIFG + rSTN  + rGPe ))

  # classic hyper-indirect model (hindi)
  nostnhindi <- makeMG(dg = DAG(
          rCaudate  ~ rPreSMA + rIFG,
          rGPe      ~ rCaudate,
          rGPi      ~ rGPe  +  rPreSMA + rIFG,
          rThalamus ~ rGPi
          ),ug = UG(  ~ rPreSMA*rIFG))

  # direct (dir)
  nostndir <- makeMG(dg = DAG(
          rCaudate  ~ rPreSMA  +  rIFG,
          #rSTN     ~ rPreSMA + rIFG,
          #rGPe     ~ rCaudate,
          rGPi      ~ rCaudate,
          rThalamus ~ rGPi
          ),ug = UG( ~ rPreSMA*rIFG + rPreSMA + rIFG + rGPe ))

  # Package
  #models <- list(ag0, ag1, ag2, ag3)
  #names(models) <- c('hindi', 'hyp', 'indir', 'dir')
  models <- list(nostnhindi, nostndir)
  names(models) <- c('hindi', 'dir')

  return(models)
}
