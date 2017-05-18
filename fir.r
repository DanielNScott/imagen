fir <- function(onsets, durations, activations, n_scans, n_conds, n_voxels) {

# Scale parameter - sets resolution
scale <- 10

# Rescale onsets
onsets    <- onsets    * scale
durations <- durations * scale

n_scans <- n_scans * scale
res     <- TR / scale

#if (type == "user") 
#   shrf <- sum(hrf(0:(ceiling(scans) - 1)/scale))
#   no  <- length(onsets)
#if (length(durations) == 1) {
#  durations <- rep(durations, no)
#}
#else if (length(durations) != no) {
#  stop("Length of duration vector does not match the number of onsets!")
#}
stims <- rep(0, ceiling(scans))

## ESTIMATE HRF USING FIR BASIS SET

# CREATE FIR DESIGN MATRIX
# WE ASSUME HRF IS 16 TRS LONG
hrf_len <- 16
 
# BASIS SET FOR EACH CONDITOIN IS A TRAIN OF INPULSES
fir_bases <- zeros(n_scans, hrf_len*n_conds)
 
for (cond in 1:n_conds) {
    col_subset <- ((cond - 1)* hrf_len + 1):(cond*hrf_len)
    
    for (onset in 1:numel(onsets[,cond]) ) {
        #impulse_times <- onsets(onset):onsets(onset) + hrf_len - 1;
        impulse_times <- seq(onsets(onset), onsets(onset) + hrf_len*scale - 1, scale)
        
        for (impulse in 1:numel(impulse_times)) {
            fir_bases(impulse_times(impulse), col_subset(impulse)) <- 1;
        }
    }
}

# ESTIMATE HRF FOR EACH CONDITION AND VOXEL
fir_hrf_est <- pseudoinverse(fir_bases %*% fir_bases) * fir_bases %*% activations
 
# RESHAPE HRFS
hHatFIR <- reshape(fir_hrf_est, hrf_len, n_conds, n_voxels)

}