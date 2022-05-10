# LHS model to sample from ILI parameter ranges
# Tanner Shea and Julie Spencer
require(lhs)

# Set seed in case you need to re-generate same sample set, change seed if you want different sample set
set.seed(3333)

# Generate nsamples for nparam different parameters (should be 10000)
nsamples <- 10000
nparam <- 5
lhs_raw <- randomLHS(nsamples,nparam) 
# Can switch to maximinLHS(nsamples,nparam) for less correlation between parameter sets but takes more time

## Transform raw sampling values into uniform distributions for each parameter range

# Set minimum and maximum points for each parameter
beta.min <- .00002
beta.max <- .0002
gamma1.min <- 1/14.7
gamma1.max <- 1/1.9
gamma2.min <- 1/8.56
gamma2.max <- 1/8.1
gamma3.min <- 1/11
gamma3.max <- 1/1.5
gamma4.min <- 1/11
gamma4.max <- 1/1.5
# p1.min <- 0
# p1.max <- 1
# p2.min <- 0
# p2.max <- 1

lhs_xform <- matrix(nrow = nrow(lhs_raw), ncol = ncol(lhs_raw))
lhs_xform[,1] <- qunif(lhs_raw[,1], min = beta.min, max = beta.max)
lhs_xform[,2] <- qunif(lhs_raw[,2], min = gamma1.min, max = gamma1.max)
lhs_xform[,3] <- qunif(lhs_raw[,3], min = gamma2.min, max = gamma2.max)
lhs_xform[,4] <- qunif(lhs_raw[,4], min = gamma3.min, max = gamma3.max)
lhs_xform[,5] <- qunif(lhs_raw[,5], min = gamma4.min, max = gamma4.max)
#lhs_xform[,6] <- qunif(lhs_raw[,6], min = p1.min, max = p1.max)
#lhs_xform[,7] <- qunif(lhs_raw[,7], min = p2.min, max = p2.max)

lhs_xform

write.csv(lhs_xform, file = "LHS_samples.csv",row.names=FALSE)
