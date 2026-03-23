## -----------------------------------------------------------------------------
# xy_train <- read.csv("P:/BMS/pacakge_data/samples/train_xy.csv")
# obs_train <- read.csv("P:/BMS/pacakge_data/samples/train_obs.csv")
# pred_train <- read.csv("P:/BMS/pacakge_data/samples/train_pred.csv")
# 
# full_data <- cbind(xy_train, obs_train, pred_train)


## ----warning=FALSE, message=FALSE, fig.width=5, fig.height=4------------------
# ref_density <- ClassicHCAS::ref_density(
#   data = full_data,
#   radius_km = 1000,
#   bin_width = 0.03
# )
# 
# plot(ref_density)


## ----fig.width=5, fig.height=4------------------------------------------------
# norm <- ClassicHCAS::normalise(x = ref_density)
# 
# plot(norm)


## -----------------------------------------------------------------------------



