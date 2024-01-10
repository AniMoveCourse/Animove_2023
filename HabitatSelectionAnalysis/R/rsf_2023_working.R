#' # Resource selection functions
#' Björn Reineking, 2023-08-05
#' 
#' Model the relative density of animals (also called range distribution or utilisation distribution) as a function of environmental predictors. 
#' 
# We will use the buffalo data set.
# Install animove metapackage
# https://github.com/AniMoveCourse/animove_R_package
# install.packages("remotes")
# remotes::install_github("AniMoveCourse/animove_R_package")
# 
# This code needs ctmm version 1.1.1 (or newer)
# detach("package:ctmm", unload = TRUE)
# remotes::install_github("ctmm-initiative/ctmm")

#'Loading packages
#+  results='hide', message=FALSE, warning=FALSE
library(animove)
library(ctmm)
library(sf)
library(mgcv)
library(mvtnorm)
library(lubridate)

#' ## Load buffalo data
#' See Kami's slides for how we got here...
data(buffalo_utm)

#' ## Environmental data: topography, waterways, and NDVI
data(buffalo_env)
raster::plot(buffalo_env)

#' ## Animals and elevation
raster::plot(raster(buffalo_env, 1))
points(buffalo_utm)

#' To speed up analyses, we will only work with one individual, Cilla
cilla <- buffalo_utm[["Cilla"]]
raster::plot(raster(buffalo_env, 1), ext = extent(cilla) * 2)
lines(cilla)

#' The first day of Cilla shows some unrealistic movements. We remove those observations
cilla <- cilla[timestamps(cilla) > min(timestamps(cilla)) + days(1), ]

raster::plot(raster(buffalo_env, 1), ext = extent(cilla) * 2)
lines(cilla)

#' # Minimal example of rsf.fit

#' Create telemetry object - a ctmm object
cilla_telemetry <- as.telemetry(cilla)

#' Fit ctmm model
cilla_guess <- ctmm.guess(cilla_telemetry, CTMM=ctmm(isotropic = TRUE), interactive = FALSE)
cilla_select <- ctmm.select(cilla_telemetry, cilla_guess)

#' Fit akde
cilla_akde <- akde(cilla_telemetry, cilla_select)
plot(cilla_akde)

#' Create named list of rasters
be <- list("elev" = raster(buffalo_env, "elev"),
           "slope" = raster(buffalo_env, "slope"),
           "var_NDVI" = raster(buffalo_env, "var_NDVI"))

#' The integrator = "Riemann" option is still in testing, we use it here because it is much faster
cilla_rsf_riemann <- rsf.fit(cilla_telemetry, cilla_akde, R = be, integrator = "Riemann")

summary(cilla_rsf_riemann)

#' A suitability map
suitability_riemann <- ctmm::suitability(cilla_rsf_riemann, R=be, grid=crop(be[[1]], extent(cilla) * 2))
raster::plot(suitability_riemann)

#' Range distribution (includes the ranging behaviour)
my_grid <- crop(be[[1]], extent(cilla) * 2)
agde_cilla <- agde(CTMM = cilla_rsf_riemann, R = be, grid = my_grid)
plot(agde_cilla)

# Plot raster of range distribution
agde_raster <- ctmm::raster(agde_cilla, DF = "PMF")
plot(agde_raster)


#' # Traditional RSF with downweighted Poisson regression
#' Functions to generate quadrature points and predict with the model
rsf_points <- function(x, UD, R = NULL, n = 1e5, k = 1e6, type = "Riemann", 
                       rmax =6*sqrt(UD@CTMM$sigma[1,1]),
                       interpolation = FALSE) {
  # Samples background points from a 2D normal distribution fitted to the relocation data, and extracts environmental
  # information from a raster object "R"
  # x: telemetry object
  # UD: UD object
  # R: raster* object
  # n: number of background points to sample
  # k: weight of presence points
  # rmax: maximum distance for Riemann-type integration
  # interpolation: do interpolation when sampling the grid
  # When type 0´= "MonteCarlo", importance sampling is done
  stopifnot(UD@CTMM$isotropic)
  stopifnot(type %in% c("Riemann", "MonteCarlo"))
  if (type == "Riemann") {
    quadrature_pts <- getValues(R)
    xy <- as.data.frame(xyFromCell(R, 1:ncell(R)))
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(R))
    xy <- st_transform(xy, CRS(UD@info$projection))
    xy <- st_coordinates(xy)
    r <- sqrt(((xy[,1] - UD@CTMM$mu[1]))^2 + ((xy[,2] - UD@CTMM$mu[2]))^2)
    bg <- data.frame(case_ = 0,
                     x_ = xy[r<rmax,1], y_ = xy[r<rmax,2],  w_ = prod(res(R)), k_ = k)
    bg <- cbind(bg, quadrature_pts[r<rmax,])
    bg <- sf::st_as_sf(bg, coords = c("x_", "y_"), crs = CRS(UD@info$projection))
    xx <- data.frame(case_ = 1, x_ = x$x, y_ = x$y,
                     w_ = 1/k * UD$weights * mean(UD$DOF.area), k_ = k)
    xx <- sf::st_as_sf(xx, coords = c("x_", "y_"), crs = CRS(UD@info$projection))
    xx[names(R)] <- as.data.frame(raster::extract(R, st_transform(xx, crs(R)), method = ifelse(interpolation, "bilinear", "simple")))
    xx <- rbind(bg, xx)
  } else {
    quadrature_pts <- MASS::mvrnorm(n, mu = UD@CTMM$mu, Sigma = UD@CTMM$sigma)
    xx <- data.frame(case_ = 0, x_ = quadrature_pts[, 1], y_ = quadrature_pts[, 2], w_ = UD@CTMM$sigma[1,1]/n, k_ = k)
    xx <- rbind(xx, data.frame(case_ = 1, x_ = x$x, y_ = x$y,
                               w_ = 1/k * UD$weights * mean(UD$DOF.area), k_ = k 
    ))
    xx <- sf::st_as_sf(xx, coords = c("x_", "y_"), crs = CRS(UD@info$projection))
    xx[names(R)] <- as.data.frame(raster::extract(R, st_transform(xx, crs(R)), method = ifelse(interpolation, "bilinear", "simple")))
  }
  xy <- st_coordinates(xx)
  colnames(xy) <- c("x_", "y_")
  sd <- sqrt(UD@CTMM$sigma[1,1])
  xy[,1] <- (xy[,1] - UD@CTMM$mu[1])/sd
  xy[,2] <- (xy[,2] - UD@CTMM$mu[2])/sd
  xx <- cbind(xy, xx)
  xx
}

predict_rsf <- function(model, UD, object, include_avail = TRUE, data_crs = UD@info$projection) {
  if(include_avail) {
    xy <- as.data.frame(xyFromCell(object, 1:ncell(object)))
    xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs(object))
    xy <- st_transform(xy, data_crs)
    xy <- st_coordinates(xy)
    colnames(xy) <- c("x_", "y_")
    
    sd <- sqrt(UD@CTMM$sigma[1,1])
    xy[,1] <- (xy[,1] - UD@CTMM$mu[1])/sd
    xy[,2] <- (xy[,2] - UD@CTMM$mu[2])/sd
    
  } else {
    xy <- cbind("x_" = rep(0, ncell(object)), "y_" = rep(0, ncell(object)))
  }
  newdata <- as.data.frame(cbind(xy, getValues(object)))
  lambda <- as.numeric(predict(model, newdata, type = "link"))
  r <- raster(object, 1)
  r[] <- exp(lambda - max(lambda, na.rm = TRUE))
  r <- r / sum(getValues(r), na.rm = TRUE)
  r
}

#' ## A minimal "classic" example
#' Generate quadrature points ("background points")
set.seed(2)
rsf_cilla_df <- rsf_points(cilla_telemetry, cilla_akde, buffalo_env, interpolation = TRUE)
rsf_cilla_df <- rsf_cilla_df[!is.na(rsf_cilla_df$slope),]
#' Fit a downweighted Poisson regression
m_rsf_cilla <- glm(case_*k_ ~ x_ + y_ + I(-(x_^2 + y_^2)/2) + elev + slope + var_NDVI, 
                  family = poisson(), data= rsf_cilla_df, weights = w_)
#' Summary of model and confidence intervals for parameter estimates
#+ message=FALSE, warning=FALSE
summary(m_rsf_cilla)
confint(m_rsf_cilla)

#' Map of suitability, including the home ranging behaviour
suitability_glm <- predict_rsf(m_rsf_cilla, cilla_akde, crop(buffalo_env, extent(cilla) * 2))
raster::plot(suitability_glm)

#' Map of suitability, without the home ranging behaviour
suitability_no_avail_glm <- predict_rsf(m_rsf_cilla, cilla_akde, crop(buffalo_env, extent(cilla) * 2), include_avail = FALSE)
raster::plot(suitability_no_avail_glm)

#' Comparison of estimates from the two approaches (ctmm and "classic" via glm)
vars <- c("elev", "slope", "var_NDVI")

# Estimates are not identical but similar
cilla_rsf_riemann$beta[vars]
coef(m_rsf_cilla)[vars]

# Confidence intervals are of similar width, but location is different
#+ message=FALSE, warning=FALSE
ci_glm <- confint(m_rsf_cilla)
ci_glm[vars, ]
summary(cilla_rsf_riemann)$CI[sapply(vars, function(x) grep(x, rownames(summary(cilla_rsf_riemann)$CI))), ]

#' Create maps of suitability based on parameter estimates
M <- getValues(buffalo_env)
eta_glm <- M[,vars] %*% coef(m_rsf_cilla)[vars]
eta_rsf <- M[,vars] %*% cilla_rsf_riemann$beta[vars]

suit_r_glm <- raster(buffalo_env, 1)
suit_r_glm[] <- exp(eta_glm)
suit_r_rsf <- raster(buffalo_env, 1)
suit_r_rsf[] <- exp(eta_rsf)

suit_raster <- raster::stack(suit_r_rsf, suit_r_glm)
names(suit_raster) <- c("ctmm", "classic")
raster::plot(suit_raster)
raster::plot(suit_raster, ext = extent(cilla) * 2)

#' ## Multiple animals
#' 
#' - with rsf.fit: you can use the mean function on a list of rsf.fit objects
#' - "classic" approach: use glmmTMB and a mixed-effects Poisson regression: Muff, Signer & Fieberg (2020) J Anim Ecol 89: 80-92.
#' 

