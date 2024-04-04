library(pacman)

# Load necessary packages
p_load(tidyverse, car, ggplot2, sf, glmulti, spdep, sp, INLA, spatialreg, tmap)

# Load the data
data <- read_csv("../../data/data_round5_with_covariates.csv")

map <- st_read("../../data/shapefiles/STzVill_NBS2012_data_covariates_round5.shp")

# Scale the data
s_data <- map %>% 
  mutate_at(vars(pvrty_n, elevatn, pop_den, Ar_kmsq, dist_km), scale) %>%
  mutate_at(vars(pvrty_n), exp)

# Fit a binomial regression model
model2 <- glm(cbind(collars, nocllrs) ~ nwctgry + Ar_kmsq + pvrty_n + elevatn,
              data = s_data, family = "binomial",
              na.action = na.exclude)

# View model summary
summary(model2)


# Calculate Variance Inflation Factors (VIF)
vif(model2)

# Fit multivariable model using glmulti function
formula <- cbind(collars, nocllrs) ~ nwctgry + Ar_kmsq + pvrty_n + elevatn + pop_den + dist_km
res <- glmulti(formula, data = s_data, method = "h",
               level = 1, crit = "aicc", family = binomial)

# Print result
print(res)
plot(res)

# Select top models
top <- weightable(res)
top <- top[top$aicc <= min(top$aicc) + 2, ]
top



# Checking for spatial correlation
s_data$res <- as.numeric(residuals(model2))
nb <- poly2nb(s_data)
wts <- nb2listw(neighbours = nb, style = 'B', zero.policy = T)
mor.i <- moran.test(x = s_data$res, listw = wts,  alternative = "greater",
                    zero.policy = T, na.action = na.exclude)
mor.i



######################### Fitting non-spatial model
formula <- collars ~  nwctgry + Ar_kmsq + pvrty_n  + elevatn + dist_km + dist_km*nwctgry
res0 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link=1), verbose = T,
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))


summary(res0)


##### non-spatial structure but with independent random effect
formula <- collars ~  nwctgry + Ar_kmsq + pvrty_n  + elevatn + dist_km + dist_km*nwctgry + f(ID, model = "iid") 
res1 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link=1), verbose = T, 
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                    return.marginals.predictor=TRUE))


summary(res1)

# Fit spatial model (Besags)
formula <- collars ~ nwctgry + Ar_kmsq + pvrty_n + elevatn + dist_km +
  f(ID, model = "besag", graph = g) 
res2 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link = 1), verbose = T, 
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                    return.marginals.predictor = TRUE))
summary(res2)

######## spatial model bym
formula <- collars ~  nwctgry + Ar_kmsq + pvrty_n  + elevatn + dist_km + dist_km*nwctgry + 
  f(ID, model = "bym", graph = g
  ) 
res3 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link=1), verbose = T, 
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                    return.marginals.predictor=TRUE),
             control.family = list(control.link=list(model="logit")))
summary(res3)

# Fit spatial model (BYM2)

formula <- collars ~ nwctgry + Ar_kmsq + pvrty_n + elevatn + dist_km + dist_km * nwctgry + 
  f(ID, model = "bym2", graph = g,  hyper = prior) 
res4 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link=1), verbose = T, 
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                    return.marginals.predictor=TRUE),
             control.family = list(control.link=list(model="logit")))
summary(res4)

####################### Leroux
Q <- Diagonal(x = sapply(nb, length))
Q <- Q - spatialreg::as_dsTMatrix_listw(nb2listw(nb, style = "B", zero.policy = T))
C <- Diagonal(x = 1, n = nrow(s_data)) - Q
formula <- collars ~  nwctgry + Ar_kmsq + pvrty_n + elevatn + dist_km + dist_km*nwctgry + f(ID, model = "generic1", Cmatrix = C)
res5 <- inla(formula = formula, 
             data = s_data, Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE,link=1), verbose = T,
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                    return.marginals.predictor=TRUE))


summary(res5)



# Fit geostatistical model
coo <- s_data %>% 
  st_centroid() %>% 
  st_geometry() %>% 
  st_coordinates()
max.edge <- as.numeric(diff(range(coo[,1])) / (3 * 5))
bound.outer <- as.numeric(diff(range(coo[,1])) / 3)

mesh <- inla.mesh.2d(loc = coo,
                     max.edge = c(1, 2) * max.edge,
                     offset = c(max.edge, bound.outer),
                     cutoff = max.edge / 5)
plot(mesh)
points(coo, col = "red")

spde <- inla.spde2.matern(mesh = mesh, alpha = 3/2)
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 3/2, constr = TRUE,
                            prior.range = c(100000, 0.5),
                            prior.sigma = c(1, 0.01))

indexs <- inla.spde.make.index("s", spde$n.spde)
A <- inla.spde.make.A(mesh = mesh, loc = coo)

stk.e <- inla.stack(
  tag = "est",
  data = list(y = s_data$collars),
  A = list(1, A),
  effects = list(data.frame(b0 = 1, data.frame(s_data %>% st_set_geometry(NULL))), s = indexs)
)

formula <- collars ~ nwctgry + Ar_kmsq + pvrty_n + elevatn + dist_km + dist_km * nwctgry + 
  f(s, model = spde) + f(ID, model = "iid")
res6 <- inla(formula = formula, 
             data = inla.stack.data(stk.e), Ntrials = tdogs,
             family = "binomial",
             control.predictor = list(compute = TRUE, link = 1), verbose = T, 
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                    return.marginals.predictor = TRUE),
             control.family = list(control.link = list(model = "logit")))
summary(res6)

# Output results
res.field <- inla.spde2.result(res6, name = "s", spde, do.transf = TRUE)
res.field$summary.hyperpar
inla.emarginal(function(x) x, res.field$marginals.kappa[[1]])
inla.emarginal(function(x) x, res.field$marginals.variance.nominal[[1]])
inla.emarginal(function(x) x, res.field$marginals.range.nominal[[1]])
