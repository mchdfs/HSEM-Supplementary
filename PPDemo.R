
data <- dynr.data(out, id="id", time="time", observed=c("y1","y2"))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.measurement(
  values.load=diag(1, 2),
  obs.names = c('y1', 'y2'),
  state.names=c('prey', 'predator'))

# Initial conditions on the latent state and covariance
initial <- prep.initial(
  values.inistate=c(5,3),
  params.inistate=c("fixed", "fixed"),
  values.inicov=diag(c(0.25,0.25)), 
  params.inicov=diag(c("fixed","fixed"))
)

# Measurement and dynamics errors
mdcov <- prep.noise(
  values.latent=diag(0, 2),
  params.latent=diag(c("fixed","fixed"), 2),
  values.observed=diag(rep(1,2)),
  params.observed=diag(c("var_1","var_2"),2)
)

# Model dynamics
formula=list(list(prey~ r*prey - c*prey*predator,
                  predator~ -m*predator + b*prey*predator))
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r = 2, c = 2, b = 1, m = 2),
                           isContinuousTime=TRUE)

# (optional)
# Parameter transformations 
# (Exponential transformation will ensure positive parameters)
trans<-prep.tfun(formula.trans=list(r~exp(r), 
                                    c~exp(c),
                                    b~exp(b),
                                    m~exp(m)),
                 formula.inv=list(r~log(r),
                                  c~log(c),
                                  b~log(b),
                                  m~log(m)))

#----  Put all the recipes together into an overall model  ----
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial,
                    transform=trans, 
                    data=data,
                    outfile="PP.c")

# (optional)
# If we know a range of possible parameter values, we can help the optimization
# algorithm for parameter estimation by bounding the search area
model@ub <- rep(1,6)
model@lb <- rep(-5,6)


#---- Model fitting for the maximum likelihood parameter estimates  ----


# Sometimes the optimization algorithm gets stuck
# Changing the starting values for the parameters may help
repeat{
  model@xstart <- runif(6,-1,1)
  # Estimate free parameters
  res <- try(dynr.cook(dynrModel=model),silent = TRUE)
  # This checks whether 1. the model converged and 
  # 2. whether the optimizer is stuck on a boundary (which we do not want)
  if(class(res)!="try-error"&&sum(res@bad.standard.errors)==0&&
     all(round(res@fitted.parameters,6)!=model@ub)&&
     all(round(res@fitted.parameters,6)!=model@lb)){break}
}

# Examine results
summary(res)
coef(res)
