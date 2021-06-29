data <- dynr.data(out, id="id", time="time", observed=c("y1","y2"))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.measurement(
  values.load=diag(1, 2),
  obs.names = c('y1', 'y2'),
  state.names=c('x1', 'x2'))

# Initial conditions on the latent state and covariance
initial <- prep.initial(
  values.inistate=c(.2,.5),
  params.inistate=c("fixed", "fixed"),
  values.inicov=diag(c(1e-4,1e-4)), 
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
formula=list(list(x1~r1/K1*x1*(K1-x1+alpha12*x2),
                  x2~r2/K2*x2*(K2-x2+alpha21*x1)))
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r1 = 2, r2 = 2,K1 = 5, K2=5,
                                      alpha12 = 1,alpha21 = 2),
                           isContinuousTime=TRUE)
# (optional)
# Parameter transformations 
# (Exponential transformation will ensure positive parameters)
trans<-prep.tfun(formula.trans=list(r1~exp(r1), 
                                    r2~exp(r2),
                                    K1~exp(K1),
                                    K2~exp(K2),
                                    alpha12~exp(alpha12),
                                    alpha21~exp(alpha21)),
                 formula.inv=list(r1~log(r1),
                                  r2~log(r2),
                                  K1~log(K1),
                                  K2~log(K2),
                                  alpha12~log(alpha12),
                                  alpha21~log(alpha21)))

#----  Put all the recipes together into an overall model  ----
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial,
                    transform=trans, 
                    data=data,
                    outfile="Comp.c")

# (optional)
# If we know a range of possible parameter values, we can help the optimization
# algorithm for parameter estimation by bounding the search area
model@ub <- rep(5,8)
model@lb <- rep(-5,8)


#---- Model fitting for the maximum likelihood parameter estimates  ----

# Sometimes the optimization algorithm gets stuck
# Changing the starting values for the parameters may help
repeat{
  model@xstart <- runif(8,-2,2)
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


# For the competition model, change the formula argument to:
formula=list(list(x1~r1/K1*x1*(K1-x1-alpha12*x2),
                  x2~r2/K2*x2*(K2-x2-alpha21*x1)))
