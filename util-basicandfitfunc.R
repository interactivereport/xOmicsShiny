
# this file is from https://cran.r-project.org/web/packages/DRomics/


###########################################
# Basic functions defining models 
# for plot, fit, BMD calculation
###########################################

flin <- function(x, b, d)
  # x the dose
  # d the y value for a null dose
  # b the slope
{
  d  + b*x
}
## inverse lin (X for an y value)
invlin <- function(y, b, d)
{
  (y - d) /b
}

### Expo model and starting values
formExp3p <- as.formula(signal ~ d + b * (exp(dose/e)  - 1) )
startvalExp3pnls.1 <- function(xm, ym, increase, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1] 
  # initial value of absolute valeur of e
  eabs <- 0.1*max(xm)
  
  if ((increase & !Ushape) | (!increase & Ushape))
  {
    # initial value of c even if it does not appear as a parameter
    c <- ym[which.max(xm)] # mean of the response at the highest dose
    b <- d - c
    # initial value of e that corresponds to its value in the 2 parameter model
    e <- -eabs
  } else
  {
    # initial value of e that corresponds to its value in the 2 parameter model
    e <- eabs
    # initial value of b 
    reg <- lm(ym ~ exp(xm / e))
    b <- coef(reg)[2]
  }
  startval <- list(b = b, d = d, e = e)
}

startvalExp3pnls.2 <- function(xm, ym, increase, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1] 
  # initial value of absolute valeur of e
  eabs <- max(xm)     # WHAT CHANGES IN THE SECOND TRIAL !!!!!!!!!!!
  
  if ((increase & !Ushape) | (!increase & Ushape))
  {
    # initial value of c even if it does not appear as a parameter
    c <- ym[which.max(xm)] # mean of the response at the highest dose
    b <- d - c
    # initial value of e that corresponds to its value in the 2 parameter model
    e <- -eabs
  } else
  {
    # initial value of e that corresponds to its value in the 2 parameter model
    e <- eabs
    # initial value of b 
    reg <- lm(ym ~ exp(xm / e)) 
    b <- coef(reg)[2]
  }
  startval <- list(b = b, d = d, e = e)
}


# for plot
fExpo <- function(x, b, d, e)
{
  d + b * (exp(x/e)  - 1)
}
# x the dose
# d the y value for a null dose
# b a shape parameter 
# e another shape parameter 
# if e < 0 the model tends to a non infinite limit (d - b) when x tends toward infinity
# the model is increasing if e*b > 0

## inverse Expo (X for an y value)
invExpo <- function(y, b, d, e)
{
  if ( ((e < 0) & (b < 0) & (y > d - b)) | ((e < 0) & (b > 0) & (y < d - b)) )
    return(NaN) else
      return(e * log(1 + (y - d) / b))
}

### Hill model and starting values
formHill <- as.formula(signal ~ c + (d - c) / (1 + (dose/e)^b ) )
startvalHillnls2 <- function(x, y, xm, ym, increase) # requires the definition of increase from min and max values
  # inputs
  # - x values of the dose
  # - y values the corresponding signal
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # 
{
  maxi <- max(y, na.rm = TRUE)
  mini <- min(y, na.rm = TRUE)
  ampl <- maxi - mini
  # inflate maxi and mini so as all values are strictly inside the interval [mini; maxi]
  maxi <- maxi + 0.001 * ampl
  mini <- mini - 0.001 * ampl
  # initial value of c
  c <- ifelse(increase, maxi, mini) 
  # initial value of d
  d <-ifelse(increase, mini, maxi) 
  # initial value of e and b from regression
  yreg <- log((d - c) / (y[x!=0] - c) - 1)
  xreg <- log(x[x!=0])
  reg <- lm(yreg ~ xreg)
  b <- reg$coefficients[2]
  e <- reg$coefficients[1] / (-b)
  startval <- list(b = b, c = c, d = d, e = e)
}


# for plot
fHill <- function(x, b, c, d, e)
{
  c + (d - c) / (1 + (x/e)^b)
}

# x the dose
# c the y value for a high dose
# d the y value for a null dose
# b a shape parameter > 0
# e the dose at the inflexion point 

## inverse Hill (X for an y value)
invHill <- function(y, b, c, d, e)
{
  if ( ((d < c) & (y > c)) | ((d > c) & (y < c)) )
    return(NaN) else
      return(e * ((d - y) / (y - c))^(1/b))
}

### Gaussian model 5 p and starting values
formGauss5p <- as.formula(signal ~ f * exp(-0.5 * ((dose-e)/b)^2) + d + (c - d) * pnorm((dose-e)/b)) 
startvalGauss5pnls <- function(xm, ym, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1]# mean of the response at the first dose
  # initial value of c
  c <- ym[which.max(xm)]# mean of the response at the highest dose
  # initial value of f
  yextremum <- ifelse(Ushape, min(ym), max(ym))
  f <-  yextremum - (c + d) / 2 # amplitude of the gaussian part
  # initial value of b (standard error) assuming (after having looked many curves) the
  # bell shape extends to the whole dose range so dose max = 4 sd = 4 * b
  # 
  b <- max(xm) / 4
  # initial value of e (dose corresponding to the maximal (or minimal) signal)
  xextremum <- median(xm[which(ym == yextremum)]) # just in case there is more than one dose at which ym == yextremum
  e <- min(xextremum - (c - d)*b/(f*sqrt(2*pi)), 1e-6) 
  startval <- list(b = b, c = c, d = d, e = e, f = f)
}

## simplified version for cases where Gauss 5p does not converge (c = d)
formGauss4p <- as.formula(signal ~ f * exp(-0.5 * ((dose-e)/b)^2) + d ) # without c and without probit part
startvalGauss4pnls <- function(xm, ym, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[which.max(xm)]# mean of the response at the highest dose (supposed identical to the other asymptote)
  # initial value of f
  yextremum <- ifelse(Ushape, min(ym), max(ym))
  f <-  yextremum - d # amplitude of the gaussian part
  # initial value of b (standard error) assuming (after having looked many curves) the
  # bell shape extends to the whole dose range so dose max = 4 sd = 4 * b
  # with 2sd between 0 and e -> b = e /2
  b <- max(xm) / 4
  # initial value of e (dose corresponding to the maximal (or minimal) signal)
  xextremum <- median(xm[which(ym == yextremum)] )# just in case there is more than one dose with ym = yextremum
  e <- min(xextremum, 1e-6)
  startval <- list(b = b, d = d, e = e, f = f)
}

# for plot

fGauss5p <- function(x, b, c, d, e, f)
{
  f * exp(-0.5 * ((x-e)/b)^2) + # gaussian part
    d + (c - d) * pnorm((x-e)/b) # probit part
}

fGauss5pBMR <- function(x, b, c, d, e, g, threshold)
{
  g * exp(-0.5 * ((x-e)/b)^2) + # gaussian part
    d + (c - d) * pnorm((x-e)/b) - threshold # probit part
  
}

fGauss5pBMR_xinlog <- function(xinlog, b, c, d, e, g, threshold)
{
  x <- exp(xinlog)
  g * exp(-0.5 * ((x-e)/b)^2) + # gaussian part
    d + (c - d) * pnorm((x-e)/b) - threshold # probit part
  
}


### Gaussian model 5 p and starting values in log scale
formLGauss5p <- as.formula(signal ~ f * exp(-0.5 * (log(dose/e)/b)^2) + d + (c - d) * pnorm(log(dose/e)/b)) 
startvalLGauss5pnls <- function(xm, ym, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1]# mean of the response at the first dose
  # initial value of c
  c <- ym[which.max(xm)]# mean of the response at the highest dose
  # initial value of f
  yextremum <- ifelse(Ushape, min(ym), max(ym))
  f <-  yextremum - (c + d) / 2 # amplitude of the gaussian part
  # initial value of b (standard error) assuming (after having looked many curves) the
  # bell shape extends to the whole dose range so dose max = 4 sd = 4 * b IN LOG SCALE !!!
  b <- (log(max(xm)) - log(min(xm[xm!=0]))) / 4
  # initial value of e (dose corresponding to the maximal (or minimal) signal)
  xextremum <- median(xm[which(ym == yextremum)]) # median just in case there is more than one dose with ym = yextremum
  e <- xextremum * exp(- (c - d)*b/(f*sqrt(2*pi))) 
  startval <- list(b = b, c = c, d = d, e = e, f = f)
}


## simplified version for cases where Gauss 5p dose not converge (c = d)
formLGauss4p <- as.formula(signal ~ f * exp(-0.5 * (log(dose/e)/b)^2) + d) # without c and without probit part
startvalLGauss4pnls <- function(xm, ym, Ushape)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1]# mean of the response at the smallest dose
  # initial value of f
  yextremum <- ifelse(Ushape, min(ym), max(ym))
  f <-  yextremum - d # amplitude of the gaussian part
  # initial value of b (standard error) assuming (after having looked many curves) the
  # bell shape extends to the whole dose range so dose max = 4 sd = 4 * b IN LOG SCALE !!!
  b <- (log(max(xm)) - log(min(xm[xm!=0]))) / 4
  # initial value of e (dose corresponding to the maximal (or minimal) signal)
  xextremum <- median(xm[which(ym == yextremum)]) # just in case there is more than one dose with ym = yextremum
  e <- xextremum 
  startval <- list(b = b, d = d, e = e, f = f)
}


### log-probit model and starting values in log scale
formLprobit <- as.formula(signal ~ d + (c - d) * pnorm(log(dose/e)/b)) 
startvalLprobitnls1 <- function(xm, ym) # to suppress
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # - Ushape TRUE if U shape FALSE if umbrella shape
{
  # initial value of d
  d <- ym[1]# mean of the response at the first dose
  # initial value of c
  c <- ym[which.max(xm)]# mean of the response at the highest dose
  # initial value of b (standard error) assuming (after having looked many curves) the
  # bell shape extends to the whole dose range so dose max = 4 sd = 4 * b IN LOG SCALE !!!
  ldosemin <- log(min(xm[xm!=0]))
  ldosemax <- log(max(xm))
  b <- (ldosemax - ldosemin) / 4
  # initial value of e (in the middle of the dose range in log)
  e <- exp( ldosemin + (ldosemax - ldosemin) /2 )
  startval <- list(b = b, c = c, d = d, e = e)
}

startvalLprobitnls2 <- function(x, y, xm, ym, increase) # requires the definition of increase from min and max values
  # inputs
  # - x values of the dose
  # - y values the corresponding signal
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # 
{
  maxi <- max(y, na.rm = TRUE)
  mini <- min(y, na.rm = TRUE)
  ampl <- maxi - mini
  # inflate maxi and mini so as all values are strictly inside the interval [mini; maxi]
  maxi <- maxi + 0.001 * ampl
  mini <- mini - 0.001 * ampl
  # initial value of c
  c <- ifelse(increase, maxi, mini) 
  # initial value of d
  d <-ifelse(increase, mini, maxi) 
  # initial value of e and b from regression
  Y <- (y[x!=0] - d) / (c - d)
  yreg <- qnorm( Y )
  xreg <- log(x[x!=0])
  reg <- lm(yreg ~ xreg)
  b <- 1/ reg$coefficients[2]
  e <- reg$coefficients[1] * (-b)
  startval <- list(b = b, c = c, d = d, e = e)
}


# for plot

fLGauss5p <- function(x, b, c, d, e, f)
{
  f * exp(-0.5 * (log(x/e)/b)^2) + # gaussian part
    d + (c - d) * pnorm(log(x/e)/b) # probit part
}

fLGauss5pBMR <- function(x, b, c, d, e, g, threshold)
{
  g * exp(-0.5 * (log(x/e)/b)^2) + # gaussian part
    d + (c - d) * pnorm(log(x/e)/b) - threshold # probit part
  
}

fLGauss5pBMR_xinlog <- function(xinlog, b, c, d, e, g, threshold)
{
  g * exp(-0.5 * ((xinlog - log(e))/b)^2) + # gaussian part
    d + (c - d) * pnorm((xinlog - log(e))/b) - threshold # probit part
  
}


fLprobit <- function(x, b, c, d, e)
{
  d + (c - d) * pnorm(log(x/e)/b)
}

## inverse Lprobit (X for an y value)
invLprobit <- function(y, b, c, d, e)
{
  if ( ((d < c) & (y > c)) | ((d > c) & (y < c)) )
    return(NaN) else
      return(e * exp(qnorm((y - d) / (c - d)) *b))
}

