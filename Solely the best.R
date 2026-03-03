





#AFI34806: Modelling Marine Socio-ecological systems

#title: Group work: North Sea Flatfish
#group 3
#authors: Eden Baas, Eline van Gelderen & Manu Krishnadath 
#date: 03/02/2026

### 1) Load packages
library(tidyverse)




## 2.1) Add length from age to the model with the von Bertalanffy Growth curve

Lasym <- 54.5   #asympotic length (cm) for plaice
k     <- 0.1   #growth coefficient
t0    <- 0     #theoretical age at zero length

a <- 0.0089   # scaling coefficient
b <- 3.04     # scaling exponent

# Length-at-age
F_vonBertalanffy <-
  function(age, k) {
    Lasym*(1-exp(-k*(age-t0)))
  }

# Weight-at-age (standard fish allometry: W = a*L^b)
F_LengthWeight <- 
  function(L) {
    a*L^b
  }

## 2.2) Parameterization of the main model

ages          <- 1:50                                                 # Number of age cohorts 
years         <- 1980:2100                                            # number of years
N             <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for population N
catches_N     <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for catches 
catches_W     <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for the weight of catches
discards_W    <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for the weight of discards
landings_W    <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for the weight of landings
landings_V    <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for landings value
landings_tot  <- matrix(NA,nrow=length(ages), ncol=length(years))     # creation empty matrix for total landings value
rents         <- numeric(length(years))                               # creation empty vector for rents
effort        <- numeric(length(years))                               # creation empty vector for effort
Q             <- numeric(length(years))                               # creation empty vector for licenses
discards_j    <- numeric(length(years))                               # creation empty vector for the amount of juvenile discards
juvenile_discards_frac <- numeric(length(years))                      # creation empty vector for the fraction of juvenile discards (as function of total biomass)

row.names(N)            <- paste0("age", ages)
colnames(N)             <- paste0(" ", years)
row.names(catches_N)    <- paste0("age", ages)
colnames(catches_N)     <- paste0(" ", years)
row.names(catches_W)    <- paste0("age", ages)
colnames(catches_W)     <- paste0(" ", years)
row.names(discards_W)   <- paste0("age", ages)
colnames(discards_W)    <- paste0(" ", years)
row.names(landings_W)   <- paste0("age", ages)
colnames(landings_W)    <- paste0(" ", years)
row.names(landings_V)   <- paste0("age", ages)
colnames(landings_V)    <- paste0(" ", years)
row.names(landings_tot) <- paste0("age", ages)
colnames(landings_tot)  <- paste0(" ", years)


R                 <- 4477139      # recruitment per year, independent of biomass
M                 <- 0.1          # natural mortality
q                 <- 0.00082      # catchability
costs_ue          <- 693750       # costs per unit of effort per year (€)
gamma             <- 0.000001     # sensitivity of effort to rents. If zero, no change in effort     
Q_Low             <- 80           # amount of licenses when juvenile discards fraction > discard threshold
Q_High            <- 220          # amount of licenses when juvenile discards fraction < discard threshold
discard_threshold <- 0.003        # threshold for the percentage 

effort[1]         <- 150          # days at sea for year 1
N[,1]             <- 4477139 * exp(-( (q*effort[1])+M) * (ages - 1))   # declining age structure

standardDeviation <- 0.2  # st dev of disturbance term: low to reduce environmental noise


F_BiomassGrowth <-
  function(S) {
    r <- rep(growthRate,length(S))
    # In those instances where S < S_Lim, the grow rate is proportional to escapement
    r[S<S_Lim] <- (growthRate/S_Lim)*S[S<S_Lim]
    # Next year's biomass is equal to:
    Bnext = S + r*S*(1-S/maximumStock)*
      rlnorm(1, meanlog = 0, sdlog = standardDeviation)
    #Biological safeguard, so Bnext > 0 
    Bnext <- pmax(Bnext, 0.01)
    return(Bnext)
  }

## 2.3) Plot von Bertalanffy Growth curve

plot(ages, F_vonBertalanffy(ages, k), type="b",
     xlab="Age", ylab="Length (cm)",
     main="Length-at-age (discrete ages)")


## 2.4) Make a data set from the von Bertalanffy Growth curve
lengths <- F_vonBertalanffy(ages,k)

growth_df <- data.frame(
  age = ages, length_cm = lengths
)

growth_df


## 2.5) Baranov Catch Equation
baranov_fraction <- function(F, M, N){
  (1 - exp(-(F+M) * (F/(F+M)))) * N
}


# Inital biomass matrix
biomass <- matrix(NA,nrow=length(ages), ncol=length(years))
wts <- F_LengthWeight(F_vonBertalanffy(ages, k)) # Initial weights
biomass[,1] <- N[,1]*wts

## 2.6) Inital biomass matrix
biomass <- matrix(NA,nrow=length(ages), ncol=length(years))
wts <- F_LengthWeight(F_vonBertalanffy(ages, k)) # Initial weights
biomass[,1] <- N[,1]*wts


## 2.7) Discards fraction per age 
#disc_frac <- c(1,0.8,0.3,0.05,0,0,0,0,0,0)

disc_frac <- c(
  rep(1,6),
  seq(0.9,0.3,length.out=6),
  seq(0.25,0.05,length.out=8),
  seq(0.04,0.01,length.out=10),
  seq(0.008,0.001,length.out=20)
)

## 2.8) Price per age
#price_age <- c(0.8, 1.3, 2.3, 3.5, 4.6, 5.8, 6.5, 6.9, 6.2, 5.0)

price_raw <- lengths^1.2
price_age <- price_raw * (2.57 / mean(price_raw))


## 3) Main model

## 3.1) Yearly biomass loop

for(ii in 1:(length(years)-1)){
  
  # limit effort by licences
  effort[ii] <- min(effort[ii], Q[ii])
  
  # calculate F as function of q and effort
  Fi <- q * effort[ii]
  
  # calculate N from previous, takes two steps: add recruitment (age 1), have mortality for other ages
  N[1, ii+1] <-  R *  rlnorm(1, meanlog = 0, sdlog = standardDeviation)
  N[2: length(ages), ii+1] <- N[1:(length(ages)-1), ii] * exp(- (Fi +  M))
  
  # now calculate biomass
  biomass[,ii] <- N[,ii] * wts
  
  # calculate catches in numbers
  catches_N[,ii] <- baranov_fraction(Fi,M,N[,ii]) 
  
  # calculate catches weight
  catches_W[,ii] <- catches_N[,ii] * wts
  
  # calculate discards weight
  discards_W[,ii] <- catches_W[,ii] * disc_frac
  
  # calculate landings weight
  landings_W[,ii] <- catches_W[,ii] - discards_W[,ii]
  
  # calculate landings value
  landings_V[,ii] <- landings_W[,ii] * price_age
  
  # calculate rents
  rents[ii] <- sum(landings_V[,ii]) - effort[ii] * costs_ue
  
  # update effort given rents
  effort[ii+1] <- max(c(effort[ii] + (gamma * rents[ii])) , 0.001)
  
  # calculate the number of discards per year 
  discards_j[ii] <- sum(discards_W[1:2, ii])
  
  # calculate unstructured biomass per year
  S_year  <- sum(biomass[,ii])
  
  # calculate juvenile discards as a fraction of total biomass per year
  juvenile_discards_frac[ii] <- discards_j[ii]/S_year
  
  # calculate the number of licences under juvenile discard rules
  if(juvenile_discards_frac[ii] > discard_threshold) {
    Q[ii+1] <- Q_Low
  } else {
    Q[ii+1] <- Q_High
  }
  
}


## 3.1) Calculate unstructured biomass per year
S  <- colSums(biomass, na.rm = TRUE)
plot(S, type="b")
plot(N[1,], type="b")
hist(N[1,],50)


## 3.2) Calculate total landings per year
landings_tot <- apply(landings_V, 2, "sum") #total value of the landings in each year
plot(landings_tot, type="b")


## 4) Plot results

## 4.1) Plot the socio-ecological dynamics

par(mfrow = c(3,2), mar=c(4,4,2,1))

#plot(S[-length(years)], type="l", lwd=2, col="darkgreen",
#     ylab="Total biomass", xlab="Year",
#     main="Stock biomass")

plot(N[1,], type="l", lwd=2, col="steelblue",
     ylab="Recruitment (age 1)", xlab="Year",
     main="Recruitment")

plot(effort[-length(years)], type="l", lwd=2, col="firebrick",
     ylab="Effort", xlab="Year",
     main="Fishing effort")

#plot(Q[-length(years)], type="l", lwd=2, col="purple",
#     ylab="Licences", xlab="Year",
#     main="Licence dynamics")

plot(rents[-length(years)], type="l", lwd=2, col="black",
     ylab="Rents (€)", xlab="Year",
     main="Economic rents")

plot(juvenile_discards_frac[-length(years)], type="l", lwd=2, col="orange",
     ylab="Juvenile discard fraction", xlab="Year",
     main="Discard pressure")


## 4.2) Plot juvenile discards against the thresholds: shows when juvenile discards 
##      trigger license reductions

par(mfrow = c(1,1))

plot(juvenile_discards_frac, type="l", lwd=2,
     ylab="Juvenile discard fraction", xlab="Year",
     main="Juvenile discard rule")
abline(h = discard_threshold, col="red", lwd=2, lty=2)
legend("topright",  
       legend=c("Observed","Threshold"),
       col=c("black","red"),
       lty=c(1,2), lwd=2)


##4.3) Plot biomass against licences to get information about yearly management feedback  

scale_factor <- max(S) / max(Q)

par(mar = c(5,4,4,10))

plot(S, type="l", lwd=2, col="darkgreen",
     ylab="Biomass", xlab="Year",
     main="Management feedback")

par(new=TRUE)

plot(Q * scale_factor, type="l", col="purple", lwd=2,
     axes=FALSE, xlab="", ylab="")

axis(side=4,
     at = pretty(Q) * scale_factor,
     labels = pretty(Q))

mtext("Licences", side=4, line=3)

par(xpd=NA)
legend("topright", inset=c(-0.55,0),
       legend=c("Biomass","Licences"),
       col=c("darkgreen","purple"), lwd=2, bty="n")
par(xpd=FALSE)


## 4.4) Plot recruitment to check if the variability is reasonable 

plot(years, N[1,], type="l", lwd=2, col="steelblue",
     ylab="Recruitment (age 1)", xlab="Year",
     main="Recruitment variability")

hist(N[1,], breaks=40, col="lightblue",
     main="Recruitment distribution",
     xlab="Age-1 abundance")


## 5) Performance indicators

## 5.1) Calculate percentage of years that juvenile discards exceed the threshold
perc_exceedThreshold <- round(mean(juvenile_discards_frac > discard_threshold, na.rm = TRUE),3) * 100
print(paste("Percentage juvenile discards above threshold: ", perc_exceedThreshold, sep=""))

## 5.2) Calculate the number of years that juvenile discards exceed the threshold
numberofyears_exceed <- round(sum(juvenile_discards_frac > discard_threshold, na.rm = TRUE),3)
print(paste("Number of years that juvenile discards exceed threshold: ", numberofyears_exceed, sep=""))

## 5.3) Calculate in which years juvenile discards exceed the threshold
whichyears_exceed <- years[juvenile_discards_frac > discard_threshold]
print(paste("The years in which juvenile discards exceed threshold: ", whichyears_exceed, sep=""))
