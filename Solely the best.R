
#AFI34806: Modelling Marine Socio-ecological systems

#title: Group work: North Sea Flatfish
#group 3
#authors: Eden Baas, Eline van Gelderen & Manu Krishnadath 
#date: 03/02/2026
  
### 1) Load packages


### 2) Make a population N of 10 age cohorts for years 2080-2030. 
ages          <- 1:10     #10 population age classes 
years         <- 1980:2030
N             <- matrix(NA,nrow=length(ages), ncol=length(years))
row.names(N)  <- paste0("age", ages)
colnames(N)   <- paste0(" ", years)
R             <- 100
M             <- 0.15
F             <- 0.1      #later F=q*effort
N[,1]         <- 10 * exp(-0.5 * (ages - 1))   # declining age structure

wts1 <- c(0.1,0.2,0.3,0.4,0.45,0.5,0,5,0.5,0.5)

biomass1 <- sweep(N,1,wts, FUN= "*")


###### -----------
ages <- 1:10
years <- 1980:2030
N <- matrix(NA,nrow=length(ages), ncol=length(years))
R <- 100
M <- 0.1
N[1,] <- R
F <- 0.1  #later F=q*effort
for(ii in 2:length(years)){
  for(aa in 2:length(ages))
    N[aa,ii] <- N[aa-1,ii-1]*exp(-(M + F))
}
wts <- c(0.1,0.2,0.3,0.4,0.45,0.5,0,5,0.5,0.5)

biomass1 <- sweep(N,1,wts, FUN= "*")

########3-----------







## 2.1) Add density-dependent growth to the model with the von Bertalanffy Growth curve

Lasym <- 100   #asympotic length (cm) for plaice
k     <- 0.2   #growth coefficient
t0    <- 0     #theoretical age at zero length

a <- 0.009
b <- 3

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

#Density-dependent growth: adjust k by total biomass
k_eff <- 
  function(B) {
    k*exp(-0.002*B)
  }


# Biomass growth / recruitment
growthRate         <- 0.1  # Maximum growth rate
maximumStock       <- 100  # Carrying capacity
S_Lim              <- 30   # Threshold escapement
standardDeviation  <- 0.2  # St dev of disturbance term: low to reduce environmental noise

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

# Inital biomass matrix
biomass <- matrix(NA,nrow=length(ages), ncol=length(years))
wts <- F_LengthWeight(F_vonBertalanffy(ages, k)) # Initial weights
biomass[,1] <- N[,1]*wts

#biomass <- sweep(N,1,wts, FUN= "*") ??

# Yearly biomass loop
for(ii in 2:length(years)){
  
  # Total biomass from previous year
  S_prev <- sum(biomass[,ii-1], na.rm=TRUE)
  
  # Density-dependent growth
  k_year <- k_eff(S_prev)
  
  #Update weights at age for this year
  wts <- F_LengthWeight(F_vonBertalanffy(ages, k_year))
  
  # Biomass-at-age last year using current weights
  biomass_prev <- N[,ii-1] * wts
  
  # Biomass growth 
  B_next <- F_BiomassGrowth(sum(biomass_prev, na.rm = TRUE))
  
  # Recruitment from biomass 
  recruitment_scale <- 0.05
  N[1,ii] <- max(recruitment_scale * B_next / sum(wts), 0.01)
  
  # Survival of older ages (result of natural mortality and fishing)
  for(aa in 2:length(ages))
    N[aa,ii] <- N[aa-1,ii-1]*exp(-(M + F))
  
  # Update biomass matrix
  biomass[,ii] <- N[,ii]*wts
}




# Set the number of licences under TAC
F_MNG_propRule <-
  function(B){
    # In principle the number of licenses is high
    Q <- rep(Q_High,length(B))
    # When B < B_Trigger, the number of licenses is low
    Q[B<=B_Trigger] <- Q_Low
  }

Q_High             <- Inf    # Number of licenses if B > B_Trigger
B_Trigger          <- 0      # Precautionary biomass
Q_Low              <- Inf    # Number of licenses if B <= B_Trigger


# 4) Fishing behaviour
F_FLT_entryExit_C <-
  function(B,E,R,Q) {
    # Effort changes according to last year's rents
    out_effort                 <- E+entryExitParm*R
    # 1. If no effort last year -> effort = 'entryFromZero'
    out_effort[E == 0]         <- entryFromZero
    # 2. If effort is negative -> set E=0
    out_effort[out_effort<0]   <- 0
    # 3. If effort > licenses -> set effort = licenses
    out_effort[out_effort>Q] <- Q[out_effort>Q]
    return(out_effort)
  }

entryFromZero       <- 10
entryExitParm       <- 0.5





# -------------------------------
# 1) Parameters
# -------------------------------
ages        <- 1:10
years       <- 1980:2030
n_ages      <- length(ages)
n_years     <- length(years)

# Initial population
N <- matrix(NA, nrow=n_ages, ncol=n_years)
row.names(N) <- paste0("age", ages)
colnames(N) <- years

N[,1] <- 10 * exp(-0.5 * (ages-1))  # declining initial age structure

# Von Bertalanffy growth
Lasym <- 100   # cm
k_base <- 0.2
t0 <- 0

# Weight-at-age
a <- 0.01
b <- 3

# Mortality
M <- 0.1
F <- 0.1

# Density-dep growth
k_eff <- function(B){
  k_min <- 0.05
  k_adj <- k_base * exp(-0.002 * B)
  return(max(k_adj, k_min))
}

# Weight-at-age function
F_vonBertalanffy <- function(age, k){
  Lasym * (1 - exp(-k * (age - t0)))
}

F_LengthWeight <- function(L){
  a * L^b
}

# Biomass growth function
growthRate        <- 0.1
maximumStock      <- 100
standardDeviation <- 0.2

F_BiomassGrowth <- function(S){
  growth_term <- growthRate * S * (1 - S / maximumStock)
  growth_term <- max(growth_term, 0)
  Bnext <- S + growth_term * rlnorm(1, meanlog=0, sdlog=standardDeviation)
  Bnext <- max(Bnext, 0.01)
  return(Bnext)
}

# Recruitment scaling
recruitment_scale <- 0.1
recruitment_max   <- 500

# -------------------------------
# 2) Initialize biomass matrix
# -------------------------------
biomass <- matrix(NA, nrow=n_ages, ncol=n_years)
wts <- F_LengthWeight(F_vonBertalanffy(ages, k_base))
biomass[,1] <- N[,1] * wts

# -------------------------------
# 3) Yearly loop
# -------------------------------
for(ii in 2:n_years){
  
  # Total biomass previous year
  S_prev <- sum(biomass[,ii-1], na.rm=TRUE)
  S_prev <- min(S_prev, maximumStock)  # cap at carrying capacity
  
  # Density-dependent growth
  k_year <- k_eff(S_prev)
  
  # Update weight-at-age
  wts <- F_LengthWeight(F_vonBertalanffy(ages, k_year))
  
  # Biomass last year at current weights
  biomass_prev <- N[,ii-1] * wts
  
  # Biomass growth
  B_next <- F_BiomassGrowth(sum(biomass_prev))
  
  # Recruitment (age-1)
  recruit <- recruitment_scale * B_next / sum(wts)
  recruit <- max(min(recruit, recruitment_max), 0.01)
  N[1,ii] <- recruit
  
  # Survival of older ages
  for(aa in 2:n_ages){
    N[aa,ii] <- N[aa-1,ii-1] * exp(-(M + F))
  }
  
  # Update biomass matrix
  biomass[,ii] <- N[,ii] * wts
}

# -------------------------------
# 4) Quick check
# -------------------------------
head(N)
head(biomass)
