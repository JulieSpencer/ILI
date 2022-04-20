############################
# ILI Modeling Project 
# Julie A. Spencer
# This code will solve SEIR a system of ODEs aimed at modeling the outbreak curve of 
# five common upper respiratory viruses within a contained population.
# 
# The output of this code is 
# an epidemic curve to show timing and magnitude of the epidemic peak for each virus
############################

# starting with SEIR model from Keeling & Rohani (2008)

# assumptions: 
# (1) from initial infection, individuals progress to (a) hospital, (b) recovery
# (2) from the hospital, individuals progress to (a) death, (b) recovery
# (3) everyone who recovers gains full immunity
# (4) total infectious population = I + H
# (5) no coinfection is assumed
# (6) homogeneous mixing

# 1. INFLUENZA ##################################################################
# Set up a time vector
time = seq(0, 365, by=1) # in days

############################################### PARAMETERS BEGIN ###############################
# Set parameter values
beta = .0001 #  basic transmission rate dimless  (citation) 
c = 1 # transmission reduction in hospital (lit citation)
gamma1= 1/1 # rate of progression through exposed class (1/days)
gamma2 = 1/7 # rate of progression through infectious class (1/days)
gamma3 = 1/2 # rate of progression through hospitalized class (1/days)
p1 = 0.1 # proportion of infectious class that is hospitalized
p2 = 0.1 # proportion of hospitalized class that dies


#Vector of hard coded or inputed parameter values
params = c(beta, c, gamma1, gamma2, gamma3, p1, p2)
############################################ PARAMETERS END ####################################

############################################ INITIAL CONDITIONS BEGIN ####################################
vacc = 0
vacc.effect = 0
S = 9999 
SV = S*vacc*vacc.effect
Sinit = ceiling(S - SV)
Einit = 0
Iinit = 1 
Rinit = 0
Hinit = 0
Dinit = 0

xstart = c(S = Sinit, E = Einit, I = Iinit, R = Rinit, H = Hinit, D = Dinit)
############################################ INITIAL CONDITIONS END ####################################

############################################ SEIR ODE SYSTEM SET UP BEGIN #############################
# Function that evaluates the SEIR model
ILI_model <- function (t, x, params) {
  # Extract the state variables
  S = x['S'] # susceptible state
  E = x['E'] # exposed, not infectious state
  I = x['I'] # initially infectious state
  R = x['R'] # recovered state
  H = x['H'] # hospitalized state
  D = x['D'] # dead state
  
  N = S + E + I + R + H + D
  
  # Extract the parameter values
  beta = params[1]
  c = params[2]
  gamma1 = params[3]
  gamma2 = params[4]
  gamma3 = params[5]
  p1 = params[6]
  p2 = params[7]
  
  # Evaluate the ordinary differential equations at time t
  dSdt = -beta*S*(I+H) # transmission by all infected classes
  dEdt = beta*S*(I+H) - gamma1*E
  dIdt = gamma1*E - gamma2*I
  dRdt = gamma2*I*(1-p1) + gamma3*H*(1-p2)
  dHdt = gamma2*p1*I - gamma3*H
  dDdt = gamma3*p2*H
  
  # Combine into a single vector
  dxdt = c(dSdt, dEdt, dIdt, dRdt, dHdt, dDdt)
  
  # Return as a list (required by ODE solver)
  return(list(dxdt))
}
############################################ SEIR ODE SYSTEM SET UP END #############################

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R BEGIN ##############################
solns = as.data.frame(ode(xstart, time, ILI_model, params))

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R END ##############################

# 2. RESPIRATORY SYNCYTIAL VIRUS (RSV) ###################################################################
##################################################################################
# Set up a time vector
timeb = seq(0, 365, by=1) # in days

############################################### PARAMETERS BEGIN ###############################
# Set parameter values
betab = .0001 #  basic transmission rate dimless  (citation) 
cb = 1.2 # transmission reduction in hospital (lit citation)
gamma1b= 1/4 # rate of progression through exposed class (1/days)
gamma2b = 1/8 # rate of progression through infectious class (1/days)
gamma3b = 1/5 # rate of progression through hospitalized class (1/days)
p1b = 0.05 # proportion of infectious class that are hospitalized
p2b = 0.05 # proportion of hospitalized class that die


#Vector of hard coded or inputed parameter values
paramsb = c(betab, cb, gamma1b, gamma2b, gamma3b, p1b, p2b)
############################################ PARAMETERS END ####################################

############################################ INITIAL CONDITIONS BEGIN ####################################
vaccb = 0
vacc.effectb = 0
Sb = 9999 
SVb = Sb*vaccb*vacc.effectb
Sinitb = ceiling(Sb - SVb)
Einitb = 0
Iinitb = 1 
Rinitb = 0
Hinitb = 0
Dinitb = 0

xstartb = c(Sb = Sinitb, Eb = Einitb, Ib = Iinitb, Rb = Rinitb, Hb = Hinitb, Db = Dinitb)
############################################ INITIAL CONDITIONS END ####################################

############################################ SEIR ODE SYSTEM SET UP BEGIN #############################
# Function that evaluates the SEIR model
ILI_modelb <- function (t, x, paramsb) {
  # Extract the state variables
  Sb = x['Sb'] # susceptible state
  Eb = x['Eb'] # exposed, not infectious state
  Ib = x['Ib'] # initially infectious state
  Rb = x['Rb'] # recovered state
  Hb = x['Hb'] # hospitalized state
  Db = x['Db'] # dead state
  
  Nb = Sb + Eb + Ib + Rb + Hb + Db
  
  # Extract the parameter values
  betab = paramsb[1]
  cb = paramsb[2]
  gamma1b = paramsb[3]
  gamma2b = paramsb[4]
  gamma3b = paramsb[5]
  p1b = paramsb[6]
  p2b = paramsb[7]
  
  # Evaluate the ordinary differential equations at time t
  dSdtb = -betab*Sb*(Ib+Hb) # transmission by all infected classes
  dEdtb = betab*Sb*(Ib+Hb) - gamma1b*Eb
  dIdtb = gamma1b*Eb - gamma2b*Ib
  dRdtb = gamma2b*Ib*(1-p1b) + gamma3b*Hb*(1-p2b)
  dHdtb = gamma2b*p1b*Ib - gamma3b*Hb
  dDdtb = gamma3b*p2b*Hb
  
  # Combine into a single vector
  dxdtb = c(dSdtb, dEdtb, dIdtb, dRdtb, dHdtb, dDdtb)
  
  # Return as a list (required by ODE solver)
  return(list(dxdtb))
}
############################################ SEIR ODE SYSTEM SET UP END #############################

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R BEGIN ##############################
solnsb = as.data.frame(ode(xstartb, timeb, ILI_modelb, paramsb))

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R END ##############################

# 3. RHINOVIRUS ###################################################################
##################################################################################
# Set up a time vector
timec = seq(0, 365, by=1) # in days

############################################### PARAMETERS BEGIN ###############################
# Set parameter values
betac = .0001 #  basic transmission rate dimless  (citation) 
cc = 1 # transmission reduction in hospital (lit citation)
gamma1c= 1/2.5 # rate of progression through exposed class (1/days)
gamma2c = 1/3 # rate of progression through infectious class (1/days)
gamma3c = 1/2 # rate of progression through hospitalized class (1/days)
p1c = 0.01 # proportion of infectious class that are hospitalized
p2c = 0.01 # proportion of hospitalized class that die


#Vector of hard coded or inputed parameter values
paramsc = c(betac, cc, gamma1c, gamma2c, gamma3c, p1c, p2c)
############################################ PARAMETERS END ####################################

############################################ INITIAL CONDITIONS BEGIN ####################################
vaccc = 0
vacc.effectc = 0
Sc = 10000 
SVc = Sc*vaccc*vacc.effectc
Sinitc = ceiling(Sc - SVc)
Einitc = 0
Iinitc = 1 
Rinitc = 0
Hinitc = 0
Dinitc = 0

xstartc = c(Sc = Sinitc, Ec = Einitc, Ic = Iinitc, Rc = Rinitc, Hc = Hinitc, Dc = Dinitc)
############################################ INITIAL CONDITIONS END ####################################

############################################ SEIR ODE SYSTEM SET UP BEGIN #############################
# Function that evaluates the SEIR model
ILI_modelc <- function (t, x, paramsc) {
  # Extract the state variables
  Sc = x['Sc'] # susceptible state
  Ec = x['Ec'] # exposed, not infectious state
  Ic = x['Ic'] # initially infectious state
  Rc = x['Rc'] # recovered state
  Hc = x['Hc'] # hospitalized state
  Dc = x['Dc'] # dead state
  
  Nc = Sc + Ec + Ic + Rc + Hc + Dc
  
  # Extract the parameter values
  betac = paramsc[1]
  cc = paramsc[2]
  gamma1c = paramsc[3]
  gamma2c = paramsc[4]
  gamma3c = paramsc[5]
  p1c = paramsc[6]
  p2c = paramsc[7]
  
  # Evaluate the ordinary differential equations at time t
  dSdtc = -betac*Sc*(Ic+Hc) # transmission by all infected classes
  dEdtc = betac*Sc*(Ic+Hc) - gamma1c*Ec
  dIdtc = gamma1c*Ec - gamma2c*Ic
  dRdtc = gamma2c*Ic*(1-p1c) + gamma3c*Hc*(1-p2c)
  dHdtc = gamma2c*p1c*Ic - gamma3c*Hc
  dDdtc = gamma3c*p2c*Hc
  
  # Combine into a single vector
  dxdtc = c(dSdtc, dEdtc, dIdtc, dRdtc, dHdtc, dDdtc)
  
  # Return as a list (required by ODE solver)
  return(list(dxdtc))
}
############################################ SEIR ODE SYSTEM SET UP END #############################

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R BEGIN ##############################
solnsc = as.data.frame(ode(xstartc, timec, ILI_modelc, paramsc))

############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R END ##############################

# # 4. HUMAN CORONAVIRUS (HCV) ###################################################################
# ##################################################################################
#Set up a time vector
timed = seq(0, 365, by=1) # in days
# 
# ############################################### PARAMETERS BEGIN ###############################
# # Set parameter values
betad = .008 #  basic transmission rate dimless  (citation) 
cd = 1 # transmission reduction in hospital (lit citation)
gamma1d= 1/4 # rate of progression through exposed class (1/days)
gamma2d = 1/8 # rate of progression through infectious class (1/days)
gamma3d = 1/5 # rate of progression through hospitalized class (1/days)
p1d = 0.05 # proportion of infectious class that are hospitalized
p2d = 0.05 # proportion of hospitalized class that die
# 
# 
# #Vector of hard coded or inputed parameter values
paramsd = c(betad, cd, gamma1d, gamma2d, gamma3b, p1d, p2d)
# ############################################ PARAMETERS END ####################################
# 
# ############################################ INITIAL CONDITIONS BEGIN ####################################

vaccd = 0
vacc.effectd = 0
Sd = 10000 
SVd = Sd*vaccd*vacc.effectd
Sinitd = ceiling(Sd - SVd)
Einitd = 0
Iinitd = 1 
Rinitd = 0
Hinitd = 0
Dinitd = 0
 
xstartd = c(Sd = Sinitd, Ed = Einitd, Id = Iinitd, Rd = Rinitd, Hd = Hinitd, Dd = Dinitd)
############################################ INITIAL CONDITIONS END ####################################
 
############################################ SEIR ODE SYSTEM SET UP BEGIN #############################
# Function that evaluates the SEIR model
ILI_modeld <- function (t, x, paramsd) {
   # Extract the state variables
   Sd = x['Sd'] # susceptible state
   Ed = x['Ed'] # exposed, not infectious state
   Id = x['Id'] # initially infectious state
   Rd = x['Rd'] # recovered state
   Hd = x['Hd'] # hospitalized state
   Dd = x['Dd'] # dead state
   
   Nd = Sd + Ed + Id + Rd + Hd + Dd
   
# Extract the parameter values
   betad = paramsd[1]
   cd = paramsd[2]
   gamma1d = paramsd[3]
   gamma2d = paramsd[4]
   gamma3d = paramsd[5]
   p1d = paramsd[6]
   p2d = paramsd[7]
   
# Evaluate the ordinary differential equations at time t
   dSdtd = -betad*Sd*(Id+Hd) # transmission by all infected classes
   dEdtd = betad*Sd*(Id+Hd) - gamma1d*Ed
   dIdtd = gamma1d*Ed - gamma2d*Id
   dRdtd = gamma2d*Id*(1-p1d) + gamma3d*Hd*(1-p2d)
   dHdtd = gamma2d*p1d*Id - gamma3d*Hd
   dDdtd = gamma3d*p2d*Hd
   
# Combine into a single vector
   dxdtd = c(dSdtd, dEdtd, dIdtd, dRdtd, dHdtd, dDdtd)
   
# Return as a list (required by ODE solver)
   return(list(dxdtd))
 }
# ############################################ SEIR ODE SYSTEM SET UP END #############################
 
# ############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R BEGIN ##############################
solnsd = as.data.frame(ode(xstartd, timed, ILI_modeld, paramsd))
 
# ############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R END ##############################

# # 5. ADENOVIRUS ###################################################################
# ##################################################################################
# # Set up a time vector
timee = seq(0, 365, by=1) # in days
 
# ############################################### PARAMETERS BEGIN ###############################
 # Set parameter values
betae = .0001 #  basic transmission rate dimless  (citation) 
ce = 1 # transmission reduction in hospital (lit citation)
gamma1e= 1/5 # rate of progression through exposed class (1/days)
gamma2e = 1/12 # rate of progression through infectious class (1/days)
gamma3e = 1/5 # rate of progression through hospitalized class (1/days)
p1e = 0.1 # proportion of infectious class that are hospitalized
p2e = 0.2 # proportion of hospitalized class that die
# 
# 
#Vector of hard coded or inputed parameter values
paramse = c(betae, ce, gamma1e, gamma2e, gamma3e, p1e, p2e)
# ############################################ PARAMETERS END ####################################
 
# ############################################ INITIAL CONDITIONS BEGIN ####################################
vacce = 0
vacc.effecte = 0
Se = 10000 
SVe = Se*vacce*vacc.effecte
Sinite = ceiling(Se - SVe)
Einite = 0
Iinite = 1 
Rinite = 0
Hinite = 0
Dinite = 0
 
xstarte = c(Se = Sinite, Ee = Einite, Ie = Iinite, Re = Rinite, He = Hinite, De = Dinite)
# ############################################ INITIAL CONDITIONS END ####################################

# ############################################ SEIR ODE SYSTEM SET UP BEGIN #############################
# Function that evaluates the SEIR model
ILI_modele <- function (t, x, paramse) {
# Extract the state variables
   Se = x['Se'] # susceptible state
   Ee = x['Ee'] # exposed, not infectious state
   Ie = x['Ie'] # initially infectious state
   Re = x['Re'] # recovered state
   He = x['He'] # hospitalized state
   De = x['De'] # dead state
   
   Ne = Se + Ee + Ie + Re + He + De
   
# Extract the parameter values
   betae = paramse[1]
   ce = paramse[2]
   gamma1e = paramse[3]
   gamma2e = paramse[4]
   gamma3e = paramse[5]
   p1e = paramse[6]
   p2e = paramse[7]
   
#   Evaluate the ordinary differential equations at time t
    dSdte = -betae*Se*(Ie+He) # mass action contact by all infected classes
    dEdte = betae*Se*(Ie+He) - gamma1e*Ee
    dIdte = gamma1e*Ee - gamma2e*Ie
    dRdte = gamma2e*Ie*(1-p1e) + gamma3e*He*(1-p2e)
    dHdte = gamma2e*p1e*Ie - gamma3e*He
    dDdte = gamma3e*p2e*He
   
#   Combine into a single vector
    dxdte = c(dSdte, dEdte, dIdte, dRdte, dHdte, dDdte)
   
# Return as a list (required by ODE solver)
    return(list(dxdte))
 }
# ############################################ SEIR ODE SYSTEM SET UP END #############################
 
# ############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R BEGIN ##############################
solnse = as.data.frame(ode(xstarte, timee, ILI_modele, paramse))
 
# ############################## SOLVING THE MODEL WITH THE ODE SOLVER BUILT IN R END ##############################
 


############################ PLOTS FOR ALL VIRUSES #####################################

plot((I+H)~time, data=solns, type='l', main=paste('Infected Individuals'),
     xlab='Time in Days', ylab='Population', col="red",lwd=2.0, xlim=c(0,100), ylim=c(0,10000))
lines((Ib+Hb)~timeb, data=solnsb,type='l', col="blue",lwd=2.0)
lines((Ic+Hc)~timec, data=solnsc, type='l', col="orange", lwd=2.0)
lines((Id+Hd)~timed, data=solnsd, type='l', col="green", lwd=2.0)
lines((Ie+He)~timee, data=solnse,type='l', col="purple",lwd=2.0)
legend('topright', c('Influenza', 'RSV', 'Rhinovirus', 'HCV', 'Adenovirus'), lty=c(1,1,1,1,1), lwd=c(2.5,2.5,2.5,2.5,2.5), col=c('red','blue','orange','green','purple'))


