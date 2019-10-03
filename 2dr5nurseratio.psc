
###########################################################
# Dynamic Transmission Model of MRSA in a Meta-population #
# Leaky-Pod Model					  #
# Queue-based Steady State Populations             	  #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)     #
###########################################################

# Descriptive Information for PML File
Modelname: MRSA Meta-population Expanded
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
#new iteration#
# Trying the 5 pt to nurse ratio
#and 2 doctor right now
# 15 pts = 3 cohorts/3nurses

Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 2)
	
R3: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 2)
		
R4:
	N_c1 > N_u1
	N_c1 * iota_N
	
R5:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R6:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 2)
	
R7:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 2)
	


# Reactions Governing Movement of Nurses (N) Cohort 2 #

R8:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R9: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 2)
	
R10: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 2)
	
R11:
	N_c2 > N_u2
	N_c2 * iota_N

R12:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R13:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 2)
	
R14:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 2)


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R15:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R16: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 2)
	
R17: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 2)

R18:
	N_c3 > N_u3
	N_c3 * iota_N

R19:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma
	
R20:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 2)
	
R21:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 2)

########################

# Reactions Governing Movement of the Doctor 

R22:
	D_u1 > D_c1
	rho_D * sigma * D_u1 * (P_c1 + P_c2 + P_c3 / (P_c1 + P_c2 + P_c3 + P_u1 + P_u2 + P_u3))

R23:
	D_c1 > D_u1
	D_c1 * iota_D

R24:
	D_c1 > D_u1
	D_c1 * tau_D * (P_c1 + P_c2 + P_c3 / (P_c1 + P_c2 + P_c3 + P_u1 + P_u2 + P_u3))

R25:
	D_u2 > D_c2
	rho_D * sigma * D_u2 * (P_c1 + P_c2 + P_c3 / (P_c1 + P_c2 + P_c3 + P_u1 + P_u2 + P_u3))
R26:
	D_c2 > D_u2
	D_c2 * iota_D

R27:
	D_c2 > D_u2
	D_c2 * tau_D * (P_c1 + P_c2 + P_c3 / (P_c1 + P_c2 + P_c3 + P_u1 + P_u2 + P_u3))
	
########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R28:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R29:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 2)	

R30:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 2)	

R31:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))
			
R32:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R33:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R34:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R35:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 2)

R36:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 2)	

R37:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))	

R38:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R39:	
	P_u2 > P_c2
	theta * P_u2 * nu
	

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R40:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R41:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 2)	

R42:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 2)	

R43:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))	

R44:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R45:	
	P_u3 > P_c3
	theta * P_u3 * nu


########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R46:
    P_c1 > P_u1
    mu * P_c1
    
R47:
	P_c1 > P_c1
	theta * P_c1 * nu

R48:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R49:
    P_c2 > P_u2
    mu * P_c2
    
R50:
	P_c2 > P_c2
	theta * P_c2 * nu

R51:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R52:
    P_c3 > P_u3
    mu * P_c3
    
R53:
	P_c3 > P_c3
	theta * P_c3 * nu

R54:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0

D_u1 = 1
D_u2 = 1

D_c1 = 0
D_c2 = 0

P_u1 = 5
P_u2 = 5
P_u3 = 5

P_c1 = 0
P_c2 = 0
P_c3 = 0


Acquisition = 0

# Contact Rates and Contamination Probabilities #
rho_N = 3.973 # nurse direct care tasks per patient per hour
rho_D = 0.181 # doctor direct care tasks per patient per hour 
sigma = 0.054 # hand contamination probability
psi = 0.0464 # successful colonization of an uncolonized patient probability

#gamma = 1.0 proportion of time a nurse is in the assigned/original cohort
gamma = 0.85

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 # proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 6.404 #11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 1.748 #3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 2.728 #3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.744 #0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 # natural decolonization rate median 20 days per Star*ICU trial
