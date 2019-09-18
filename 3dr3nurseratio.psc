###########################################################
# Dynamic Transmission Model of MRSA in a Meta-population #
# Leaky-Pod Model		                          #
# Queue-based Steady State Populations             	  #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)     #
###########################################################

# Descriptive Information for PML File
Modelname: MRSA Meta-population Expanded
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####
# 15 bed ICU
# 3:1 nurse pt ratio (5 nurses) and 2 doctors
# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R3: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R4: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R5: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)
			
R6:
	N_c1 > N_u1
	N_c1 * iota_N
	
R7:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R8:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R9:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R10:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R11:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)
	

# Reactions Governing Movement of Nurses (N) Cohort 2 #

R12:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R13: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R14: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R15: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)
	
R16: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R17:
	N_c2 > N_u2
	N_c2 * iota_N

R18:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R19:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R20:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R21:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R22:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R23:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R24: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R25: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)
	
R26: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R27: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)

R28:
	N_c3 > N_u3
	N_c3 * iota_N

R29:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma
	
R30:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R31:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R32:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R33:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)	


# Reactions Governing Movement of Nurses (N) Cohort 4 #

R34:
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4)) * gamma

R35: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)
	
R36: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R37: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R38: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)

R39:
	N_c4 > N_u4
	N_c4 * iota_N

R40:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c4 / (P_c4 + P_u4)) * gamma

R41:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R42:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R43:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R44:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 4)


# Reactions Governing Movement of Nurses (N) Cohort 5 #

R45:
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5)) * gamma
	
R46: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
R47: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R48: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R49: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)

R50:
	N_c5 > N_u5
	N_c5 * iota_N

R51:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c5 / (P_c5 + P_u5)) * gamma

R52:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 4)
	
R53:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 4)
	
R54:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 4)
	
R55:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 4)
	
########################

# Reactions Governing Movement of the Doctor 

R56:
	D_u1 > D_c1
	rho_D * sigma * D_u1 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5))

R57:
	D_c1 > D_u1
	D_c1 * iota_D

R58:
	D_c1 > D_u1
	D_c1 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5))

R59:
	D_u2 > D_c2
	rho_D * sigma * D_u2 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5))

R60:
	D_c2 > D_u2
	D_c2 * iota_D

R61:
	D_c2 > D_u2
	D_c2 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5))

R62:
	D_u3 > D_c3
	rho_D * sigma * D_u3 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5))

R63:
	D_c3 > D_u3
	D_c3 * iota_D

R64:
	D_c3 > D_u3
	D_c3 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_u1 + P_u2 + P_u3 + P_u4 +P_u5))

	
########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R65:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R66:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 4)	

R67:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 4)	

R68:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 4)	

R69:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 4)	

R70:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))

R71:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R72:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R73:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R74:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 4)

R75:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 4)	

R76:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 4)	

R77:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 4)	

R78:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))

R79:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R80:	
	P_u2 > P_c2
	theta * P_u2 * nu
	

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R80:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R81:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 4)	

R82:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 4)	

R83:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 4)

R84:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 4)	

R85:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))

R86:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R87:	
	P_u3 > P_c3
	theta * P_u3 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #

R88:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4)) * gamma

R89:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 4)	

R90:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 4)	

R91:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 4)	

R92:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 4)	

R93:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R94:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R95:	
	P_u4 > P_c4
	theta * P_u4 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #

R96:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5)) * gamma

R97:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 4)	

R98:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 4)	

R99:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 4)	

R100:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 4)	

R101:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))
	
R102:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R103:	
	P_u5 > P_c5
	theta * P_u5 * nu



########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R104:
    P_c1 > P_u1
    mu * P_c1
    
R105:
	P_c1 > P_c1
	theta * P_c1 * nu

R106:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R107:
    P_c2 > P_u2
    mu * P_c2
    
R108:
	P_c2 > P_c2
	theta * P_c2 * nu

R109:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R110:
    P_c3 > P_u3
    mu * P_c3
    
R111:
	P_c3 > P_c3
	theta * P_c3 * nu

R112:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #

R113:
    P_c4 > P_u4
    mu * P_c4
    
R114:
	P_c4 > P_c4
	theta * P_c4 * nu

R115:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #

R116:
    P_c5 > P_u5
    mu * P_c5
    
R117:
	P_c5 > P_c5
	theta * P_c5 * nu

R118:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1
N_u4 = 1
N_u5 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0

D_u1 = 1
D_u2 = 1
D_u3 = 1
D_c1 = 0
D_c2 = 0
D_c3 = 0

P_u1 = 3
P_u2 = 3
P_u3 = 3
P_u4 = 3
P_u5 = 3

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0

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
