
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
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####
## 15 Bed ICU ##
## 2.5 ratio to nurses to pts##
## 3 pods of 3 patients, 3 of 2 pts##
## 3 Doctors##

# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R3: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R4: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R5: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R6: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
		
R7:
	N_c1 > N_u1
	N_c1 * iota_N
	
R8:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R9:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R10:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R11:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R12:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R13:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
	

# Reactions Governing Movement of Nurses (N) Cohort 2 #

R14:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R15: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R16: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R17: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R18: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
	
R19: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R20:
	N_c2 > N_u2
	N_c2 * iota_N

R21:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R22:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R23:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R24:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R25:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R26:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R27:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R28: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R29: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R30: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
	
R31: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R32: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)

R33:
	N_c3 > N_u3
	N_c3 * iota_N

R34:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma
	
R35:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R36:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R37:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R38:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R39:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)	


# Reactions Governing Movement of Nurses (N) Cohort 4 #

R40:
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4)) * gamma

R41: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R42: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
	
R43: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R44: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R45: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)

R46:
	N_c4 > N_u4
	N_c4 * iota_N

R47:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c4 / (P_c4 + P_u4)) * gamma

R48:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R49:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R50:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R51:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R52:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)


# Reactions Governing Movement of Nurses (N) Cohort 5 #

R53:
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5)) * gamma

R54: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)
	
R55: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R56: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R57: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R58: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)

R59:
	N_c5 > N_u5
	N_c5 * iota_N

R60:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c5 / (P_c5 + P_u5)) * gamma

R61:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R62:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R63:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R64:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R65:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 5)


# Reactions Governing Movement of Nurses (N) Cohort 6 #

R66:
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c6 / (P_c6 + P_u6)) * gamma

R67: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)
	
R68: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R69: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R70: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R71: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)

R72:
	N_c6 > N_u6
	N_c6 * iota_N

R73:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c6 / (P_c6 + P_u6)) * gamma

R74:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 5)
	
R75:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 5)
	
R76:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 5)
	
R77:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 5)
	
R78:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 5)


########################

# Reactions Governing Movement of the Doctor 

R79:
	D_u1 > D_c1
	rho_D * sigma * D_u1 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))

R80:
	D_c1 > D_u1
	D_c1 * iota_D

R81:
	D_c1 > D_u1
	D_c1 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))

R82:
	D_u2 > D_c2
	rho_D * sigma * D_u2 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))

R83:
	D_c2 > D_u2
	D_c2 * iota_D

R84:
	D_c2 > D_u2
	D_c2 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))

R85:
	D_u3 > D_c3
	rho_D * sigma * D_u3 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))

R86:
	D_c3 > D_u3
	D_c3 * iota_D

R87:
	D_c3 > D_u3
	D_c3 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6))
	
########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R88:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R89:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 5)	

R90:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 5)	

R91:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 5)	

R92:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 5)	

R93:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 5)	

R94:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))* (D_c3 / (D_c3 + D_u3))
			
R95:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R96:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R97:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R98:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 5)

R99:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 5)	

R100:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 5)	

R101:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 5)	

R102:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 5)
	
R103:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R104:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R105:	
	P_u2 > P_c2
	theta * P_u2 * nu
	

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R106:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R107:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 5)	

R108:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 5)	

R109:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 5)

R110:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 5)	

R111:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 5)

R112:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		

R113:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R114:	
	P_u3 > P_c3
	theta * P_u3 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #

R115:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4)) * gamma

R116:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 5)	

R117:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 5)	

R118:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 5)	

R119:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 5)	

R120:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 5)

R121:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	
	
R122:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R123:	
	P_u4 > P_c4
	theta * P_u4 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #

R124:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5)) * gamma

R125:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 5)	

R126:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 5)	

R127:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 5)	

R128:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 5)	

R129:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 5)

R130:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		
	
R131:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R132:	
	P_u5 > P_c5
	theta * P_u5 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #

R133:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c6 / (N_c6 + N_u6)) * gamma

R134:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 5)	

R135:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 5)	

R136:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 5)	

R137:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 5)	

R138:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 5)

R139:
	P_u6 > P_c6 + Acquisition
	rho_D * psi * P_u6 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R140:
	P_u6 > P_u6
	theta * P_u6 * (1-nu)

R141:	
	P_u6 > P_c6
	theta * P_u6 * nu


########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R142:
    P_c1 > P_u1
    mu * P_c1
    
R143:
	P_c1 > P_c1
	theta * P_c1 * nu

R144:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R145:
    P_c2 > P_u2
    mu * P_c2
    
R146:
	P_c2 > P_c2
	theta * P_c2 * nu

R147:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R148:
    P_c3 > P_u3
    mu * P_c3
    
R149:
	P_c3 > P_c3
	theta * P_c3 * nu

R150:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #

R151:
    P_c4 > P_u4
    mu * P_c4
    
R152:
	P_c4 > P_c4
	theta * P_c4 * nu

R153:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #

R154:
    P_c5 > P_u5
    mu * P_c5
    
R155:
	P_c5 > P_c5
	theta * P_c5 * nu

R156:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #

R157:
    P_c6 > P_u6
    mu * P_c6
    
R158:
	P_c6 > P_c6
	theta * P_c6 * nu

R159:
	P_c6 > P_u6
	theta * P_c6 * (1-nu)


########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1
N_u4 = 1
N_u5 = 1
N_u6 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0

D_u1 = 1
D_u2 = 1
D_u3 = 1

D_c1 = 0
D_c2 = 0
D_c3 = 0

P_u1 = 3
P_u2 = 3
P_u3 = 3
P_u4 = 2
P_u5 = 2
P_u6 = 2

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0

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
