
###########################################################
# Dynamic Transmission Model of MRSA in a Meta-population #
# Leaky-Pod Model					  #
# Queue-based Steady State Populations             	  #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)     #
###########################################################

# Descriptive Information for PML File
#Modelname: MRSA Meta-population Expanded
#Description: PML Implementation of MRSA transmission model 

## Most crazy iteration##
##Doing a 1:1 nurse ratio
##1 doctor, 15 nurses, with 15 cohorts##

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R3: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R4: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R5: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R6: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R7: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R8: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R9: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R10: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R11: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R12: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R13: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R14: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R15: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)
		
R16:
	N_c1 > N_u1
	N_c1 * iota_N
	
R17:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R18:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R19:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R20:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R21:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R22:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R23:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R24:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R25:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R26:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R27:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R28:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R29:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R30:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R31:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)
	
	

# Reactions Governing Movement of Nurses (N) Cohort 2 #

R32:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R33: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R34: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R35: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R36: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R37: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R38: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R39: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R40: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R41: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R42: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R43: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R44: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R45: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R46: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)
	
R47:
	N_c2 > N_u2
	N_c2 * iota_N

R48:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R49:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R50:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R51:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R52:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R53:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R54:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R55:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R56:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R57:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R58:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R59:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R60:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R61:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R62:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R63:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R64: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R65: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R66: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R67: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R68: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)

R69: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R70: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R71: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R72: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R73: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R74: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R75: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R76: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R77: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R78:
	N_c3 > N_u3
	N_c3 * iota_N

R79:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma
	
R80:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R81:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R82:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R83:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R84:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R85:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R86:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R87:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R88:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R89:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R90:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R91:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R92:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R93:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)	


# Reactions Governing Movement of Nurses (N) Cohort 4 #

R94:
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4)) * gamma

R95: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R96: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R97: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R98: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R99: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)

R100: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R101: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R102: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R103: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R104: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R105: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R106: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R107: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R108: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R109:
	N_c4 > N_u4
	N_c4 * iota_N

R110:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c4 / (P_c4 + P_u4)) * gamma

R112:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R113:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R114:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R115:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R116:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R117:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R118:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R119:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R120:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R121:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R122:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R123:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R124:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R125:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)	



# Reactions Governing Movement of Nurses (N) Cohort 5 #

R126:
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5)) * gamma

R127: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R128: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R129: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R130: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R131: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)

R132: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R133: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R134: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R135: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R136: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R137: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R138: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R139: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R140: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


R141:
	N_c5 > N_u5
	N_c5 * iota_N

R142:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c5 / (P_c5 + P_u5)) * gamma

R143:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R144:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R145:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R146:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R147:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)

R148:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R149:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R150:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R151
	N_c5 > N_u5
	N_c5 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R152:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R153:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R154:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R155:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R156:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 6 #

R157:
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c6 / (P_c6 + P_u6)) * gamma

R158: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R159: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R160: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R161: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R162: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R163: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R164: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R165: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R166: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R167: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R168: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R169: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R170: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R171: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R172:
	N_c6 > N_u6
	N_c6 * iota_N

R173:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c6 / (P_c6 + P_u6)) * gamma

R174:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R175:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R176:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R177:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R178:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R179:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R180:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R181:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R182:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R183:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R184:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R185:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R186:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R187:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 7 #

R188:
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c7 / (P_c7 + P_u7)) * gamma

R189: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R190: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R191: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R192: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R193: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R194: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R195: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R196: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R197: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R198: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R199: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R200: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R201: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R202: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R203:
	N_c7 > N_u7
	N_c7 * iota_N

R204:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c7 / (P_c7 + P_u7)) * gamma

R205:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R206:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R207:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R208:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R209:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R210:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R211:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R212:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R213:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R214:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R215:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R216:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R217:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R218:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)



# Reactions Governing Movement of Nurses (N) Cohort 8 #

R219:
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c8 / (P_c8 + P_u8)) * gamma

R220: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R221: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R222: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R223: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R224: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R225: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R226: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R227: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R228: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R229: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R230: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R231: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R232: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R233: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R234:
	N_c8 > N_u8
	N_c7 * iota_N

R235:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c8 / (P_c8 + P_u8)) * gamma

R236:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R237:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R238:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R239:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R240:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R241:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R242:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R243:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R244:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R245:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R246:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R247:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R248:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R249:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 9 #

R250:
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c9 / (P_c9 + P_u9)) * gamma

R251: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R252: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R253: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R254: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R255: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R256: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R257: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R258: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R259: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R260: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R261: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R262: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R263: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R264: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R265:
	N_c9 > N_u9
	N_c9 * iota_N

R266:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c9 / (P_c9 + P_u9)) * gamma

R267:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R268:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R269:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R270:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R271:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R272:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R273:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R274:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R275:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
	
R276:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R277:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R278:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R279:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R280:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 10 #

R281:
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c10 / (P_c10 + P_u10)) * gamma

R282: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R283: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R284: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R285: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R286: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R287: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R288: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R289: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R290: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R291: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
		
R292: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R293: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R294: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R295: 
	N_u10 > N_c10
	rho_N * sigma * N_u10 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R296:
	N_c10 > N_u10
	N_c10 * iota_N

R297:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c10 / (P_c10 + P_u10)) * gamma

R298:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R299:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R300:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R301:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R302:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R303:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R304:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R305:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R306:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R307:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)

R308:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R309:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R310:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R311:
	N_c10 > N_u10
	N_c10 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 11 #

R312:
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c11 / (P_c11 + P_u11)) * gamma

R313: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R314: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R315: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R316: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R317: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R318: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R319: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R320: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R321: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R322: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
		
R323: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R324: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R325: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R326: 
	N_u11 > N_c11
	rho_N * sigma * N_u11 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R327:
	N_c11 > N_u11
	N_c11 * iota_N

R328:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c11 / (P_c11 + P_u11)) * gamma

R329:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R330:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R331:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R332:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R333:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R334:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R335:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R336:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R337:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R338:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)

R339:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R340:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R341:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R342:
	N_c11 > N_u11
	N_c11 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 12 #

R343:
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c12 / (P_c12 + P_u12)) * gamma

R344: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R345: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R346: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R347: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R348: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R349: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R350: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R351: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R352: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R353: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
		
R354: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R355: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R356: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R357: 
	N_u12 > N_c12
	rho_N * sigma * N_u12 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R358:
	N_c12 > N_u12
	N_c12 * iota_N

R359:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c12 / (P_c12 + P_u12)) * gamma

R360:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R361:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R362:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R363:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R364:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R365:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R366:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R367:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R368:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R369:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)

R370:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R371:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R372:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R373:
	N_c12 > N_u12
	N_c12 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 13 #

R374:
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c13 / (P_c13 + P_u13)) * gamma

R375: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R376: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R377: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R378: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R379: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R380: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R381: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R382: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R383: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R384: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
		
R385: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R386: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R387: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R388: 
	N_u13 > N_c13
	rho_N * sigma * N_u13 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R389:
	N_c13 > N_u13
	N_c13 * iota_N

R390:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c13 / (P_c13 + P_u13)) * gamma

R391:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R392:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R393:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R394:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R395:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R396:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R397:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R398:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R399:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R400:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)

R401:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R402:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R403:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)
	
R404:
	N_c13 > N_u13
	N_c13 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 14 #

R405:
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c14 / (P_c14 + P_u14)) * gamma

R406: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R407: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R408: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R409: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R410: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R411: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R412: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R413: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R414: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R415: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
		
R416: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R417: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R418: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R419: 
	N_u14 > N_c14
	rho_N * sigma * N_u14 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)

R420:
	N_c14 > N_u14
	N_c14 * iota_N

R421:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c14 / (P_c14 + P_u14)) * gamma

R422:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R423:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R424:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R425:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R426:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R427:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R428:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R429:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R430:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R431:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)

R432:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R433:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R434:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R435:
	N_c14 > N_u14
	N_c14 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 14)


# Reactions Governing Movement of Nurses (N) Cohort 15 #

R436:
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c15 / (P_c15 + P_u15)) * gamma

R437: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)
	
R438: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R439: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R440: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R441: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)

R442: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R443: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R444: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R445: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R446: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)
		
R447: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R448: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R449: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R450: 
	N_u15 > N_c15
	rho_N * sigma * N_u15 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)

R451:
	N_c15 > N_u15
	N_c15 * iota_N

R452:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c15 / (P_c15 + P_u15)) * gamma

R453:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 14)
	
R454:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 14)
	
R455:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 14)
	
R456:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 14)
	
R457:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 14)

R458:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 14)
	
R459:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 14)
	
R460:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 14)
	
R461:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 14)
	
R462:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 14)

R463:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 14)
	
R464:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 14)
	
R465:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 14)
	
R466:
	N_c15 > N_u15
	N_c15 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 14)

########################

# Reactions Governing Movement of the Doctor 

R467:
	D_u > D_c
	rho_D * sigma * D_u * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15/ (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15))

R468:
	D_c > D_u
	D_c * iota_D

R469:
	D_c > D_u
	D_c * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15/ (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15))

	
########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R470:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R471:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R472:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R473:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R474:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)	

R475:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R476:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R477:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R478:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R479:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R480:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R481:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R482:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R483:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R484:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	


R485:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c / (D_c + D_u))
			
R486:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R487:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R488:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R489:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)

R490:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R491:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R492:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)	

R493:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)

R494:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R495:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R496:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R497:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R498:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R499:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R500:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R501:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R502:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	
	
R503:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c / (D_c + D_u))	

R504:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R505:	
	P_u2 > P_c2
	theta * P_u2 * nu
	

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R506:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R507:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R508:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R509:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)

R510:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)	

R511:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)

R512:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R513:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R514:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R515:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R516:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R517:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R518:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R519:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R520:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R521:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c / (D_c + D_u))	

R522:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R523:	
	P_u3 > P_c3
	theta * P_u3 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #

R524:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4)) * gamma

R525:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R526:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R527:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R528:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)	

R529:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)

R530:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R531:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R532:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R533:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R534:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R535:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R536:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R537:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R538:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R539:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c / (D_c + D_u))
	
R540:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R541:	
	P_u4 > P_c4
	theta * P_u4 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #

R542:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5)) * gamma

R543:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R544:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R545:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R546:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R547:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)

R548:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R549:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R550:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R551:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R552:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R553:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R554:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R555:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R556:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R557:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c / (D_c + D_u))	
	
R558:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R559:	
	P_u5 > P_c5
	theta * P_u5 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #

R560:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c6 / (N_c6 + N_u6)) * gamma

R561:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R562:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R563:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R564:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)	

R565:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)

R566:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R567:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R568:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R569:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R570:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R571:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R572:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R573:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R574:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R575:
	P_u6 > P_c6 + Acquisition
	rho_D * psi * P_u6 * (D_c / (D_c + D_u))	

R576:
	P_u6 > P_u6
	theta * P_u6 * (1-nu)

R577:	
	P_u6 > P_c6
	theta * P_u6 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 7 #

R578:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c7 / (N_c7 + N_u7)) * gamma

R579:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R580:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R581:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R582:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R583:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R584:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R585:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R586:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R587:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R589:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R590:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R591:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R592:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R593:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R594:
	P_u7 > P_c7 + Acquisition
	rho_D * psi * P_u7 * (D_c / (D_c + D_u))	

R595:
	P_u7 > P_u7
	theta * P_u7 * (1-nu)

R596:	
	P_u7 > P_c7
	theta * P_u7 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 8 #

R597:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c8 / (N_c8 + N_u8)) * gamma

R598:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R599:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R600:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R601:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R602:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R603:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R604:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R605:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R606:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R607:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R608:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R609:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R610:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R611:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R612:
	P_u8 > P_c8 + Acquisition
	rho_D * psi * P_u8 * (D_c / (D_c + D_u))	

R613:
	P_u8 > P_c8
	theta * P_u8 * (1-nu)

R614:	
	P_u8 > P_c8
	theta * P_u8 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 9 #

R615:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c9 / (N_c9 + N_u9)) * gamma

R616:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R617:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R618:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R619:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R620:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R621:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R622:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R623:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R624:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)	

R625:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R626:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R627:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R628:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R629:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R630:
	P_u9 > P_c9 + Acquisition
	rho_D * psi * P_u9 * (D_c / (D_c + D_u))	

R631:
	P_u9 > P_c9
	theta * P_u9 * (1-nu)

R632:	
	P_u9 > P_c9
	theta *P_u9 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 10 #

R633:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c10 / (N_c10 + N_u10)) * gamma

R634:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R635:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R636:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R637:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R638:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R639:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R640:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R641:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R642:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R643:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)

R644:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R645:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R646:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R647:
	P_u10 > P_c10 + Acquisition
	rho_N * psi * P_u10 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R648:
	P_u10 > P_c10 + Acquisition
	rho_D * psi * P_u10 * (D_c / (D_c + D_u))	

R649:
	P_u10 > P_c10
	theta * P_c10 * (1-nu)

R650:	
	P_u10 > P_c10
	theta * P_c10 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 11 #

R651:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c11 / (N_c11 + N_u11)) * gamma

R652:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R653:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R654:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R655:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R656:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R657:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R658:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R659:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R660:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R661:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)

R662:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R663:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R664:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R665:
	P_u11 > P_c11 + Acquisition
	rho_N * psi * P_u11 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R666:
	P_u11 > P_c11 + Acquisition
	rho_D * psi * P_u11 * (D_c / (D_c + D_u))	

R667:
	P_u11 > P_c11
	theta * P_c11 * (1-nu)

R668:	
	P_u11 > P_c11
	theta * P_c11 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 12 #

R669:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c12 / (N_c12 + N_u12)) * gamma

R670:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R671:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R672:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R673:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R674:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R675:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R676:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R677:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R678:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R679:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)

R680:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)	

R681:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R682:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R683:
	P_u12 > P_c12 + Acquisition
	rho_N * psi * P_u12 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R684:
	P_u12 > P_c12 + Acquisition
	rho_D * psi * P_u12 * (D_c / (D_c + D_u))	

R685:
	P_u12 > P_c12
	theta * P_c12 * (1-nu)

R686:	
	P_u12 > P_c12
	theta * P_c12 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 13 #

R687:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c13 / (N_c13 + N_u13)) * gamma

R688:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R689:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R690:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R691:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R692:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R693:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R694:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R695:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R696:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R697:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)

R698:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)	

R699:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R700:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c14 / (N_c14 + N_u14)) * ((1 - gamma) / 14)	

R701:
	P_u13 > P_c13 + Acquisition
	rho_N * psi * P_u13 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R702:
	P_u13 > P_c13 + Acquisition
	rho_D * psi * P_u13 * (D_c / (D_c + D_u))	

R703:
	P_u13 > P_c13
	theta * P_c13 * (1-nu)

R704:	
	P_u13 > P_c13
	theta * P_c13 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 14 #

R705:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c14 / (N_c14 + N_u14)) * gamma

R706:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R707:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R708:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R709:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R710:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R711:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R712:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R713:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R714:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R715:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)

R716:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)	

R717:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R718:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R719:
	P_u14 > P_c14 + Acquisition
	rho_N * psi * P_u14 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R720:
	P_u14 > P_c14 + Acquisition
	rho_D * psi * P_u14 * (D_c / (D_c + D_u))	

R721:
	P_u14 > P_c14
	theta * P_c14 * (1-nu)

R722:	
	P_u14 > P_c14
	theta * P_c14 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 15 #

R723:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c15 / (N_c15 + N_u15)) * gamma

R724:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 14)	

R725:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 14)	

R726:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 14)	

R727:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 14)	

R728:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 14)

R729:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 14)	

R730:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 14)	

R731:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 14)	

R732:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 14)	

R733:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 14)

R734:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 14)	

R735:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 14)	

R736:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 14)	

R737:
	P_u15 > P_c15 + Acquisition
	rho_N * psi * P_u15 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 14)	

R738:
	P_u15 > P_c15 + Acquisition
	rho_D * psi * P_u15 * (D_c / (D_c + D_u))	

R739:
	P_u15 > P_c15
	theta * P_c15 * (1-nu)

R740:	
	P_u15 > P_c15
	theta * P_c15 * nu
########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R741:
    P_c1 > P_u1
    mu * P_c1
    
R742:
	P_c1 > P_c1
	theta * P_c1 * nu

R743:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R744:
    P_c2 > P_u2
    mu * P_c2
    
R745:
	P_c2 > P_c2
	theta * P_c2 * nu

R746:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R747:
    P_c3 > P_u3
    mu * P_c3
    
R748:
	P_c3 > P_c3
	theta * P_c3 * nu

R749:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #

R750:
    P_c4 > P_u4
    mu * P_c4
    
R751:
	P_c4 > P_c4
	theta * P_c4 * nu

R752:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #

R753:
    P_c5 > P_u5
    mu * P_c5
    
R754:
	P_c5 > P_c5
	theta * P_c5 * nu

R755:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #

R756:
    P_c6 > P_u6
    mu * P_c6
    
R757:
	P_c6 > P_c6
	theta * P_c6 * nu

R758:
	P_c6 > P_u6
	theta * P_c6 * (1-nu)
# Reactions Involving Contaminated Patients (P_c) Cohort 7 #

R759:
    P_c7 > P_u7
    mu * P_c7
    
R760:
	P_c7 > P_c7
	theta * P_c7 * nu

R761:
	P_c7 > P_u7
	theta * P_c7 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 8 #

R762:
    P_c8 > P_u8
    mu * P_c8
    
R763:
	P_c8 > P_c8
	theta * P_c8 * nu

R764:
	P_c8 > P_u8
	theta * P_c8 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 9 #

R765:
    P_c9 > P_u9
    mu * P_c9
    
R766:
	P_c9 > P_c9
	theta * P_c9 * nu

R767:
	P_c9 > P_u9
	theta * P_c9 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 10 #

R768:
    P_c10 > P_u10
    mu * P_c10
    
R769:
	P_c10 > P_c10
	theta * P_c10 * nu

R770:
	P_c10 > P_u10
	theta * P_c10 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 11 #

R771:
    P_c11 > P_u11
    mu * P_c11
    
R772:
	P_c11 > P_c11
	theta * P_c11 * nu

R773:
	P_c11 > P_u11
	theta * P_c11 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 12 #

R774:
    P_c12 > P_u12
    mu * P_c12
    
R775:
	P_c12 > P_c12
	theta * P_c12 * nu

R776:
	P_c12 > P_u12
	theta * P_c12 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 13 #

R777:
    P_c13 > P_u13
    mu * P_c13
    
R778:
	P_c13 > P_c13
	theta * P_c13 * nu

R779:
	P_c13 > P_u13
	theta * P_c13 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 14 #

R780:
    P_c14 > P_u14
    mu * P_c14
    
R781:
	P_c14 > P_c14
	theta * P_c14 * nu

R782:
	P_c14 > P_u14
	theta * P_c14 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 15 #

R783:
    P_c15 > P_u15
    mu * P_c15
    
R784:
	P_c15 > P_c15
	theta * P_c15 * nu

R785:
	P_c15 > P_u15
	theta * P_c15 * (1-nu)




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
N_u7 = 1
N_u8 = 1
N_u9 = 1
N_u10 = 1
N_u11 = 1
N_u12 = 1
N_u13 = 1
N_u14 = 1
N_u15 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0
N_c7 = 0
N_c8 = 0
N_c9 = 0
N_c10 = 0
N_c11 = 0
N_c12 = 0
N_c13 = 0
N_c14= 0
N_c15 = 0

D_u = 1
D_c = 0

P_u1 = 1
P_u2 = 1
P_u3 = 1
P_u4 = 1
P_u5 = 1
P_u6 = 1
P_u7 = 1
P_u8 = 1
P_u9 = 1
P_u10 = 1
P_u11 = 1
P_u12 = 1
P_u13 = 1
P_u14 = 1
P_u15 = 1

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0
P_c7 = 0
P_c8 = 0
P_c9 = 0
P_c10 = 0
P_c11 = 0
P_c12 = 0
P_c13 = 0
P_c14= 0
P_c15 = 0

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
