#trying to run my first psc iteration file
#we'll see how this goes
# Module Imports
import os
import stochpy 
import numpy as numpy
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(80209)
start_time = 0.0
end_time = 8760
n_runs = 15
# Model output storage arrays
acquisitions = numpy.empty([n_runs,2])

# Metapopulation Model Run
def MetaPop(iteration):
	model = stochpy.SSA()
	model.Model(model_file="MattsOriginalModel.psc", dir='C:/Users/ssjoh/Documents/Lofgren/Iterations')
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[16][0][-1]
	acquisitions[iteration, 0] = Incident
   
# 1dr 2.5nurse ratio Model Run    
def FirstIt(iteration):
	model = stochpy.SSA()
	model.Model(model_file="2dr3nurseratiov2.psc", dir='C:/Users/ssjoh/Documents/Lofgren/Iterations')
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[16][0][-1]
	acquisitions[iteration, 1] = Incident

     
for i in range(0,n_runs):
    print("*** Iteration %i of %i ***" % (i+1,n_runs))
    MetaPop(i)
    FirstIt(i)

    
numpy.savetxt('blah.csv',acquisitions,delimiter=',',header="MetaPop, 1stIt",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")