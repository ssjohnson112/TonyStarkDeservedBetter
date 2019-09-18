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
n_runs = 100
# Model output storage arrays
acquisitions = numpy.empty([n_runs,4])

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
	model.Model(model_file="1dr3nurseratio.psc", dir='C:/Users/ssjoh/Documents/Lofgren/Iterations')
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[14][0][-1]
	acquisitions[iteration, 1] = Incident

def SecondIt(iteration):
	model = stochpy.SSA()
	model.Model(model_file="2dr3nurseratiov2.psc", dir='C:/Users/ssjoh/Documents/Lofgren/Iterations')
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[16][0][-1]
	acquisitions[iteration, 2] = Incident

def ThirdIt(iteration):
	model = stochpy.SSA()
	model.Model(model_file="3dr3nurseratio.psc", dir='C:/Users/ssjoh/Documents/Lofgren/Iterations')
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[18][0][-1]
	acquisitions[iteration, 3] = Incident

     
for i in range(0,n_runs):
    print("*** Iteration %i of %i ***" % (i+1,n_runs))
    MetaPop(i)
    FirstIt(i)
    SecondIt(i)
    ThirdIt(i)

    
numpy.savetxt('3nursemodelsgamma.85psi.0464.csv',acquisitions,delimiter=',',header="MetaPop, 1stIt, 2ndIt, 3rdIt",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")