# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:12:24 2019

@author: ssjoh
"""

# Module Imports
import os
import stochpy 
import numpy as numpy
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(80209) 
start_time = 0.0 
end_time = 8760 #43800 is the number of hours in 5 years
n_runs = 10  #Tells it to run the model 1000 times, essentially 1000 years

# Model output storage arrays
acquisitions = numpy.empty([n_runs,1]) #The second number is the number of models being compared

def Metapop_exp(iteration): 
    model = stochpy.SSA()
    model.Model(model_file='1dr_3nurseratio.psc.txt', dir=workingdir)
    model.Endtime(end_time)
    
    model.DoStochSim()
    model.ShowSpecies()