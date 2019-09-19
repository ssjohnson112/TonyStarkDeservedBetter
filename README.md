Iterations for Rotations

Matt's original code.psc- What it sounds like. Matt's orignial .psc file (gamma.psc) to play around with for all the iterations

1dr_3nurseratio .psc - psc file for first iteration. Basically taking out the 6th cohort from Matt's orignial file (going down to a 15 bed ICU instead of 18)

2dr_3nurseratio.psc - file for 2 docs with 3 nurses. 2 docs use seperate rxns to decribe movements with them and patients.

2dr_3nurseratiov2.psc - second version. 1 rxn to desribe them with pt movements, that multiplys the 2 together

3dr_3nurseratio.psc - psc file for 3docs/3nurses:pt. use rxns from v2 to describe pt and dr interaction

1dr_2.5nurseratio.psc - file with 6 pt cohorts-3x3pts, 3x2pts- and 1 doctor

Modelruns.py - running the three above models with matt's orginial code for 100 runs

1dr3nurseruns.py - running the 1dr_3nurseratio.psc with matt's code to make sure that the program to run things is working correctly. 

outcomes.py - Finding where the variable 'acquisitions' is so that each model can have it correctly listed in outcomes when models are run.
