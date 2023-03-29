# Code for running the different IBM models

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the codes to perform and plot the simulations of different individual based models.
All the code is made on Julia (*1.7.3*).

Here is the different models currently available:

1. Kefi model for vegetation dynamics in drylands (Kefi et al., TPB 2007)
2. Schneider extension of the Kefi model (Schneider et al., TE 2016)
3. Contact process with positive feedback (Eby et al., GEB 2017)
4. Guichard model for mussel disturbance (Guichard et al., Am' Nat' 2004)


```julia

include("IBM_models.jl") #loading the functions



model = "Eby" #other include keif, Eby and schneider. 
param = Get_parameters(model) #getting the parameters. It returns a dictionary

landscape = Get_initial_lattice(model,size_landscape=75) #initial landscape

param["p"]=.4
param["q"]=.4
dyn, land = Run_model(model=model, param=copy(param),landscape=copy(landscape),
    intensity_feedback = 6) #running the dynamics

#Ploting functions
Plot_dynamics(model, dyn)
Plot_landscape(model, land)

```
