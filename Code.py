
from netCDF4 import Dataset    
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import cartopy.crs as ccrs

# mask water fluxes 
Pobs, LEobs, dSobs, Resobs = endo_fluxes_EXT(basin_mask)

# mask energy fluxes
DSRobs, DLRobs, USWobs, ULWobs, SHobs, NETobs = endo_energy_fluxesEXT(basin_mask)

#GRACE storage 
som_stor2020 = som_stor2020(basin_mask)

# P - E obs
PEobs = Pobs - LEobs
# Radiative fluxes
Rnobs  = DSRobs + DLRobs -USWobs -ULWobs

# FLux-inferred energy storage with yearly balance constraint
FISeYC = FISeYC(NETobs)


x = np.arange(2002,2020,1/12) # 2002 to 2020
x12 = np.arange(12)
months= ['J','F','M','A','M','J','J','A','S','O','N','D']


# Budget Closure
D = np.zeros(216) # water 
E = np.zeros(216) # energy 

opt_sol= []# array for optimised solution 
dS_opt = [] # array to append optimised storage change 
NET_opt = [] # array ro append optimised net 

# initialise with dS1 and E1
D[0] = (som_stor2020[1] - som_stor2020[0])/3.046 # unit cm 
E[0] = (FISeYC[1] - FISeYC[0])/3.046

# Loop over months to solve for NET and dS to use in D and E array 
for i in range(214):
    NET_opt = np.append(NET_opt, optimisation(i)[2])
    dS_opt = np.append(dS_opt, optimisation(i)[1])
    # start at 1 since we initialised already 
    D[i+1] = (som_stor2020[i+2] - som_stor2020[0])/3.046 - sum(dS_opt[0:i+1])
    E[i+1] = (FISeYC[i+2] - FISeYC[0])/3.046- sum(NET_opt[0:i+1])

# loop over 120 months 2002-2012
for i in range(216):    
    opt_sol = np.append(opt_sol, optimisation(i)[0])

# reshape array into the differnt fluces, note: last 2 columns are tehe lagrange multipliers    
opt = opt_sol.reshape(216,11)

# units in mm/day
Popt = opt[:,0] # precipitation
LEopt = opt[:,1] # latent heat
dSopt  = opt[:,2] # storage change
DSRopt = opt[:,3] # downwards shortwave
DLRopt = opt[:,4] # downwards longwave
USWopt = opt[:,5] # upwards shortwave
ULWopt = opt[:,6] # upwards longwave
SHopt=  opt[:,7] # sensible heat 
NETopt = opt[:,8] # NET 
PEopt =Popt-LEopt
Resopt = Popt - LEopt -dSopt
Rnopt = DSRopt +DLRopt -USWopt -ULWopt

# produce flux inferred storages 
FIS_PEobs = FISd(PEobs)
FIS_LEobs = FISd(-LEobs)
FIS_Pobs = FISd(Pobs)

FIS_Popt = FISd(Popt)
FIS_LEopt = FISd(-LEopt)
FIS_PEopt = FISd(PEopt)

FISg = FISd(dSobs)

# Desesonalise  and detrend
som_storDS = deseason(som_stor2020 )
FIS_LEobsDS  = deseason(FIS_LEobs)
FIS_PobsDS  = deseason(FIS_Pobs)
FIS_PEobsDS = deseason(FIS_PEobs)

# mean seasonal cycle - mean removed 
mmPobs = monthly_mean(Pobs - np.mean(Pobs))
mmPopt = monthly_mean(Popt - np.mean(Popt))
mmLEobs = monthly_mean(LEobs - np.mean(LEobs))
mmLEopt = monthly_mean(LEopt- np.mean(LEopt))

mmPEobs = monthly_mean(PEobs- np.mean(PEobs))
mmPEopt = monthly_mean(PEopt- np.mean(PEopt))


SIM = 2.628e6
# CONVERT UNITS
energy_Stor_opt = FISend(NETopt*28.9*SIM) 
energy_Stor_obs = FISend(NETobs*28.9*SIM)



print(scipy.stats.pearsonr(FISgDS, FIS_PoptDS))
print(scipy.stats.pearsonr(FISgDS, FIS_PobsDS))

print(scipy.stats.pearsonr(FISgDS, FIS_LEoptDS))
print(scipy.stats.pearsonr(FISgDS, FIS_LEobsDS))

print(scipy.stats.pearsonr(FISgDS, FIS_PEoptDS))
print(scipy.stats.pearsonr(FISgDS, FIS_PEobsDS))

print(scipy.stats.pearsonr(FIS_LEoptDS, FIS_PoptDS))
print(scipy.stats.pearsonr(FIS_LEobsDS, FIS_PobsDS))


Pcmyear  = np.zeros(18)
dScmyear = np.zeros(18)

for i in range(18):  
    Pcmyear[i] = sum(Pobs[i*12:i*12 + 12]*3.046)
    dScmyear[i] = sum(dSobs[i*12:i*12 + 12]*3.046)

coeffsP = np.polyfit(XY2, Pcmyear , 1)
coeffsdS = np.polyfit(XY2, dScmyear , 1)

Ptrend  = XY2*coeffsP[0] + coeffsP[1]
dStrend  = XY2*coeffsdS[0] + coeffsdS[1]



