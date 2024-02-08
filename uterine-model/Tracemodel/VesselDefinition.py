import numpy as np

#########################
#########Fixed parameters
###########################
#Distance from 'inlet' to insonation site
Dist2Inson=90 #mm
#Blood viscosity (Pa.s)
mu=3.4e-3
#Blood density (g/mm3)
rho=1.05e-03
#Vascular elastance parameters
h=0.1 #no units

#########################
#########Variable parameters
###########################
#Vessel elastance
E=1.50e+06 #Pa
#Definition of geometry
#Generation |Number of vessels at this level | Vessel Radius (mm) | Vessel length (mm) | 
#To eliminate anastomoses simply give them zero length (note you also need to do this in baseline case)
vessels = np.array([(1, 1, 2.0, 100.0,'Uterine'),(2, 1, 2.0, 18.0,'Arcuate'),(3, 50, 0.2, 6.0,'Radial'),(4, 50, 0.2, 9.0,'Anastomose')],
                  dtype=[('generation', 'i4'),('number', 'i4'),('radius', 'f8'),('length', 'f8'),('vessel_type', 'U10')])
#spirals and IVS are defined by  resistance (Pa.s/mm) and compliance (/Pa) [R|C|0=off 1=on]
#To remove these from the model (eg post partum) set third parameter to be zero, otherwise set as 1 
SA_IVS = np.array([1.6,1e-8,1]);

#steady flow component (baseline) in ml/min - will be scaled with resistance
SteadyFlow=249.0


#Option to plot waveform to screen
plotv='y'


##BASELINE VALUES FOR COMPARISON - Not necessary to change unless you change the structure of the geometry(!)
#Definition of geometry
#Generation |Number of vessels at this level | Vessel Radius (mm) | Vessel length (mm) | 
vessels_bl = np.array([(1, 1, 2.0, 100.0,'Uterine'),(2, 1, 2.0, 18.0,'Arcuate'),(3, 50, 0.2, 6.0,'Radial'),(4, 50, 0.2, 9.0,'Anastomose')],
                  dtype=[('generation', 'i4'),('number', 'i4'),('radius', 'f8'),('length', 'f8'),('vessel_type', 'U10')])
#spirals and IVS are defined by  resistance and compliance [R|C]
SA_IVS_bl = np.array([1.6,1e-8,1]);

