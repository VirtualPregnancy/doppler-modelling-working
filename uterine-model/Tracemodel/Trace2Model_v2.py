import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import VesselDefinition as params
import FunctionDefinitions as funcs
from scipy import interpolate
from scipy.optimize import curve_fit, minimize,least_squares
from scipy.interpolate import  UnivariateSpline, Akima1DInterpolator, PchipInterpolator



def run_model(params): 
    ## Calculate total resistance of the system and compare to baseline (flow decreases by this factor as resistance increases assuming a constant driving pressure)
    [TotalResistance,VenousResistance]=funcs.total_resistance(params.vessels,params.SA_IVS,params)
    [BaselineResistance,BaselineVenous]=funcs.total_resistance(params.vessels_bl,params.SA_IVS_bl,params)
    SteadyFlow=params.SteadyFlow*BaselineResistance/TotalResistance

    ## Calculate characteristic admittance,wave propogation constant, and uterine artery compliance
    [CharacteristicAdmittance,WaveProp,UtCompliance]=funcs.characteristic_admittance(params.vessels,params.SA_IVS,params)
    #Calculate effective admittance and reflection coefficient for each vessel
    [EffectiveAdmittance,ReflectionCoefficients]=funcs.effective_admittance(params.vessels,params.SA_IVS,CharacteristicAdmittance,WaveProp,VenousResistance,params)

    #Calculate phase and amplitude offset for waveform at insonation site vessel    
    #Output is [[amplitude],[phase]]
    reflect_coeff=np.transpose(np.column_stack((np.absolute(ReflectionCoefficients[0][:]),np.angle(ReflectionCoefficients[0][:]))))
    char_admit=np.transpose(np.column_stack((np.absolute(CharacteristicAdmittance[0][:]),np.angle(CharacteristicAdmittance[0][:]))))

    #Real and imaginary parts of wave propogation constant  
    wave_prop_constant=np.transpose(np.column_stack((WaveProp[0][:].real,WaveProp[0][:].imag)))

    ## Incident flow profile, insonation site flow (incident, reflected, total), insonation site pressure, insonation site velocity
    [InsonationSiteVelocity,time]=funcs.timecourse(params.StartTime,params.EndTime,params.dt,reflect_coeff,char_admit,wave_prop_constant,SteadyFlow,UtCompliance,'Uterine',params)
    return time, InsonationSiteVelocity
    
def minimise_function(variable_params,params,data):
    #data is the x_coordinates and y_coordinates of the ultrasound trace

    # params.SteadyFlow = variable_params[0]
    # params.E = variable_params[1]
    params.vessels[0][2] = variable_params[0]
    params.vessels[0][3] = variable_params[1]
    params.vessels[1][2] = variable_params[2]
    params.vessels[1][3] = variable_params[3]
    params.vessels[2][2] = variable_params[4]
    params.vessels[2][3] = variable_params[5]
    params.vessels[3][2] = variable_params[6]    
    params.vessels[3][3] = variable_params[7]    
    #run model
    time, InsonationSiteVelocity = run_model(params)

    #Interpolate model predictions to data locations
    f = interpolate.interp1d(time, InsonationSiteVelocity)

    y_new = f(data[0]) #grab the model y values that correspond to the x_coordinates(ultrasound trace)
    err = data[1]-y_new #compute the difference between y_coordinates (ultrasound trace) and model y values
    
    return err
    

# Load the image containing the waveform
IMG='myImage.png'
image = cv2.imread(IMG, cv2.IMREAD_GRAYSCALE)

# Reverse the y-coordinates to match Matplotlib's coordinate system
image = np.flipud(image)

# Threshold the image to separate the waveform from the background
_, thresholded = cv2.threshold(image, 128, 255, cv2.THRESH_BINARY)

# Find contours to extract the waveform coordinates
contours, _ = cv2.findContours(thresholded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

if contours:
    # Assuming the largest contour represents the waveform
    largest_contour = max(contours, key=cv2.contourArea)
    
    x_coordinates_pixels = [point[0][0] for point in largest_contour]
    y_coordinates_pixels = [point[0][1] for point in largest_contour]

    # # Load the x and y scales computed in MATLAB
    with open('scales.txt', 'r') as file:
        lines = file.readlines()

    if len(lines) >= 2:
        time_scale = float(lines[0].strip())
        cms_scale = float(lines[1].strip())
    else:
        print('Error: The file does not contain two lines with values.')

    # Convert x-coordinates from pixels to real-world units (seconds)
    x_coordinates_seconds = [x * time_scale for x in x_coordinates_pixels]
    x_coordinates = x_coordinates_seconds

    y_coordinates_cms = [y * cms_scale for y in y_coordinates_pixels]
    y_coordinates = y_coordinates_cms


                
## Define time points at which you want to plot flows
params.dt=0.01 #time step for plotting

## Calculate total duration of the heart beat
total_duration = max(x_coordinates)-min(x_coordinates) + params.dt #adding an extra time step to this to fix an error in interpolation, may need to fix later
print('Total duration = ', total_duration, 's')
heart_rate = 60/total_duration
print('Estimated heart rate = ',  heart_rate, 'bpm')
params.HeartRate = heart_rate
#Number of flow harmonics
params.NHar=10
#First ten harmonics of incident waveform [[omega],[A_n],[Phi_n]]. See Mo et al. A transmission line modelling approach to the interpretation of uterine doppler waveforms. Ultrasound in Medicine and Biology, 1988. 14(5): p. 365-376
params.IWavHar=np.array([[1.0*params.HeartRate/60.0, 2.0*params.HeartRate/60.0, 3.0*params.HeartRate/60.0, 4.0*params.HeartRate/60.0, 5.0*params.HeartRate/60.0, 6.0*params.HeartRate/60.0, 7.0*params.HeartRate/60.0, 8.0*params.HeartRate/60.0, 9.0*params.HeartRate/60.0, 10.0*params.HeartRate/60.0],
                [64.32,  43.08,  21.48,   7.68,   2.64,  1.8,    0.96,   0.84,   0.84,   0.48],
                [-1.375319451, -2.138028334,-2.998475655,-3.548254369,-3.394665395,-3.185225885,-3.131120678,-2.448696941,-2.602285915,-2.441715624]])
params.StartTime=0.0
params.EndTime=60.0/params.HeartRate #end of a single beat

#Initial run to allow offset of data             
#time, InsonationSiteVelocity = run_model(params)

# Calculate the time offset to align the two traces
time_offset = params.StartTime - min(x_coordinates)

# Shift the x-coordinates and time to align the traces
#x_coordinates = [x - time_offset for x in x_coordinates]
x_coordinates= [t + time_offset for t in x_coordinates] 
# Sort the extracted data by "time" or x_coordinate
x_coordinates = np.array(x_coordinates)
y_coordinates = np.array(y_coordinates)
sort = np.argsort(x_coordinates)
x_coordinates = x_coordinates[sort]
y_coordinates = y_coordinates[sort]


# test variable parameters
SteadyFlow = params.SteadyFlow
E=params.E

#variable parameters
uterine_radius = params.vessels[0][2]
uterine_length = params.vessels[0][3]
arcuate_radius = params.vessels[1][2]
arcuate_length = params.vessels[1][3]
radial_radius = params.vessels[2][2]
radial_length = params.vessels[2][3]
anast_radius = params.vessels[3][2]
anast_length = params.vessels[3][3]

# variable_params = [SteadyFlow,E,uterine_radius]
variable_params = [uterine_radius,uterine_length,arcuate_radius,arcuate_length,radial_radius,radial_length,anast_radius,anast_length]

#Call the minimisation function
lb=[1,50,0.5,18,0.1,5,0.1,5]
ub=[3,150,4,180,0.25,10,3,10]
opt=least_squares(minimise_function,variable_params,bounds=(lb,ub),args=(params,[x_coordinates,y_coordinates]),xtol=1e-6)

# opt=least_squares(minimise_function,variable_params,args=(params,[x_coordinates,y_coordinates]),xtol=1e-6)

variable_params = opt.x
err = minimise_function(variable_params,params,[x_coordinates,y_coordinates])

#Run model for plotting
time, InsonationSiteVelocity = run_model(params)

#extract Image
image = mpimg.imread(IMG)

#Output to screen key properties of the flow velocity waveform    
funcs.flow_velocity_properties(InsonationSiteVelocity)    
#If selected as an option in VesselDefinition.py, plot flow velocity waveform
if(params.plotv=="y"):

    # Calculate the y-offset to align the y-coordinates
    #y_offset = min(y_coordinates) - min(InsonationSiteVelocity)

    # Shift the y-coordinates to align the traces
    #y_coordinates = [y - y_offset for y in y_coordinates]

    fig, ax1 = plt.subplots()

    ax1.plot(x_coordinates, y_coordinates, color='r',marker='+',ls='None')
    ax1.set_xlabel("Time (seconds)")
    ax1.set_ylabel("Model prediction", color="b")

    # Create a secondary y-axis for the insonation site velocity
    #ax2 = ax1.twinx()
    # Plot the insonation site velocity on the secondary y-axis
    ax1.plot(time, InsonationSiteVelocity, color='b')
    #ax1.plot(x_coordinates, y_new, color='k',marker='*',ls='None')
    
    #ax2.set_yticks([])

    plt.show()
    
#plot error points
plt.plot(x_coordinates,err,marker='+',ls='None')
plt.xlabel('Time')
plt.ylabel('Error (cm/s)')
plt.show()

# print(str(err))

#print variable parameters

print("Uterine radius " + str(params.vessels[0][2]))
print("Uterine length " + str(params.vessels[0][3]))
print("Arcuate radius " + str(params.vessels[1][2]))
print("Arcuate length " + str(params.vessels[1][3]))
print("Radial radius " + str(params.vessels[2][2]))
print("Radial length " + str(params.vessels[2][3]))
print("Anastomose radius " + str(params.vessels[3][2]))
print("Anastomose length " + str(params.vessels[3][3]))
