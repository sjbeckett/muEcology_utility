import numpy as np
import math

# functions to help evaluate properties of microbial communities
# inputs typically as microns, outputs typically at metres, seconds
#
# constants -- some constants.
# ViralQuota -- molar quotas of viruses from size
# GrazerQuota -- molar quotas of grazers from size
# dynamicViscosity -- computation of dynamic viscosity of seawater
# densitySW -- computation of seawater density
# kinematicViscosity -- computation of kinematic viscosity of seawater
# Stokes_Einstein_Sutherland -- translational diffusion
# Stokes_Einstein_Debye -- rotational diffusion
# SmoluchowskiCoagulation -- encounter via diffusion of spheres of differing radii.
# SwimmingSpeed -- microbial swimming speed as characterized by Kiørboe, 2011.
# RelativeSwimmingEncounter -- encounter via swimming of individuals of differing radii.
# Encounter -- encounter rates following equations from Talmy et al. 2019.

class Constants: 
    def __init__(constants): 
        constants.Avogadro = 6.02214076*10**23 #Avogadro's number [mol^(-1)]
        constants.Boltzmann = 1.380649*10**(-23) #Boltzmann constant [J/K = m^2.kg/(s^2.K)]
        constants.MM_C = 12.011 #molar mass carbon [g/mol]
        constants.MM_Fe = 55.845 #molar mass iron [g/mol]
        constants.MM_N = 14.007 #molar mass nitrogen [g/mol]
        constants.MM_O = 15.999 #molar mass oxygen [g/mol]
        constants.MM_P = 30.974 #molar mass phosphorous [g/mol]
        constants.MM_Si = 28.085 #molar mass silica [g/mol]


def ViralQuota(rv):
    #inputs: viral capsid radius rv in micron
    #outputs: quotas for C, N, P [mmol X per indiv] # milimole X per indiv
    rv = rv*1000 #convert from micron to nm
    #Following Jover et al. 2014 https://doi.org/10.1038/nrmicro3289
    x = (10**3/Constants().Avogadro)
    #Quotas from Jover et al. 2014: https://doi.org/10.1038/nrmicro3289
    # all in micromole X per virus
    QC = x*(41*(rv-2.5)**3 + 130*(7.5*rv**2 - 18.75*rv +15.63))
    QN = x*(16*(rv-2.5)**3 + 36*(7.5*rv**2 - 18.75*rv +15.63))
    QP = x*4.2*(rv-2.5)**3
    return QC, QN, QP


def radius_2_volume(radius):
    V = 4/3*math.pi*radius**3
    return V


def GrazerQuota(rg):
    #inputs: grazer radius rg in micron
    #output: quota for C, N [mmol X per indiv] #milimole X per indiv
    #Following dinoflagellate relationship in Menden-Deuer and Lessard, 2000. https://doi.org/10.4319/lo.2000.45.3.0569 (note may differ for other groups: e.g., https://doi.org/10.1002/lno.12284)
    mgC = 10**-9 * 0.76* radius_2_volume(rg)**0.819 #miligram C per cell
    mgN = 10**-9 * 0.118* radius_2_volume(rg)**0.849 #miligram N per cell
    QC = mgC/Constants().MM_C # milimolar C per cell
    QN = mgN/Constants().MM_N # milimolar N per cell
    return QC, QN
    

def dynamicViscosity(Temperature=15, Salinity=35):
    #inputs: Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #ouputs: dynamic viscocity of seawater [kg/m-s]
    
    # Following Sharqawy et al. 2010  https://doi.org/10.5004/dwt.2010.1079
    #from MIT seawater:  http://web.mit.edu/seawater/
    # VALIDITY: 0 < T < 180 C and 0 < S < 150 g/kg; ACCURACY: 1.5%
    
    S = Salinity
    T = Temperature
    S = S/1000
    a = np.array([1.5700386464E-01,6.4992620050E+01,-9.1296496657E+01,4.2844324477E-05,1.5409136040E+00,1.9981117208E-02,-9.5203865864E-05,7.9739318223E+00,-7.5614568881E-02,4.7237011074E-04])
    mu_w = a[3] + 1/(a[0]*(T+a[1])**2+a[2])
    A  = a[4] + a[5] * T + a[6] * T**2
    B  = a[7] + a[8] * T + a[9]* T**2
    mu = mu_w*(1 + A*S + B*S**2) # [kg/m-s]
    return mu


def densitySW(Temperature=15, Salinity=35, Pressure = 0.101325):
    #inputs: Temperature as Celcius; Salinity as g/kg; Pressure as MPa; defined for 0 < T < 180 C and 0 < S < 150 g/kg and 0 < P < 12 MPa.
    #output: denity of seawater [kg/m^3]
    #Following Nayar et al. 2016. http://dx.doi.org/10.1016/j.desal.2016.02.024
    #from MIT seawater:  http://web.mit.edu/seawater/
    T = Temperature
    S = Salinity/1000
    P = Pressure
    P0 = 0.101325 # atmospheric pressure [MPa]
    a = np.array([9.999*10**2, 2.034*10**-2, -6.162*10**-3, 2.261*10**-5, -4.657*10**-8])
    b = np.array([8.02*10**2, -2.001, 1.677*10**-2, -3.06*10**-5, -1.613*10**-5])
    psw_P0 = (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4) + (b[0]*S + b[1]*S*T + b[2]*S*T**2 + b[3]*S*T**3 + b[4]*S**2*T**2)
    cc = np.array([5.0792*10**-4, -3.4168*10**-6, 5.6931*10**-8, -3.7263*10**-10, 1.4465*10**-12, -1.7058*10**-15, -1.3389*10**-6, 4.8603*10**-9, -6.8039*10**-13])
    d = np.array([-1.1077*10**-6, 5.5585*10**-9, -4.2539*10**-11, 8.3702*10**-9])
    Fp = exp((P-P0)*(cc[0]+cc[1]*T+cc[2]*T**2+cc[3]*T**3+cc[4]*T**4+cc[5]*T**5 + S*(d[0]+d[1]*T+d[2]*T**2)) + 0.5*(P**2-P0**2)*(cc[6]+cc[7]*T+cc[8]*T**3+d[3]*S))
    rho = Fp*psw_P0
    return rho
	

def kinematicViscosity(Temperature=15, Salinity=35, Pressure = 0.101325):
    #inputs: Temperature as Celcius; Salinity as g/kg; Pressure as MPa; defined for 0 < T < 180 C and 0 < S < 150 g/kg and 0 < P < 12 MPa.
    #output: Kinematic viscosity [m^2/s]
    v = dynamicViscosity(Temperature, Salinity)
    rho = densitySW(Temperature, Salinity, Pressure)
    k = v/k
    return k


def Stokes_Einstein_Sutherland(radius, Temperature=15, Salinity=35): #translational diffusion - movement across space
    #inputs: radius as micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #output: diffusion as [m^2/s]
    DynamicViscosity = dynamicViscosity(Temperature,Salinity) # [kg/m-s]
    Diff_t = (Constants().Boltzmann*Celcius2Kelvin(Temperature))/( 6*math.pi*micron_2_metre(radius)*DynamicViscosity ) # [m^2/s]
    return Diff_t


def Stokes_Einstein_Debye(radius, Temperature=15, Salinity=35): #rotational diffusion - angular rotation
    #inputs: radius as micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #output: diffusion as [m^2/s]
    DynamicViscosity = dynamicViscosity(Temperature,Salinity) # [kg/m-s]
    Diff_r = (Constants().Boltzmann*Celcius2Kelvin(Temperature))/( 8*math.pi**micron_2_metre(radius)**3*DynamicViscosity ) # [m^2/s]
    return Diff_r


def SmoluchowskiCoagulation(radius1,radius2, Temperature=15, Salinity=35):
    #inputs:  radius1, radius2 in micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #ouput: encounter rate via diffusion between spheres of radius1 and radius2 [m^3/s]
    SC = 4*math.pi*(micron_2_metre(radius1) + micron_2_metre(radius2))*(Stokes_Einstein_Sutherland(radius1,Temperature,Salinity) + Stokes_Einstein_Sutherland(radius2,Temperature,Salinity)) # [m^3/s]
    return SC
    

def SwimmingSpeed(radius):
    #inputs: radius as micron.  Valid for 1 micron to 1 cm.
    #output: swimming speed as  m per second
    
    #Following Talmy et al. 2019. Front. Mar. Sci. https://doi.org/10.3389/fmars.2019.00182; from Kiørboe, 2011. Biological Reviews. https://doi.org/10.1111/j.1469-185X.2010.00148.x
    u = math.exp(0.39 + 0.79*math.log((2*radius)/10**4))/10**2 # [m/s]
    return u


def micron_2_metre(x):
    return x/(10**6)


def RelativeSwimmingEncounter(radius1, radius2, fdetect=3, Temperature=15, Salinity=35):
    #inputs:  radius1, radius2 as micron; fdetect as scaler: detect = fdetect*r; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #output:  encounter rate [m^3/second]
    
    #Following Talmy et al. 2019 https://doi.org/10.3389/fmars.2019.00182
    predator = micron_2_metre(max(radius1,radius2))
    prey = micron_2_metre(min(radius1,radius2))
    swim1 = SwimmingSpeed(radius1)
    swim2 = SwimmingSpeed(radius2)
    Enc = ((swim1**2+swim2**2)^0.5)*math.pi*(predator*fdetect + prey)**2  # m^3 per second
    return Enc


def EncounterRate(prey_radius, predator_radius, fdetect=3, Temperature=15, Salinity=35, Enc_Type=1):
    #inputs:  prey_radius, predator_radius as micron; fdetect as scaler: detect = fdetect*r; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
    #Enc_Type -- 1), full predator-prey setup (assuming both motile) 2), predator-prey assuming both motile and predator diffusion is negligble, 3) predator-prey assuming diffusive predator, with mobile prey.
    #output:  encounter rate [m^3/second]
    #Following Talmy et al. 2019. https://doi.org/10.3389/fmars.2019.00182
    if Enc_Type==1:  #equation 4
        Encounter = RelativeSwimmingEncounter(prey_radius,predator_radius,fdetect,Temperature,Salinity) + SmoluchowskiCoagulation(prey_radius,predator_radius,Temperature,Salinity)
    elif Enc_Type==2: #equation 5 - motile predators (grazers)
        Encounter = RelativeSwimmingEncounter(prey_radius,predator_radius,fdetect,Temperature,Salinity) + 4*math.pi*(Stokes_Einstein_Sutherland(prey_radius,Temperature,Salinity))*(prey_radius+predator_radius)
    elif Enc_Type==3: # equation 6 - diffusion predators (viruses)
        Encounter = math.pi*SwimmingSpeed(prey_radius)*(prey_radius+predator_radius)**2 + SmoluchowskiCoagulation(prey_radius,predator_radius,Temperature,Salinity)
    elif Enc_Type==4: # diffusion only
        Encounter = SmoluchowskiCoagulation(prey_radius,predator_radius,Temperature,Salinity)
    elif Enc_Type==5: # swimming only
        Encounter = RelativeSwimmingEncounter(prey_radius,predator_radius,fdetect,Temperature,Salinity)
    else:
        warnings.warn("Input for Enc_Type is not documented")
    return Encounter
