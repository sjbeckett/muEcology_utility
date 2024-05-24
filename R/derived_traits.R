# functions to help evaluate properties of microbial communities
# inputs typically as microns, outputs typically at metres, seconds
#
# constants
# quotaJover2014
# dynamicViscosity
# Stokes_Einstein_Sutherland
# Stokes_Einstein_Debye
# SmoluchowskiCoagulation - RelativeDiffusingEncounter
# SwimmingSpeed
# RelativeSwimmingEncounter
# Encounter


constants <- function(){
	constants = c()
	constants$Avogadro = 6.02214076*10^23 #Avogadro's number [mol^(-1)]
	constants$Boltzmann = 1.380649*10^(-23) #Boltzmann constant [J/K = m^2.kg/(s^2.K)]
	constants$MM_C = 12.011 #molar mass carbon [g/mol]
	constants$MM_Fe = 55.845 #molar mass iron [g/mol]
	constants$MM_N = 14.007 #molar mass nitrogen [g/mol]
	constants$MM_O = 15.999 #molar mass oxygen [g/mol]
	constants$MM_P = 30.974 #molar mass phosphorous [g/mol]
	constants$MM_Si = 28.085 #molar mass silica [g/mol]
 return(constants)
 }

quotaJover2014 <- function(rv){
  #inputs: viral capsid radius rv in micron
  #outputs: quotas for C, N, P [mmol X per indiv]
  rv = rv/1000 #nm
  #Following Jover et al. 2014 https://doi.org/10.1038/nrmicro3289
  x = (10^6/constants()$Avogadro)
  #Quotas from Jover et al. 2014: https://doi.org/10.1038/nrmicro3289
  # all in micromole X per virus
     QC = x*(41*(rv-2.5)^3 + 130*(7.5*rv^2 - 18.75*rv +15.63))
     QN = x*(16*(rv-2.5)^3 + 36*(7.5*rv^2 - 18.75*rv +15.63))
     QP = x*4.2*(rv-2.5)^3
  
return(cbind(QC,QN,QP))
}

dynamicViscosity <- function(Temperature=15, Salinity=35){
	#inputs: Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
	#ouputs: dynamic viscocity of seawater [kg/m-s]
	
	# Following Sharqawy et al. 2010  https://doi.org/10.5004/dwt.2010.1079
  #from MIT seawater:  http://web.mit.edu/seawater/
  # VALIDITY: 0 < T < 180 C and 0 < S < 150 g/kg; ACCURACY: 1.5%

  S = Sal
  T = Temp
  S = S/1000
  a = c(1.5700386464E-01,6.4992620050E+01,-9.1296496657E+01,4.2844324477E-05,1.5409136040E+00,1.9981117208E-02,-9.5203865864E-05,7.9739318223E+00,-7.5614568881E-02,4.7237011074E-04)
  
  mu_w = a[4] + 1/(a[1]*(T+a[2])^2+a[3])
  A  = a[5] + a[6] * T + a[7] * T^2
  B  = a[8] + a[9] * T + a[10]* T^2
  mu = mu_w*(1 + A*S + B*S^2) # [kg/m-s]
return(mu)
}


Stokes_Einstein_Sutherland <- function(radius, Temperature=15, Salinity=35){ #translational diffusion - movement across space
	#inputs: radius as micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
	#output: diffusion as [m^2/s]
	DynamicViscosity = dynamicViscosity(Temperature,Salinity) # [kg/m-s]
	Diff_t = (constants()$Boltzmann *Celcius2Kelvin(Temperature))/( 6*pi*micron_2_metre(radius)*DynamicViscosity) # [m^2/s]
	return( Diff_t  )
}
 
Stokes_Einstein_Debye <- function(radius, Temperature=15, Salinity=35){ #rotational diffusion - angular rotation
  #inputs: radius as ; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
	#output: diffusion as
	DynamicViscosity = dynamicViscosity(Temperature,Salinity) # [kg/m-s]
	Diff_r = (constants()$Boltzmann *Celcius2Kelvin(Temperature))/( 8*pi*radius^3*DynamicViscosity) # [m^2/s]
	return( Diff_r )
 }
 
SmoluchowskiCoagulation <- function(radius1,radius2, Temperature=15, Salinity=35){
	#inputs:  radius1, radius2 in micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
	#ouput: encounter rate via diffusion between spheres of radius1 and radius2 [m^3/s]
	SC = 4*pi*(micron_2_metre(radius1) + micron_2_metre(radius2))*(Stokes_Einstein_Sutherland(radius1,Temperature,Salinity) + Stokes_Einstein_Sutherland(radius2,Temperature,Salinity)) # [m^3/s]
	return( SC )
}

SwimmingSpeed <- function(radius){
	#inputs: radius as micron.  Valid for 1 micron to 1 cm.
	#output: swimming speed as  m per second
	
	#Following Talmy et al. 2019. Front. Mar. Sci. https://doi.org/10.3389/fmars.2019.00182
	u = exp(0.39 + 0.79*log((2*radius)/10^4))/10^2 # [m/s]
	return(u)
}
 
RelativeSwimmingEncounter <- function(radius1,radius2,Temperature=15,Salinity=35){
	#inputs:  radius1, radius2 as micron; Temperature as Celcius; Salinity as g/kg; defined for 0 < T < 180 C and 0 < S < 150 g/kg.
	#output:  encounter rate [m^3/second]
	
	#Following Talmy et al. 2019 https://doi.org/10.3389/fmars.2019.00182 & Visser et al. 2009: https://doi.org/10.1093/plankt/fbn109
	fdetect = 3; #predator detection radius factor (detection = fdetect*r)
	predator = micron_2_metre(max(radius1,radius2))
	prey = micron_2_metre(min(radius1,radius2))
	swim1 = SwimmingSpeed(radius1)
	swim2 = SwimmingSpeed(radius2)
	
  Enc = ((swim1^2+swim2^2)^0.5)*pi*(predator*fdetect + prey)^2  # m^3 per second
	return(Enc)
}




