
cubicmetres_2_Litres <- function(X){
  return(X*1000.0)
}
Litres_2_cubicmetres <-function(X){
  return(X/1000.0)
}

days_2_seconds <- function(X){
return(X*86400.0)
}
seconds_2_days <- function(X){
return(X/86400.0)
}

nm_2_micron <- function(X){
return(X/1000.0)
}
micron_2_nm <- function(X){
return(X*1000.0)
}

metre_2_micron <-function(X){
return(X*10.0^6.0)
}
micron_2_metre <-function(X){
return(X/(10.0^6.0))
}

cubicmicron_2_Litres<-function(X){
return(X*(10.0^-15.0))
}
Litres_2_cubicmicron<-function(X){
return(X*(10.0^15.0))
}

mole_2_micromole<-function(X){
return(X*10.0^6.0)
}
micromole_2_mole<-function(X){
return(X*10.0^-6.0)
}

radius_2_volume <- function(X){
return( 4.0/3.0 * pi * X^3.0)
}
volume_2_radius<-function(X){
return( (3.0*X/(4.0*pi))^(1.0/3.0) )
}

IndivPerL_2_BiomassmmolPerCubicmetre<-function(X,Q){
#  [#ind/L]  ->  [mmol Q m^-3]
return( Litres_2_cubicmetres(X*Q) )
}

BiomassmmolPerCubicmetre_2_IndivPerL<-function(X,Q){
#  [mmol Q m^-3] -> [#ind/L]
return( cubicmetres_2_Litres(X/Q) )
}
