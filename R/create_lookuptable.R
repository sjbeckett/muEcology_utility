

referenceQ = c(1.7*10^-7.0, 6.1*10^-10.0 , 10^-11.0, 10^-11.0)*1000
referenceName = c("diatom biomass", "Pro.  biomass", "diatom virus  biomass", "Pro. virus  biomass")

#some reference values in indiv/L
refvals_perL = rev(10^c(0:9))

#create table
mmolcubicm_MAT = matrix(NA,length(refvals_perL),length(referenceQ))

for(aa in 1:length(refvals_perL)){
	for(bb in 1:length(referenceQ)){
		mmolcubicm_MAT[aa,bb] = IndivPerL_2_BiomassmmolPerCubicmetre(refvals_perL[aa],referenceQ[bb])
	}
}

#put into tabular format
JJ = cbind(refvals_perL,mmolcubicm_MAT)
JJ = as.data.frame(JJ)
names(JJ) = c("reference [ indiv/L ]", paste("[mmol N m^-3]",referenceName))


#write to file
write.csv(JJ,file = "lookuptable_indivperL_mmolNpercubicm.csv")
