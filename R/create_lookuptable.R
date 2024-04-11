
#reference individual level quotas (from Paul) in terms of [mmol N/indiv]
referenceQ = c(1.7*10^-7.0, 6.1*10^-10.0 , 10^-11.0, 10^-11.0)
referenceName = c("diatom cells", "Pro. cells", "diatom viruses", "Pro. viruses")

#some reference values in L/day to encompass organism densities and contact rates
refvalsLperday = rev(10^c(-15:9))

#create table
mmolcubicm_MAT = matrix(NA,length(refvalsLperday),length(referenceQ))

for(aa in 1:length(refvalsLperday)){
	for(bb in 1:length(referenceQ)){
		mmolcubicm_MAT[aa,bb] = IndivPerL_2_BiomassmmolPerCubicmetre(refvalsLperday[aa],referenceQ[bb])
	}
}

#put into tabular format
JJ = cbind(refvalsLperday,mmolcubicm_MAT)
JJ = as.data.frame(JJ)
names(JJ) = c("reference [ L/day ]", paste("[mmol N m^-3]",referenceName))

#write to file
write.csv(JJ,file = "lookuptable_Lday_mmolcubicm.csv")
