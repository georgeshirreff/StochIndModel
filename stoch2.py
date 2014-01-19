import random, math, numpy, scipy
from basic_functions import *
#from model_functions import *
from person_class import *
# from numpy import *
# from random import *
# from math import *
# from basic_functions import *
# from person_class import *

#rng_seed=2
#random.seed(rng_seed)
#numpy.random.seed(rng_seed)

param_file="param_file.txt"

extra_name=readString("extra_name",param_file,1)
out_info=open(str(extra_name)+"_info.txt","w") #opens this file and empties it
#out_info=close()
#out_info=open(str(extra_name)+"_info.txt","a")

write_info=bool(readParam("write_info",param_file,1,out_info))
write_demog=bool(readParam("write_demog",param_file,1,out_info))
write_epi=bool(readParam("write_epi",param_file,1,out_info))
write_final=bool(readParam("write_final",param_file,1,out_info))
Nsets=int(readParam("Nsets",param_file,1,out_info))
Nreps=int(readParam("Nreps",param_file,1,out_info))

dur_p=readParam("dur_p",param_file,1,out_info)
Dmax=readParam("Dmax",param_file,1,out_info)
D50=readParam("D50",param_file,1,out_info)
Dk=readParam("Dk",param_file,1,out_info)
dur_d=readParam("dur_d",param_file,1,out_info)
beta_p=readParam("beta_p",param_file,1,out_info)
betamax=readParam("betamax",param_file,1,out_info)
beta50=readParam("beta50",param_file,1,out_info)
betak=readParam("betak",param_file,1,out_info)
beta_d=readParam("beta_d",param_file,1,out_info)
rho=readParam("rho",param_file,1,out_info)
#mu_start=readParam("mu_start",param_file,1)
#mu_end=readParam("mu_end",param_file,1)
lifeExpectancy=readParam("lifeExpectancy",param_file,1,out_info)
dt=readParam("dt",param_file,1,out_info)
start_time=readParam("start_time",param_file,1,out_info)
report_from=readParam("report_from",param_file,1,out_info)
report_every=readParam("report_every",param_file,1,out_info)
end_time=readParam("end_time",param_file,1,out_info)

out_info=close() 
if write_demog | write_epi | write_final: #copy the permanent parameters into any other files which are for output
	out_info=open(str(extra_name)+"_info.txt","r")
	contents=out_info.read()
	if write_demog:
		out_demog=open(str(extra_name)+"_demog.txt","w")
		out_demog.write(contents)
	else:
		out_demog=-1
	if write_final:
		out_final=open(str(extra_name)+"_final.txt","w")	
		out_final.write(contents)
	else:
		out_final=-1
	if write_epi:
		out_epi=open(str(extra_name)+"_epi.txt","w")	
		out_epi.write(contents)
	else:
		out_epi=-1
		
outputHeader="t\tn_Uninfected_nonCore\tn_Uninfected_Core\tn_Infected_nonCore\tn_Infected_Core\tnewInfections\tViralLoad_mean\n"		
		
for set in range(1,Nsets+1):
	print("set\t"+str(set))
	out_info=open(str(extra_name)+"_info.txt","a") #open the info file again for output of changeable parameters
	out_info.write("set\t"+str(set)+"\n") 
	
	pop_size=int(readParam("pop_size",param_file,set,out_info))
	sigma_e=readParam("sigma_e",param_file,set,out_info)
	sigma_m=readParam("sigma_m",param_file,set,out_info)
	start_mu=readParam("start_mu",param_file,set,out_info)
	contact_rate=readParam("contact_rate",param_file,set,out_info)
	proportionCore=readParam("proportionCore",param_file,set,out_info)
	nonCore_rateMultiplier=readParam("nonCore_rateMultiplier",param_file,set,out_info)
	assortativity=readParam("assortativity",param_file,set,out_info)
	n_founders=readParam("n_founders",param_file,set,out_info)	
	
	out_info=close() #copy from the last instance of "set" to the end of the info file, into any other relevant output files
	if write_demog | write_epi | write_final:
		out_info=open(str(extra_name)+"_info.txt","r")
		contents=out_info.read()
		last_set=contents.rfind("set")
		end_of_file=len(contents)
		changeableParams=contents[last_set:end_of_file]
		if write_demog:
			out_demog=open(str(extra_name)+"_demog.txt","a")
			out_demog.write(changeableParams)
		if write_final:
			out_final=open(str(extra_name)+"_final.txt","a")	
			out_final.write(changeableParams)	
			#in this case the "rep" goes in the header line
			out_final.write("rep\t"+outputHeader) 
		if write_epi:
			out_epi=open(str(extra_name)+"_epi.txt","a")	
#			out_epi.write(changeableParams)
			changeableParams_list=changeableParams.split("\n")
			changeableParams_list=changeableParams_list[:-1] #takes out the last line, whichis empty				
			if set==1: #only on the first set				
				for p in range(0,len(changeableParams_list)): #writes the header line
					paramPair=changeableParams_list[p].split("\t")
					out_epi.write(paramPair[0]+"\t")
	#				out_epi.write((changeableParams_list[p].split("\t"))[1]+"\t")
				out_epi.write("p_epi\tp_epi_lo\tp_epi_hi\tmax_before_extinction\tR0\tViralLoad_mean_mean\n") #then the rest of the header line				
			for p in range(0,len(changeableParams_list)): #then write the changeable parameter values into the next line, do this for every set
				paramPair=changeableParams_list[p].split("\t")
				out_epi.write(paramPair[1]+"\t")
	
	birth_rate=1.0/lifeExpectancy*pop_size
	
	max_Infected_extinct=0
	n_successes=0.0
	ViralLoad_mean_mean=0.0
	R0=5
	for rep in range(1,Nreps+1):
#		print("rep\t"+str(rep))
		print(rep)
		
		if write_demog:
			out_demog.write("rep\t"+str(rep)+"\n")
			out_demog.write(outputHeader)																
		if write_final:
			out_final.write(str(rep)+"\t")
		
		### MAKE THE POPULATION
		PersonList=[]
		t=start_time
		Deaths=0
		
		lifespans=numpy.random.exponential(lifeExpectancy,pop_size)
		groupings=numpy.random.binomial(1,proportionCore,pop_size)
		
		for i in range(0,pop_size):
			PersonList.append(Person(t,lifespans[i],bool(groupings[i])))
		
		n_Infected=0
		lastStep=False
		postInfection=False
		max_Infected=0
		while t < end_time+0.5*dt:
#			print t
#			print "Doing births"
			DeathsThisTurn=killIfNecessary(PersonList,t)
			Deaths += DeathsThisTurn
			BirthsThisTurn=BirthStep(PersonList,t,dt,birth_rate,lifeExpectancy,proportionCore)
						
#			print "Doing extinction"
			extinct= (postInfection and n_Infected==0) or len(PersonList)==0			
			if extinct or not(t+dt < end_time+0.5*dt):
				lastStep=True
				
#			print "Doing infection"				
			if(abs(t)<1e-5):
				postInfection=True
				for i in range(0,n_founders):
					PersonList[i].infect(t,start_mu,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d) #infects patient zero and any other founders if appropriate
			#elif postInfection:	
			
			(n_Infected,ViralLoad_mean)=InfectionStep_casual(PersonList,t, dt,report_from,report_every,contact_rate,nonCore_rateMultiplier,assortativity,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d,out_demog,out_final,lastStep)
			max_Infected=max(max_Infected,n_Infected)
		
			#	this infects a single person
		
			if lastStep: # | len(PersonList)==0:
				break

			t+=dt
		#end of time loop
			
		if extinct:
			print(t)
			max_Infected_extinct=max(max_Infected_extinct,max_Infected)
		else:
			ViralLoad_mean_mean+=ViralLoad_mean
			n_successes+=1.0
			
	#end of rep loop
	p_epi=float(n_successes)/float(Nreps)
	if n_successes==0:
		p_epi_lo=0
	else:
		p_epi_lo=1-scipy.stats.beta.ppf(0.975,Nreps-n_successes+1,n_successes)
	if n_successes==1:
		p_epi_hi=1
	else:
		p_epi_hi=1-scipy.stats.beta.ppf(0.975,Nreps-n_successes,n_successes+1)
		
	ViralLoad_mean_mean/=float(n_successes)
	out_epi.write(str(p_epi)+"\t"+str(p_epi_lo)+"\t"+str(p_epi_hi)+"\t"+str(max_Infected_extinct)+"\t"+str(R0)+"\t"+str(ViralLoad_mean_mean)+"\n")
	
	out_demog.close()
 	out_final.close()
	out_epi.close()

