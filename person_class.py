import random, math, numpy
#from model_functions import *
#import model_functions
	
class Person:
	"An individual in the model"
	def __init__(self,t,lifespan,grouping):
		self.Infected=False
		self.timeBirth = t
		self.timeDeath = t+lifespan
		self.coreGroup=grouping
#		self.timeDeath = t+random.expovariate(1.0/lifespan)
#		if random.uniform(0,1) < proportionCore:
#			self.coreGroup = True
#		else:
#      			self.coreGroup = False
#		print("timeBirth: "+str(self.timeBirth))
#		print("timeDeath: "+str(self.timeDeath))		
	def infect(self,t,mu_input,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d):
		self.Infected=True
		self.mu=mu_input
		self.timeInfection=t
#		print("timeInfection: "+str(t))
#		print("mu: "+str(self.mu))
		self.n_Infected=0
# 	def setViralLoad(self,sigma_e):
		self.ViralLoad= random.normalvariate(self.mu, sigma_e)		
#		print("ViralLoad: "+str(self.ViralLoad))
# 	def setTimes(self,dur_p,Dmax,D50,Dk,dur_d,rho):
		self.timeAsymptomatic=nextStage(self.timeInfection,self.timeDeath,random.expovariate(1.0/dur_p))
#		print("timeAsymptomatic: "+str(self.timeAsymptomatic))
		
		expectedDurationAsymptomatic=Dmax*pow(D50,Dk)/(pow(10,self.ViralLoad*Dk) + pow(D50,Dk))
#		print("expectedDurationAsymptomatic: "+str(expectedDurationAsymptomatic))
		actualDurationAsymptomatic=expectedDurationAsymptomatic*pow(-math.log(random.uniform(0,1)),1.0/rho)/math.gamma(1+1.0/rho)
		
		self.timeDisease=nextStage(self.timeAsymptomatic,self.timeDeath,actualDurationAsymptomatic)
#		print("timeDisease: "+str(self.timeDisease))	
			
		self.timeDeathAIDS=nextStage(self.timeDisease,self.timeDeath,random.expovariate(1.0/dur_d))
		if not(self.timeDeathAIDS is None) and self.timeDeathAIDS < self.timeDeath:
			self.timeDeath = self.timeDeathAIDS
#		print("timeDeath: "+ str(self.timeDeath))
# 	def setInfectivity(self,beta_p,betamax,beta50,betak,beta_d):
		self.infectivityPrimary= beta_p
		self.infectivityAsymptomatic= betamax*pow(10,self.ViralLoad*betak)/(pow(10,self.ViralLoad*betak)+pow(beta50,betak))
		self.infectivityDisease= beta_d		
	def getInfectivityCurrent(self,t):
		if not(self.Infected):
			return 0
		elif t<self.timeAsymptomatic:
			return self.infectivityPrimary
		elif t<self.timeDisease:
			return self.infectivityAsymptomatic
		elif t<self.timeDeath:
			return self.infectivityDisease
		else:
			print "Error: This individual should have died at time "+str(self.timeDeath)+", but it is now "+str(t)
			quit()
   
def nextStage(timeLastStage,timeDeath,duration): 
	if timeLastStage is None or timeLastStage+duration > timeDeath:
		return None
	timeNextStage=timeLastStage+duration
	return timeNextStage
	
def sumInfectivity(PersonList,t):
	Infectivity=0
	for i in range(0,len(PersonList)):
		thisInfectivity=PersonList[i].getInfectivityCurrent(t)
		Infectivity += thisInfectivity
	return Infectivity
	
def killIfNecessary(PersonList,t): #runs through the PersonList, and if their time of death has passed, deletes them
	DeathsThisTurn=0
	for i in reversed(range(0,len(PersonList))):
		if PersonList[i].timeDeath < t:
			DeathsThisTurn+=1
			del PersonList[i]
	return DeathsThisTurn

def infect_from_to(newInfections,listUninfected,totalInfectivity,n_Infected,listInfectivity,t,PersonList,listInfected,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d):
	ViralLoad_list_addon=[]
	
	if newInfections ==0:
		return ViralLoad_list_addon
	
	
	
	if newInfections > len(listUninfected):
		print "Error: too many infections for the number of individuals.\t"+str(len(listUninfected))+"\t"+str(newInfections)
		exit()
	
	listInfectivity_normalised=[k/totalInfectivity for k in listInfectivity]
	
	infectee_list=random.sample(listUninfected,newInfections)
	infecter_index_list=numpy.random.multinomial(newInfections,listInfectivity_normalised)
	
	infecter_list=[]
	for j in range(0,len(infecter_index_list)):
		infecter_list+=[listInfected[j]]*infecter_index_list[j]
		
		

#test_list=[6,7,8,9,10]
#
#print([0]*10)
#
#test2=numpy.random.multinomial(len(test_list),(0.2,0.2,0.2,0.2,0.2))
#
#print(test2)
#test3=[]
#for i in range(0,len(test2)):
#	z=test_list[i]
#	zn=test2[i]
#	test3=test3+[z]*zn
#
#print(test3)
	
	for i in range(0,newInfections):
		PersonList[infectee_list[i]].infect(t,random.normalvariate(PersonList[infecter_list[i]].mu,sigma_m),sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
		PersonList[infecter_list[i]].n_Infected+=1

		ViralLoad_list_addon.append(PersonList[infectee_list[i]].ViralLoad)
		
#	for i in range(0,newInfections):
#		#pick an infectee
#		infectee=random.choice(listUninfected)
#		listUninfected.remove(infectee)
#            
#		#pick an infector
#		ran_temp=random.uniform(0,totalInfectivity) #pick a random number from zero to total infectivity
#		
#		infectee_list=numpy.random.multinomial(len(listInfected),listInfectivity_cumulative,size=newInfections)
#            
#		for j in range(0,n_Infected):
#			if listInfectivity_cumulative[j] > ran_temp:
#				PersonList[infectee].infect(t,random.normalvariate(PersonList[listInfected[j]].mu,sigma_m),sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
#				PersonList[listInfected[j]].n_Infected += 1
#				break
	return ViralLoad_list_addon
	
def InfectionStep_casual(PersonList,t, dt,report_from,report_every,contact_rate,nonCore_rateMultiplier,assortativity,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d,out_demog,out_final,lastStep):
	listInfectivity_Core=[] #a list of all infectious people's infectivity: Core group
#	listInfectivity_Core_cumulative=[] #a list of infectivity, but cumulative: Core group
	listInfectivity_nonCore=[] #a list of all infectious people's infectivity: nonCore group
#	listInfectivity_nonCore_cumulative=[] #a list of infectivity, but cumulative: nonCore group

	listInfected_Core=[] #a list of the indices of those infected: Core group
	listUninfected_Core=[] #and uninfected: Core group
	listInfected_nonCore=[] #a list of the indices of those infected: nonCore group
	listUninfected_nonCore=[] #and uninfected: nonCore group
	
	
#	print "Get lists"
	ViralLoad_list=[]	
	for i in range(0,len(PersonList)): #runs through the whole PersonList 
		if PersonList[i].Infected:
			temp_Infectivity=PersonList[i].getInfectivityCurrent(t)
			ViralLoad_list.append(PersonList[i].ViralLoad)
			if PersonList[i].coreGroup:
				listInfected_Core.append(i)
				listInfectivity_Core.append(temp_Infectivity) #if they're infected, add their infectivity to the list
#				if len(listInfectivity_Core_cumulative)==0:
#					listInfectivity_Core_cumulative.append(temp_Infectivity) #and to the cumulative list, starting it here if necessary				
#				else:
#					listInfectivity_Core_cumulative.append(temp_Infectivity+listInfectivity_Core_cumulative[-1])
			else:
				listInfected_nonCore.append(i)
				listInfectivity_nonCore.append(temp_Infectivity) #if they're infected, add their infectivity to the list
#				if len(listInfectivity_nonCore_cumulative)==0:
#					listInfectivity_nonCore_cumulative.append(temp_Infectivity) #and to the cumulative list, starting it here if necessary				
#				else:
#					listInfectivity_nonCore_cumulative.append(temp_Infectivity+listInfectivity_nonCore_cumulative[-1])
		else:
			if PersonList[i].coreGroup:
				listUninfected_Core.append(i)
			else:
				listUninfected_nonCore.append(i)
			
	if len(ViralLoad_list)>0:		
		ViralLoad_mean=sum(ViralLoad_list)/float(len(ViralLoad_list))
	else:
		ViralLoad_mean=0

#	print "Admin"
	totalInfectivity_Core=sum(listInfectivity_Core) #total infectivity
	totalInfectivity_nonCore=sum(listInfectivity_nonCore) #total infectivity
	
	n_Infected_Core=len(listInfected_Core)	
	n_Infected_nonCore=len(listInfected_nonCore)	
	n_Uninfected_Core=len(listUninfected_Core)	
	n_Uninfected_nonCore=len(listUninfected_nonCore)	
	n_Core=n_Infected_Core+n_Uninfected_Core
	n_nonCore=n_Infected_nonCore+n_Uninfected_nonCore
	n_Infected=n_Infected_Core+n_Infected_nonCore
	n_Uninfected=n_Uninfected_Core+n_Uninfected_nonCore
 
	p_possibleContacts_CnC=(1-assortativity)/(contact_rate*n_Core+contact_rate*nonCore_rateMultiplier*n_nonCore+1e-100) #the denominator of possible contacts between groups
	p_possibleContacts_CC=p_possibleContacts_CnC+assortativity/(contact_rate*n_Core+1e-100) #the denominator of possible contacts between Core group individuals
	p_possibleContacts_nCnC=p_possibleContacts_CnC+assortativity/(contact_rate*nonCore_rateMultiplier*n_nonCore+1e-100) #the denominator of possible contacts between nonCore group individuals

    #the rate of contact with infected persons, over all possible contacts
	p_newInfections_CC = dt*contact_rate*contact_rate*totalInfectivity_Core*p_possibleContacts_CC 
	p_newInfections_nCnC = dt*contact_rate*nonCore_rateMultiplier*contact_rate*nonCore_rateMultiplier*totalInfectivity_nonCore*p_possibleContacts_nCnC 
	p_newInfections_CnC = dt*contact_rate*contact_rate*nonCore_rateMultiplier*totalInfectivity_nonCore*p_possibleContacts_CnC 
	p_newInfections_nCC = dt*contact_rate*contact_rate*nonCore_rateMultiplier*totalInfectivity_Core*p_possibleContacts_CnC 

#	print "How many infections?"
	(newInfections_CC,newInfections_CnC,C_remainUninfected)=numpy.random.multinomial(n_Uninfected_Core,(p_newInfections_CC,p_newInfections_CnC,0)) #the number of new infections in this step
	(newInfections_nCC,newInfections_nCnC,nC_remainUninfected)=numpy.random.multinomial(n_Uninfected_nonCore,(p_newInfections_nCC,p_newInfections_nCnC,0))

#	newInfections_CC=numpy.random.binomial(n_Uninfected_Core,p_newInfections_CC) #the number of new infections in this step
#	newInfections_CnC=numpy.random.binomial(n_Uninfected_Core,p_newInfections_CnC) #the number of new infections in this step
#	newInfections_nCC=numpy.random.binomial(n_Uninfected_nonCore,p_newInfections_nCC) #the number of new infections in this step
#	newInfections_nCnC=numpy.random.binomial(n_Uninfected_nonCore,p_newInfections_nCnC) #the number of new infections in this step
	
	newInfections=newInfections_CC+newInfections_nCnC+newInfections_nCC+newInfections_CnC
	
#	print(newInfections_CC,newInfections_nCnC,newInfections_nCC,newInfections_CnC)
	
	output_str=str(t)+"\t"+str(n_Uninfected_nonCore)+"\t"+str(n_Uninfected_Core)+"\t"+str(n_Infected_nonCore)+"\t"+str(n_Infected_Core)+"\t"+str(newInfections)+"\t"+str(ViralLoad_mean)+"\n"
	
	if (t>report_from-0.5*dt) and abs(round(t/report_every)-t/report_every) < 1e-5:		
		out_demog.write(output_str)																
		
	if lastStep:
		out_final.write(output_str)																
	
	if newInfections==0:
		return (n_Infected,ViralLoad_mean)
			
#	print "Perform infections"		
	ViralLoad_list=ViralLoad_list+infect_from_to(newInfections_CC,listUninfected_Core,totalInfectivity_Core,n_Infected_Core,listInfectivity_Core,
																t,PersonList,listInfected_Core,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
	ViralLoad_list=ViralLoad_list+infect_from_to(newInfections_CnC,listUninfected_Core,totalInfectivity_nonCore,n_Infected_nonCore,listInfectivity_nonCore,
																t,PersonList,listInfected_nonCore,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
	ViralLoad_list=ViralLoad_list+infect_from_to(newInfections_nCC,listUninfected_nonCore,totalInfectivity_Core,n_Infected_Core,listInfectivity_Core,
																t,PersonList,listInfected_Core,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
	ViralLoad_list=ViralLoad_list+infect_from_to(newInfections_nCnC,listUninfected_nonCore,totalInfectivity_nonCore,n_Infected_nonCore,listInfectivity_nonCore,
																t,PersonList,listInfected_nonCore,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
	
	ViralLoad_mean=sum(ViralLoad_list)/float(len(ViralLoad_list))
	
	return (n_Infected,ViralLoad_mean)
 
# def InfectionStep(PersonList,t, dt,contact_rate,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d):
#	listInfectivity=[] #a list of all infectious people's infectivity
#	listInfectivity_cumulative=[] #a list of infectivity, but cumulative
#	listInfected=[] #a list of the indices of those infected
#	listUninfected=[] #and uninfected
#	for i in range(0,len(PersonList)): #runs through the whole PersonList 
#		if PersonList[i].Infected:      
#			listInfected.append(i)
#			temp_Infectivity=PersonList[i].getInfectivityCurrent(t)
#			listInfectivity.append(temp_Infectivity) #if they're infected, add their infectivity to the list
#			if len(listInfectivity_cumulative)==0: 
#				listInfectivity_cumulative.append(temp_Infectivity) #and to the cumulative list, starting it here if necessary
#			else:
#				listInfectivity_cumulative.append(temp_Infectivity+listInfectivity_cumulative[-1])
#		else:
#			listUninfected.append(i) #otherwise they go in the uninfected list
#	totalInfectivity=sum(listInfectivity) #total infectivity
#	
#	n_Infected=len(listInfected)	
#	n_Uninfected=len(listUninfected)	
#	
#	p_possibleContacts=1.0/(contact_rate*(n_Infected+n_Uninfected)) #the denominator of possible contacts
#	p_newInfections = dt*contact_rate*contact_rate*totalInfectivity*p_possibleContacts #the rate of contact with infected persons, over all possible contacts
#	
#	newInfections=numpy.random.binomial(n_Uninfected,p_newInfections) #the number of new infections in this step
#			
#	print(str(t)+" "+str(n_Uninfected)+" "+str(n_Infected)+" "+str(newInfections))
#		
#	if newInfections==0:
#		return 0
#		
#	temp_newInfections=newInfections
## 	temp_newInfections=5
#	
## 	print(temp_newInfections)
## 	for i in range(0,len(listInfectivity_cumulative)):
## 		listInfectivity_cumulative[i] /= totalInfectivity
## 	print(listInfectivity_cumulative)
## 	listInfectors=numpy.random.multinomial(temp_newInfections,listInfectivity_cumulative)
## 	print(listInfectors)
#	
#	
#	for i in range(0,temp_newInfections):
#		#pick an infectee
#		infectee=random.choice(listUninfected)
#		listUninfected.remove(infectee)
#            
#		#pick an infector
#		ran_temp=random.uniform(0,totalInfectivity) #pick a random number from zero to total infectivity
#            
#		for j in range(0,n_Infected):
#			if listInfectivity_cumulative[j] > ran_temp:
#				PersonList[infectee].infect(t,random.normalvariate(PersonList[listInfected[j]].mu,sigma_m),sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
#				PersonList[listInfected[j]].n_Infected += 1
#				break
#
#	
#	return newInfections
 
	
def BirthStep(PersonList,t,dt,birth_rate,lifeExpectancy,proportionCore):
	newBirths=numpy.random.poisson(lam=dt*birth_rate)
	
	lifespans=numpy.random.exponential(lifeExpectancy,newBirths)
	groupings=numpy.random.binomial(1,proportionCore,newBirths)
		
	for i in range(0,newBirths):
		PersonList.append(Person(t,lifespans[i],groupings[i]))
			
	return newBirths	
    
