import random, math, numpy


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
	Deaths=0
	for i in reversed(range(0,len(PersonList))):
		if PersonList[i].timeDeath < t:
			Deaths+=1
			del PersonList[i]
	return Deaths
	
def InfectionStep(PersonList,t, dt,contact_rate,sigma_m,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d):
	listInfectivity=[] #a list of all infectious people's infectivity
	listInfectivity_cumulative=[] #a list of infectivity, but cumulative
	listInfected=[] #a list of the indices of those infected
	listUninfected=[] #and uninfected
	for i in range(0,len(PersonList)): #runs through the whole PersonList 
		if PersonList[i].Infected:      
			listInfected.append(i)
			temp_Infectivity=PersonList[i].getInfectivityCurrent(t)
			listInfectivity.append(temp_Infectivity) #if they're infected, add their infectivity to the list
			if len(listInfectivity_cumulative)==0: 
				listInfectivity_cumulative.append(temp_Infectivity) #and to the cumulative list, starting it here if necessary
			else:
				listInfectivity_cumulative.append(temp_Infectivity+listInfectivity_cumulative[-1])
		else:
			listUninfected.append(i) #otherwise they go in the uninfected list
	totalInfectivity=sum(listInfectivity) #total infectivity
	
	n_Infected=len(listInfected)	
	n_Uninfected=len(listUninfected)	
	
	p_possibleContacts=1.0/(contact_rate*(n_Infected+n_Uninfected)) #the denominator of possible contacts
	p_newInfections = dt*contact_rate*contact_rate*totalInfectivity*p_possibleContacts #the rate of contact with infected persons, over all possible contacts
	
	newInfections=numpy.random.binomial(n_Uninfected,p_newInfections) #the number of new infections in this step
			
	print(str(t)+" "+str(n_Uninfected)+" "+str(n_Infected)+" "+str(newInfections))
		
	if newInfections==0:
		return 0
		
	temp_newInfections=newInfections
# 	temp_newInfections=5
	
# 	print(temp_newInfections)
# 	for i in range(0,len(listInfectivity_cumulative)):
# 		listInfectivity_cumulative[i] /= totalInfectivity
# 	print(listInfectivity_cumulative)
# 	listInfectors=numpy.random.multinomial(temp_newInfections,listInfectivity_cumulative)
# 	print(listInfectors)
	
	
	for i in range(0,temp_newInfections):
		#pick an infectee
		infectee=random.choice(listUninfected)
		listUninfected.remove(infectee)
            
		#pick an infector
		ran_temp=random.uniform(0,totalInfectivity) #pick a random number from zero to total infectivity
            
		for j in range(0,n_Infected):
			if listInfectivity_cumulative[j] > ran_temp:
				PersonList[infectee].infect(t,random.normalvariate(PersonList[listInfected[j]].mu,sigma_m),sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
				PersonList[listInfected[j]].n_Infected += 1
				break

	
	return newInfections
 
#import person_class #why does this have to go here, and can't just go above it

