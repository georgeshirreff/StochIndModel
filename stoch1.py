import random, math

def readParam(param_name,param_file,set):
	f=open(param_file,"r")
	param_value=float('NaN')
	while 1:
		try:
			line = f.next()
			if line.find(param_name) != -1:
				param_values=line.split("\t")
				param_value=float(param_values[set])
		except StopIteration:
			break
	
	print(param_name+"\t"+str(param_value))
	
	return param_value
	
	
def nextStage(timeLastStage,timeDeath,duration):
	if timeLastStage is None or timeLastStage+duration > timeDeath:
		return None
	timeNextStage=timeLastStage+duration
	return timeNextStage
	
def sumInfectivitySex(PersonList,t):
	femaleInfectivity=0
	maleInfectivity=0
	for i in range(0,len(PersonList)):
		thisInfectivity=PersonList[i].getInfectivityCurrent(t)
		if(PersonList[i].sex=="F"):
			femaleInfectivity += thisInfectivity
		else:
			maleInfectivity += thisInfectivity
	return (femaleInfectivity, maleInfectivity)
	
def killIfNecessary(PersonList,t):
	Deaths=0
	for i in reversed(range(0,len(PersonList))):
		if PersonList[i].timeDeath < t:
			Deaths+=1
			del PersonList[i]
	return Deaths
	
def InfectionStep(PersonList,t, dt):
	UninfectedFemales=[]
	UninfectedMales=[]
	InfectedFemales=[]
	InfectedMales=[]
	femaleInfectivity=[]
	maleInfectivity=[]
	for i in range(0,len(PersonList)):
		if PersonList[i].Infected:
			if PersonList[i].sex=="F":
				femaleInfectivity.append(PersonList[i].getInfectivityCurrent(t))
				InfectedFemales.append(i)
			else:
				maleInfectivity.append(PersonList[i].getInfectivityCurrent(t))
				InfectedMales.append(i)
		else:
			if PersonList[i].sex=="F":
				UninfectedFemales.append(i)
			else:
				UninfectedMales.append(i)
	femaleInfectivity_Total=sum(femaleInfectivity)
	maleInfectivity_Total=sum(maleInfectivity)		
	
	n_InfectedFemales=len(InfectedFemales)	
	n_UninfectedFemales=len(UninfectedFemales)	
	n_InfectedMales=len(InfectedMales)	
	n_UninfectedMales=len(UninfectedMales)	

	
	newInfectionsFtoM = len(UninfectedMales)*dt*lambda*lambda*femaleInfectivity_Total
				
	print( "IF"+str(len(InfectedFemales))+" UF"+str(len(UninfectedFemales))+" Finfect"+str(femaleInfectivity_Total)+ " ; IM"+str(len(InfectedMales))+" UM"+str(len(UninfectedMales))+" Minfect"+str(len(maleInfectivity)))
				
	# 	PersonList[0].infect(random.uniform(mu_start,mu_end),t)	
		
	
	
class Person:
	"An individual in the model"
	def __init__(self,lifespan,t):
		self.Infected=False
		self.sex=random.choice(("F","M"))
		self.timeBirth = t
		self.timeDeath = t+random.expovariate(1.0/lifespan)
		print("Sex: "+str(self.sex))
		print("timeBirth: "+str(self.timeBirth))
		print("timeDeath: "+str(self.timeDeath))		
	def infect(self,t,mu_input,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d):
		self.Infected=True
		self.mu=mu_input
		self.timeInfection=t
		print("timeInfection: "+str(t))
		print("mu: "+str(self.mu))
# 	def setViralLoad(self,sigma_e):
		self.ViralLoad= random.normalvariate(self.mu, sigma_e)		
		print("ViralLoad: "+str(self.ViralLoad))
# 	def setTimes(self,dur_p,Dmax,D50,Dk,dur_d,rho):
		self.timeAsymptomatic=nextStage(self.timeInfection,self.timeDeath,random.expovariate(1.0/dur_p))
		print("timeAsymptomatic: "+str(self.timeAsymptomatic))
		
		expectedDurationAsymptomatic=Dmax*pow(D50,Dk)/(pow(10,self.ViralLoad*Dk) + pow(D50,Dk))
		print("expectedDurationAsymptomatic: "+str(expectedDurationAsymptomatic))
		actualDurationAsymptomatic=expectedDurationAsymptomatic*pow(-math.log(random.uniform(0,1)),1.0/rho)/math.gamma(1+1.0/rho)
		
		self.timeDisease=nextStage(self.timeAsymptomatic,self.timeDeath,actualDurationAsymptomatic)
		print("timeDisease: "+str(self.timeDisease))	
			
		self.timeDeathAIDS=nextStage(self.timeDisease,self.timeDeath,random.expovariate(1.0/dur_d))
		if not(self.timeDeathAIDS is None) and self.timeDeathAIDS < self.timeDeath:
			self.timeDeath = self.timeDeathAIDS
		print("timeDeath: "+ str(self.timeDeath))
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
			
		
	
param_file="param_file.txt"

dur_p=readParam("dur_p",param_file,1)
Dmax=readParam("Dmax",param_file,1)
D50=readParam("D50",param_file,1)
Dk=readParam("Dk",param_file,1)
dur_d=readParam("dur_d",param_file,1)
beta_p=readParam("beta_p",param_file,1)
betamax=readParam("betamax",param_file,1)
beta50=readParam("beta50",param_file,1)
betak=readParam("betak",param_file,1)
beta_d=readParam("beta_d",param_file,1)
rho=readParam("rho",param_file,1)
mu_start=readParam("mu_start",param_file,1)
mu_end=readParam("mu_end",param_file,1)
lifespan=readParam("lifespan",param_file,1)
dt=readParam("dt",param_file,1)
start_time=readParam("start_time",param_file,1)
end_time=readParam("end_time",param_file,1)

set=1
pop_size=int(readParam("pop_size",param_file,set))
sigma_e=readParam("sigma_e",param_file,set)
start_mu=readParam("start_mu",param_file,set)


### MAKE THE POPULATION
PersonList=[]
t=start_time
Deaths=0
for i in range(0,pop_size):
	PersonList.append(Person(lifespan,t))

while(t<=end_time):
	DeathsThisTurn=killIfNecessary(PersonList,t)
	if DeathsThisTurn>0:
		Deaths += DeathsThisTurn	
	
	#this infects a single person	
	if(abs(t)<1e-5):
		PersonList[0].infect(t,start_mu,sigma_e,dur_p,Dmax,D50,Dk,dur_d,rho,beta_p,betamax,beta50,betak,beta_d)
	


# 	sumInfectivity=sumInfectivitySex(PersonList,t)
# 	print(sumInfectivity)
	InfectionStep(PersonList,t, dt)
# 	print(str(t)+ " "+str(len(PersonList)) + " " + str(DeathsThisTurn)+" ",str(sumInfectivity[0])+" ",str(sumInfectivity[1]))
	
	t+=dt
#end of while loop	