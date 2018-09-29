import numpy
import random
import pylab

from ps7 import *

''' 
Begin helper code
'''

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. You can use NoChildException as is, you do not need to
    modify/add any code.
    """



class ResistantVirus(SimpleVirus):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb,resistances,mutProb):

    

        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb
        self.resistances = resistances 
        self.mutProb = mutProb


    def isResistantTo(self, drug):
   

        
        return self.resistances.get(drug,False) 
    def getTotalPop(self):
        
        """
        Gets the current total virus population. 
        returns: The total virus population (an integer)
        """

        return len(self.viruses)

    
    def reproduce(self, popDensity,activeDrugs):
 

##checks the resistancy of the mother virus to all drugs in activeDrugs list

      for drugs in activeDrugs:
            if not self.isResistantTo(drugs):
                   raise NoChildException('No child created')
##if mother virus is resistant to all the drugs, reproduction will proceed

           
      maxReproduceProb = self.maxBirthProb * (1 - popDensity)
## to calculate the transfer of resistant property to child virus
        
      if random.random() < maxReproduceProb:
         resistance_trait = {}
         for drugs in self.resistances:
             if random.random() <= (1- self.mutProb):
                resistance_trait[drugs] = self.resistances[drugs]
             else:
                resistance_trait[drugs] = not self.resistances[drugs]
         childOfVirus = ResistantVirus(self.maxBirthProb,self.clearProb,resistance_trait,self.mutProb)
         return childOfVirus
      else:
          raise NoChildException()       
        
     

class Patient(SimplePatient):
    
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """
    
    def __init__(self, viruses, maxPop):
        
    
        self.viruses = viruses
                
    
        self.maxPop = maxPop
    
        self.administeredDrugs = [] 

    



    def addPrescription(self, newDrug):

        self.administeredDrugs.append(newDrug)



    def getPrescriptions(self):

       return self.administeredDrugs



    def getResistPop(self, drugResist):

    ##empty variable to hold the number of resistant virus
       virus_pop = 0 
    ##checks for each virus in viruses of a patient

       for virus in self.viruses:
            drugsResistant = 0
    ##checks for all the drugs administered to the patient

            for drugs in drugResist :
                if virus.isResistantTo(drugs):
                   drugsResistant += 1
   ##if the length of drugsResistant matches with administeredDrugs, then the virus is resistan to  all drugs.
            if drugsResistant == len(drugResist):
               virus_pop  += 1



    def update(self):

        
        popDensity = self.getTotalPop()/float(self.maxPop)
        offspringViruses = []
   
      
        for virus in self.viruses:
            if not virus.doesClear():
               offspringViruses.append(virus)
               try:
                  offspringViruses.append(virus.reproduce(popDensity,self.administeredDrugs))
               except NoChildException: pass
            
        self.viruses = self.viruses + offspringViruses
        print "ok"
        return self.getTotalPop()
      



def virusCollection(numViruses, maxBirthProb, clearProb):
    viruses = []
    for virusNum in range(numViruses):
        viruses.append(SimpleVirus(maxBirthProb, clearProb))
    return viruses


def simulationWithDrug(numTrials = 20, numTimeSteps = 150):

    """
    Run the simulation and plot the graph for problem 2 (no drugs are used,
    viruses do not have any drug resistance).

    Instantiates a patient, runs a simulation for 300 timesteps, and plots the
    total virus population as a function of time.
    """
    random.seed()

    # Virus Characteristics.
    maxPop = 1000
    numViruses = 100
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = { 'guttagonol':True }
    mutProb = 0.005
    dataMatrix1 = numpy.zeros(shape = (numTrials, numTimeSteps))
    dataMatrix2 = numpy.zeros(shape = (numTrials, numTimeSteps))
 
    
    for trial in range(numTrials): 
        viruses = ResistantVirus(maxBirthProb, clearProb,resistances, mutProb) 
        randPatientX = Patient(viruses, maxPop) 
           

         # Simulate the time-steps.
        #dataMatrix1[trial][0] = numViruses
        for time in range(1, numTimeSteps):

  
            dataMatrix1[trial][time]  = randPatientX.update()
            dataMatrix2[trial][time]  = randPatientX.getResistPop(['guttagonol'])
            if time == 150:
               randPatientX.addprescription('guttagonol')
            
            
            

      # Statistical Analysis.
    meanData1 = dataMatrix.mean(0)
    meanData2 = dataMatrix2.mean(0)
    time = numpy.arange(numTimeSteps) 
    stdData95_CI = dataMatrix1.std(0) * 2
    stdData96_CI = dataMatrix2.std(0)*2
    selectedTime = numpy.arange(0, numTimeSteps, 10)

    # Ploting.
    pylab.figure()
    pylab.plot(time, meanData1,color = 'green')
    pylab.plot(time, meanData1,color = 'red')
    #pylab.plot(time,meanData2,color = 'blue')
    pylab.errorbar(time[selectedTime], meanData[selectedTime], stdData95_CI[selectedTime], fmt = 'o')
    #pylab.savefig('prob2_withdrug.png')    
    pylab.show()
   
    
    
simulationWithDrug()
