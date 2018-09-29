import numpy
import random
import pylab

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

'''
End helper code
'''

#
# PROBLEM 1
#
class ResistantVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb,resistances,mutProb):

        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        resistances: A dictionary of drug names (strings) mapping to the
        state of this virus particle's resistance (either True or False) to each
        drug.e.g. {'guttagonol':False, 'grimpex',False}, means that this virus
        particle is resistant to neither guttagonol nor grimpex.
        mutProb: Mutation probability for this virus particle (a float). This
        is the probability of the offspring acquiring or losing resistance to a
        drug
        """

        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb
        self.resistances = resistances 
        self.mutProb = mutProb


    def isResistantTo(self, drug):
        """
Get the state of this virus particle's resistance to a drug. This
method
is called by getResistPop() in Patient to determine how many virus
particles have resistance to a drug.
drug: The drug (a string)
returns: True if this virus instance is resistant to the drug, False
otherwise.
"""

        
        return self.resistances[drug] 


    def doesClear(self):

        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.clearProb and otherwise returns
        False.
        """        
        return random.random() < self.clearProb

    
    def reproduce(self, popDensity,activeDrugs):
        """
Stochastically determines whether this virus particle reproduces at a
time step. Called by the update() method in the Patient class.
A virus particle will only reproduce if it is resistant to ALL the
drugs
in the activeDrugs list. For example, if there are 2 drugs in the
activeDrugs list, and the virus particle is resistant to 1 or no
drugs,
then it will NOT reproduce.
Hence, if the virus is resistant to all drugs
in activeDrugs, then the virus reproduces with probability:
self.maxBirthProb * (1 - popDensity).If this virus particle reproduces, then reproduce() creates and
returns
the instance of the offspring ResistantVirus (which has the same
maxBirthProb and clearProb values as its parent).
For each drug resistance trait of the virus (i.e. each key of
self.resistances), the offspring has probability 1-mutProb of
inheriting that resistance trait from the parent, and probability
mutProb of switching that resistance trait in the offspring.
For example, if a virus particle is resistant to guttagonol but not
grimpex, and `self.mutProb` is 0.1, then there is a 10% chance that
that the offspring will lose resistance to guttagonol and a 90%
chance that the offspring will be resistant to guttagonol.
There is also a 10% chance that the offspring will gain resistance to
grimpex and a 90% chance that the offspring will not be resistant to
grimpex.
popDensity: the population density (a float), defined as the current
virus population divided by the maximum population
activeDrugs: a list of the drug names acting on this virus particle
(a list of strings).
returns: a new instance of the ResistantVirus class representing the
offspring of this virus particle. The child should have the same
maxBirthProb and clearProb values as this virus. Raises a
NoChildException if this virus particle does not reproduce.
"""

        # Does the virus reproduce? 
##checks the resistancy of the mother virus to all drugs in activeDrugs list

        resistancy = []
        for drugs in activeDrugs:
            resistancy.append (self.resistances[drugs])
##if mother virus is resistant to all the drugs, reproduction will proceed

        if False not in resistancy:      
           maxReproduceProb = self.maxBirthProb * (1 - popDensity)
## to calculate the transfer of resistant property to child virus
        
           if random.random() < maxReproduceProb:
              resistance_trait = {}
              for drugs in self.resistances:
                  if random.random() <= (1- self.mutProb):
                     resistance_trait[drug] = self.resistances[drug]
                  else:
                     resistance_trait[drug] = not self.resistances[drug]
              childOfVirus = ResistantVirus(self.maxBirthProb,resistance_trait, self.clearProb)
              return childOfVirus
             
        
           else: raise NoChildException('Child not created!')

class Patient(SimplePatient):
    
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """
    
    def __init__(self, viruses, maxPop):
        
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the  maximum virus population for this patient (an integer)
        """

        errorMsg1 = 'viruses must be a list containing SimpleVirus objects'
        errorMsg2 = 'maxPop, or maximum virus population must be an integer!'
        
        if type(viruses) != list: raise ValueError(errorMsg1)
        self.viruses = viruses
                
        if type(maxPop)!= int: raise ValueError(errorMsg2)
        self.maxPop = maxPop
    
        self.administeredDrugs = [] 

    def getTotalPop(self):
        
        """
        Gets the current total virus population. 
        returns: The total virus population (an integer)
        """

        return len(self.viruses)

    def addPrescription(self, newDrug):
"""
Administer a drug to this patient. After a prescription is added, the
drug acts on the virus population for all subsequent time steps. If
the
newDrug is already prescribed to this patient, the method has no
effect.
newDrug: The name of the drug to administer to the patient (a
string).
postcondition: The list of drugs being administered to a patient is
updated
"""
        self.administeredDrugs.append(newDrug)

    def getPrescriptions(self):
"""
Returns the drugs that are being administered to this patient.
returns: The list of drug names (strings) being administered to this
patient.
"""
       return self.administeredDrugs

    def getResistPop(self, drugResist):
"""
Get the population of virus particles resistant to the drugs listed
in
drugResist.
drugResist: Which drug resistances to include in the population (a
list
of strings - e.g. ['guttagonol'] or ['guttagonol', 'grimpex'])
returns: The population of viruses (an integer) with resistances to
all
drugs in the drugResist list.
 
"""    ##empty variable to hold the number of resistant virus
       virus_pop = 0 
##checks for each virus in viruses of a patient

       for virus in self.viruses:
            drugsResistant = 0
##checks for all the drugs administered to the patient

            for drugs in administeredDrugs :
                if virus.isResistant(drugs):
                   drugsResistant += 1
##if the length of drugsResistant matches with administeredDrugs, then the virus is resistan to all drugs.
            if drugsResistant == len(administeredDrugs):
               virus_pop  += 1



    def update(self):
"""
Update the state of the virus population in this patient for a single
time step. update() should execute these actions in order:
- Determine whether each virus particle survives and update the list
ofvirus particles accordingly
- The current population density is calculated. This population
density
value is used until the next call to update().
- Determine whether each virus particle should reproduce and add
offspring virus particles to the list of viruses in this patient.
The listof drugs being administered should be accounted for in the
determination of whether each virus particle reproduces.
returns: The total virus population at the end of the update (an
integer)
"""
# TODO
        
        # Determine number of viruses to be cleaned, "stochastically".
        numRemoveVirus = 0
        for virus in self.viruses:
            if virus.doesClear():
                numRemoveVirus += 1

        # Remove numRemoveVirus from the patient's body.
        for virusNum in range(numRemoveVirus):
            self.viruses.pop()

        # Calculate population density. TO DO check (keep self!)
        popDensity = self.getTotalPop()/float(self.maxPop)
        
        if popDensity >= 1:
            print 'virus population reached maximum!'
            popDensity = 1       

        # Reproduce at a single time step.
        offspringViruses = []
        for virus in self.viruses:
            try:
                offspringViruses.append(virus.reproduce(popDensity),self.administeredDrugs)
            except NoChildException: pass
            
        self.viruses = self.viruses + offspringViruses
        
        return self.getTotalPop()

def virusCollection(numViruses, maxBirthProb, clearProb):
    viruses = []
    for virusNum in range(numViruses):
        viruses.append(SimpleVirus(maxBirthProb, clearProb))
    return viruses       

#
# PROBLEM 2
#
def simulationWithoutDrug(numTrials = 100, numTimeSteps = 500):

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
    
    dataMatrix = numpy.zeros(shape = (numTrials, numTimeSteps))    
    for trial in range(numTrials):        

        # Model a random patient with the given virus charateristics.        
        viruses = virusCollection(numViruses, maxBirthProb, clearProb)
        randPatientX = SimplePatient(viruses, maxPop)

        # Simulate the time-steps.
        dataMatrix[trial][0] = numViruses
        for time in range(1, numTimeSteps):
            dataMatrix[trial][time] = randPatientX.update()           
            
    # Statistical Analysis.
    meanData = dataMatrix.mean(0)
    time = numpy.arange(numTimeSteps) 
    stdData95_CI = dataMatrix.std(0) * 2
    selectedTime = numpy.arange(0, numTimeSteps, 10)

    # Ploting.
    pylab.plot(time, meanData)
    pylab.errorbar(time[selectedTime], meanData[selectedTime], stdData95_CI[selectedTime], fmt = 'o')    
    pylab.show()
    
simulationWithoutDrug()

# End: August 4, 2016; 9:51 am
