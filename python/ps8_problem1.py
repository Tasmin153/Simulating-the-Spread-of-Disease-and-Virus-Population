import numpy
import random
import pylab
from ps7 import SimplePatient

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
        
      

        errorMsg1 = 'viruses must be a list containing SimpleVirus objects'
        errorMsg2 = 'maxPop, or maximum virus population must be an integer!'
        
        if type(viruses) != list: raise ValueError(errorMsg1)
        self.viruses = viruses
                
        if type(maxPop)!= int: raise ValueError(errorMsg2)
        self.maxPop = maxPop
    
        self.administeredDrugs = [] 

    def getTotalPop(self):
        
       

        return len(self.viruses)

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

            for drugs in administeredDrugs :
                if virus.isResistant(drugs):
                   drugsResistant += 1
##if the length of drugsResistant matches with administeredDrugs, then the virus is resistan to all drugs.
            if drugsResistant == len(administeredDrugs):
               virus_pop  += 1



    def update(self):


        
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
