# Problem Set 7SC: Virus Population dynamics modeling.
# Example Problem: Simulating the Spread of Disease and Virus Population Dynamics 
# Name: Mohammad Ehsanul Karim
# Collaborators: None
# Start: July 31, 2016; 2:24 pm

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
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):

        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """

        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb

    def doesClear(self):

        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.clearProb and otherwise returns
        False.
        """        
        return random.random() < self.clearProb

    
    def reproduce(self, popDensity):

        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the SimplePatient and
        Patient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """

        # Does the virus reproduce?        
        maxReproduceProb = self.maxBirthProb * (1 - popDensity)
        
        if random.random() < maxReproduceProb:
            childOfVirus = SimpleVirus(self.maxBirthProb, self.clearProb)
            return childOfVirus
        
        else: raise NoChildException('Child not created!')

class SimplePatient(object):
    
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

    def getTotalPop(self):
        
        """
        Gets the current total virus population. 
        returns: The total virus population (an integer)
        """

        return len(self.viruses)

    def update(self):

        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        - Determine whether each virus particle should reproduce and add
          offspring virus particles to the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        
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
                offspringViruses.append(virus.reproduce(popDensity))
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
