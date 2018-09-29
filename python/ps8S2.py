import numpy
import random
import pylab
import collections
import matplotlib.mlab as mlab

from ps7_v1 import *

class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """
    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as
        attributes of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'grimpex',False}, means that this virus
        particle is resistant to neither guttagonol nor grimpex.

        mutProb: Mutation probability for this virus particle (a float). This
        is the probability of the offspring acquiring or losing resistance to a drug.
        """
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb
        self.resistances = resistances
        self.mutProb = mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This
        method is called by getResistPop() in Patient to determine how many virus
        particles have resistance to a drug.

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug,
        False otherwise.
        """
        # Assuming the virus is resistant to an unknown drug
        return self.resistances.get(drug, True)

    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient class.

        A virus particle will only reproduce if it is resistant to ALL the
        drugs in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no
        drugs, then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:
        self.maxBirthProb * (1 - popDensity).

        If this virus particle reproduces, then reproduce() creates and
        returns the instance of the offspring ResistantVirus (which has the same
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

        # Checks the resistance of the mother virus to all drugs in activeDrugs list
        # if mother virus is resistant to all the drugs, reproduction will proceed
        for drug in activeDrugs:
            if not self.isResistantTo(drug):
                raise NoChildException()

        maxReproduceProb = self.maxBirthProb * (1 - popDensity)

        if random.random() < maxReproduceProb:
            resistance_trait = {}
            # Calculate the transfer of resistances property to child virus
            for drug in self.resistances.keys():
                if random.random() <= (1 - self.mutProb):
                    resistance_trait[drug] = self.resistances[drug]
                else:
                    resistance_trait[drug] = not self.resistances[drug]

            childOfVirus = ResistantVirus(self.maxBirthProb, self.clearProb, resistance_trait, self.mutProb)
            return childOfVirus
        else:
            raise NoChildException()

class Patient(SimplePatient):
    """
    Representation of a patient. The patient is able to take drugs and
    his/her virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).

        viruses: The list representing the virus population (a list of
        ResistantVirus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """

        self.viruses = viruses
        self.maxPop = maxPop
        self.administeredDrugs = []

    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If
        the newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """

        self.administeredDrugs.append(newDrug)

    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this patient.
        """

        return self.administeredDrugs

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed
        in drugResist.
 
        drugResist: Which drug resistances to include in the population (a
        list of strings - e.g. ['guttagonol'] or ['guttagonol', 'grimpex'])
 
        returns: The population of viruses (an integer) with resistances to
        all drugs in the drugResist list.
        """

        # Variable to hold the number of resistant viruses
        virus_pop = 0

        # Checks each virus in viruses of a patient
        for virus in self.viruses:
            drugsResisted = 0
            for drug in drugResist:
                if virus.isResistantTo(drug):
                    drugsResisted += 1
            # If the drugsResisted equals length of drugResist,
            # then the virus is resistant to all drugs in drugResist
            if len(drugResist) != 0 and drugsResisted == len(drugResist):
                virus_pop += 1
        
        return virus_pop

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list
          virus particles accordingly

        - The current population density is calculated. This population
          density value is used until the next call to update().

        - Determine whether each virus particle should reproduce and add
          offspring virus particles to the list of viruses in this patient.
          The listof drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an integer)
        """

        # Determine how many viruses survive
        for virus in self.viruses:
            if virus.doesClear():
                self.viruses.remove(virus)

        popDensity = self.getTotalPop() / float(self.maxPop)

        offspringViruses = []
        for virus in self.viruses:
            try:
                offspringViruses.append(virus.reproduce(popDensity, self.administeredDrugs))
            except NoChildException:
                pass
        self.viruses = self.viruses + offspringViruses

        return self.getTotalPop()

def resistantVirusCollection(numViruses, maxBirthProb, clearProb, resistances, mutProb):
    '''
    Creates a list of ResistantVirus particles
    '''
    viruses = []
    for virusNum in range(numViruses):
        viruses.append(ResistantVirus(maxBirthProb, clearProb, resistances, mutProb))
    return viruses


def simulationWithDrug(numTrials=10, numTimeStepsBeforeDrug=150, numTimeStepsAfterDrug=150, drug='guttagonol'):
    '''
    Runs simulations when only one drug is administered to a 
    Patient with ResistantVirus particles.

    The simulation consists of consists of 150 time steps, 
    followed by the addition of the drug, 
    guttagonol, followed by another 150 time steps.

    Returns the final total populations from all of the trials
    '''

    # To get same random results
    random.seed(0)

    maxPop = 1000 # maximum sustainable virus population
    numViruses = 100 # initial number of viruses
    numTimeSteps = numTimeStepsBeforeDrug + numTimeStepsAfterDrug

    # Virus Characteristics
    maxBirthProb = 0.1
    clearProb = 0.05
    resistances = {'guttagonol': False}
    mutProb = 0.005

    dataMatrixTotal = numpy.zeros(shape=(numTrials, numTimeSteps))
    dataMatrixDrugResistant = numpy.zeros(shape=(numTrials, numTimeSteps))

    for trial in range(numTrials):
        # print "trial =", trial
        # Model a random patient with the given virus charateristics
        viruses = resistantVirusCollection(numViruses, maxBirthProb, clearProb, resistances, mutProb)
        randPatientX = Patient(viruses, maxPop)

        # Simulate the time-steps
        for time in range(numTimeSteps):
            if time == numTimeStepsBeforeDrug:
                randPatientX.addPrescription(drug)
            if time == 0:
                dataMatrixTotal[trial][0] = randPatientX.getTotalPop()
                dataMatrixDrugResistant[trial][0] = randPatientX.getResistPop([drug])
            else:
                dataMatrixTotal[trial][time] = randPatientX.update()
                dataMatrixDrugResistant[trial][time] = randPatientX.getResistPop([drug])
    
    finalTotalPopulations = list(dataMatrixTotal[:,-1])
    
    # Statistical Analysis.
    meanDataTotal = dataMatrixTotal.mean(0)
    meanDataDrugResistant = dataMatrixDrugResistant.mean(0)
    time = numpy.arange(numTimeSteps)
    stdDataTotal = dataMatrixTotal.std(0) * 2
    stdDataDrugResistant = dataMatrixDrugResistant.std(0) * 2
    selectedTime = numpy.arange(0, numTimeSteps, 10)

    pylab.clf()
    pylab.figure(num=0, figsize=(12, 8), dpi=72, facecolor='w', edgecolor='k')
    # Plot average total population
    totalPop, = pylab.plot(time, meanDataTotal, color='green', label='Average Total Virus Population')
    pylab.errorbar(time[selectedTime], meanDataTotal[selectedTime], \
        stdDataTotal[selectedTime], fmt='go')
    # Plot average drug resistant population
    drugResistantPop, = pylab.plot(time, meanDataDrugResistant, color='red', \
        label='Average '+ drug.capitalize() + ' Resistant Virus Population')
    pylab.errorbar(time[selectedTime], meanDataDrugResistant[selectedTime], \
        stdDataDrugResistant[selectedTime], fmt='ro')
    pylab.legend(loc = 1)
    pylab.xlabel("Number of Time Steps")
    pylab.ylabel("Average Virus Population")
    pylab.title("Virus Population Dynamics Simulation with One Drug")
    # pylab.show()
    figure_name = drug + "_nt_" + str(numTrials) + "_ntsbd_" + str(numTimeStepsBeforeDrug) + ".png"
    pylab.savefig(figure_name)

    return finalTotalPopulations

def simulationDelayedTreatment(numTrials=10, numTimeStepsBeforeDrug=150, \
    numTimeStepsAfterDrug=150, drug='guttagonol'):
    '''
    Runs the simulation for delayed treatment by the
    admninistration of only one drug.

    Returns the percentage of the cured patients on the
    assumption that a final virus population of 0-50 in
    a patient means they are cured.
    '''

    finalTotalPopulations = simulationWithDrug(numTrials, numTimeStepsBeforeDrug, numTimeStepsAfterDrug, drug)

    # Find out the percentage of cured patients
    curedPatientsVirusPopulations = [population for population in finalTotalPopulations if population <= 50]
    curedPatientPercentage = len(curedPatientsVirusPopulations) / float(len(finalTotalPopulations)) * 100

    ## Make histogram
    pylab.clf()
    mu = numpy.mean(finalTotalPopulations)
    # The histogram
    n, bins, patches = pylab.hist(finalTotalPopulations, bins=20, range=(0, 1000), facecolor='green', histtype='bar', alpha=0.75)
    pylab.xlabel('Final Virus Population')
    pylab.ylabel('Number of Patients')
    pylab.title(r'$\mathrm{Histogram\ of\ Final\ Virus\ Population\ Distribution\ for\ Drug\ Admininstation\ After\ %d\ Steps:}\  \mu=%.2f$' \
        % (numTimeStepsBeforeDrug, mu))
    pylab.xlim([0, 1000])
    pylab.ylim([0, numTrials])
    pylab.xticks(numpy.arange(0, 1001, 50))
    fileName = "hist_" + drug + "_nt_" + str(numTrials) + "_ntsbd_" + str(numTimeStepsBeforeDrug) + "_steps.png"
    pylab.savefig(fileName)

    return curedPatientPercentage

if __name__ == '__main__':
    # problems 2 and 3
    trials = 500
    BEFORE_DRUG_STEPS = [300, 150, 75, 0]

    f_name = str(trials) + "_cured_percentages.txt"
    f = open(f_name, mode='w')
    for steps in BEFORE_DRUG_STEPS:
       curedPatientPercentage = simulationDelayedTreatment(numTrials=trials, numTimeStepsBeforeDrug=steps)
       stirngToWrite = "Percentage of cured patients after {:d} steps: {:.2f}\n".format(steps, curedPatientPercentage) 
       f.write(stirngToWrite)   
    f.close()
