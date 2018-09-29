import pylab
import collections

from ps7_v1 import *
from ps8S1 import *


def problem6():
    """
    Runs simulations and make histograms for problem 6.
    Runs multiple simulations to show the relationship between administration
    of multiple drugs and patient outcome.
    Histograms of final total virus populations are displayed for lag times of
    150, 75, 0 timesteps between adding drugs (followed by an additional 150
    timesteps of simulation).
    """
    viruses = [ResistantVirus(0.1, 0.05, {'guttagonol':False, 'grimpex':False}, 0.005) for _ in range(100)]
    num_trials = 30
    print '\n%d patients -- 150 time steps, add guttagonol, 300, 150, 75 and 0 time steps, add grimpex, 150 more time steps\n' % num_trials

    delays = collections.defaultdict(list) # a dict of list, which contains the final virus population of each trial(patient)

    for steps, delay in [(600,300), (450,150), (375,75), (300,0)]:
        print '%d patients, 150 time steps, add guttagonol, %d time steps, add grimpex, 150 more time steps\n' % (num_trials, delay)
        for n in range(num_trials):
            patient = Patient(viruses, 1000)
            for i in range(steps):
                total = patient.update()
                if i == 150:
                    patient.addPrescription('guttagonol')
                if i  == delay + 150:
                    patient.addPrescription('grimpex')
            delays[delay].append(total)

    cuered_rates = {}
    for k, v in delays.items():
        cuered = [p for p in v if p <= 50]
        cuered_rates[k] = str(100 * len(cuered) / float(len(v))) + '%'

    for delay in [300, 150, 75, 0]:
        pylab.figure()
        pylab.hist(delays[delay])
        pylab.title('%d patients, 150 time steps, add guttagonol, %d time steps, add grimpex, more 150 time steps' % (delay, delay))
        pylab.xlabel('Final total virus population -- %s cuered' % cuered_rates[delay])
        pylab.ylabel('Number of patients')

    pylab.show()

problem6()

#
# PROBLEM 7
#

def problem7():
    """
    Run simulations and plot graphs examining the relationship between
    administration of multiple drugs and patient outcome.
    Plots of total and drug-resistant viruses vs. time are made for a
    simulation with a 300 time step delay between administering the 2 drugs and
    a simulations for which drugs are administered simultaneously.
    """
    viruses = [ResistantVirus(0.1, 0.05, {'guttagonol':False, 'grimpex':False}, 0.005) for _ in range(100)]

    print '150 time steps, add guttagonol, 300 time steps, add grimpex, 150 more time steps\n'
    patient = Patient(viruses, 1000)
    total_pop = []
    resist_guttagonol_pop = []
    resist_grimpex_pop = []
    resist_both_pop = []
    for i in range(600):
        total_pop.append(patient.update())
        resist_guttagonol_pop.append(patient.getResistPop(['guttagonol']))
        resist_grimpex_pop.append(patient.getResistPop(['grimpex']))
        resist_both_pop.append(patient.getResistPop(['guttagonol', 'grimpex']))
        if i == 150:
            patient.addPrescription('guttagonol')
        if i  == 450:
            patient.addPrescription('grimpex')

    pylab.figure()
    pylab.plot(range(1,601), total_pop, color='red', label='total virus population')
    pylab.plot(range(1,601), resist_guttagonol_pop, color='blue', label='guttagonol-resistant')
    pylab.plot(range(1,601), resist_grimpex_pop, color='yellow', label='grimpex-resistant')
    pylab.plot(range(1,601), resist_both_pop, color='green', label='both resistant')
    #pylab.xticks(range(0,601,100))
    pylab.xlabel('Timesteps')
    pylab.ylabel('Virus population')
    pylab.title('150 time steps, add guttagonol, 300 time steps, add grimpex, 150 more time steps')
    pylab.legend(loc='upper left')


    print '150 time steps, add guttagonol and grimpex simultaneoulsy, 150 more time steps'
    patient = Patient(viruses, 1000)
    total_pop = []
    resist_guttagonol_pop = []
    resist_grimpex_pop = []
    resist_both_pop = []
    for i in range(300):
        total_pop.append(patient.update())
        resist_guttagonol_pop.append(patient.getResistPop(['guttagonol']))
        resist_grimpex_pop.append(patient.getResistPop(['grimpex']))
        resist_both_pop.append(patient.getResistPop(['guttagonol', 'grimpex']))
        if i == 150:
            patient.addPrescription('guttagonol')
            patient.addPrescription('grimpex')

    pylab.figure()
    pylab.plot(range(1,301), total_pop, color='red', label='total virus population')
    pylab.plot(range(1,301), resist_guttagonol_pop, color='blue', label='guttagonol-resistant')
    pylab.plot(range(1,301), resist_grimpex_pop, color='yellow', label='grimpex-resistant')
    pylab.plot(range(1,301), resist_both_pop, color='green', label='both resistant')
    #pylab.xticks(range(0,601,100))
    pylab.xlabel('Timesteps')
    pylab.ylabel('Virus population')
    pylab.title('150 time steps, add guttagonol and grimpex simultaneoulsy, 150 more time steps')
    pylab.legend(loc='upper right')

    pylab.show() # show the figure


#problem7()
