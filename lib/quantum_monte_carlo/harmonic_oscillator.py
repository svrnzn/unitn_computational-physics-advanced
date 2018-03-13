import numpy
import math

np_random = numpy.random.RandomState(0)# seed set to 0 for DEBUGGING purposes

if __name__ == '__main__':
    # SETTINGS
    domain_width = 10# the domain is taken centered in 0
    N_w = 10000# number of initial walkers

    #
    delta_tau = 0.01# better naming?
    tau = 2

    # guess of the ground state energy
    trial_energy = 0.5

    '''
    # psi_T settings
    def trial_wave_function(x):
        if x >= - domain_width / 2 and x <= domain_width / 2:
            return 1 / domain_width
        else:
            return 0
            '''

    # output settings
    N_bins = 101 # MUST BE AN ODD NUMBER
    bins_width = domain_width / N_bins
    histogram = numpy.zeros(N_bins)
    
    walkers = np_random.uniform(- 0.5 * domain_width, 0.5 * domain_width, N_w)
    for walker in walkers:
        histogram[round(walker / bins_width) + int(N_bins / 2)] += 1
    histogram = histogram / (N_w * bins_width)

    output = open('tau0.dat', 'w')
    for i in range(N_bins):
        output.write(str((i - int(N_bins / 2)) * bins_width) + '\t' + str(histogram[i]) + '\n')
    output.close()

    for i in range(int(tau / delta_tau) + 1):# + ?
        walkers = walkers + np_random.normal(0, math.sqrt(delta_tau), N_w)
        walkers_multiplicity = numpy.empty(N_w)
        for j in range(N_w):
            walkers_multiplicity[j] = int(math.exp(- (0.5 * walkers[j]**2 - trial_energy) * delta_tau) + np_random.random_sample())

#        energy_numerator += numpy.sum()
#        energy_denominator += numpy.sum()

        new_walkers = numpy.empty(numpy.sum(walkers_multiplicity))
        '''
        k = 0
        for j in range(N_w):
            new_walkers[k:k + walkers_multiplicity[j]] = walkers[j] * numpy.ones(walkers_multiplicity[j])
            k += walkers_multiplicity[j]
            '''
        l = 0
        for j in range(N_w):
            for k in range(int(walkers_multiplicity[j])):
                new_walkers[l] = walkers[j]
                l += 1

        walkers = new_walkers
        N_w = len(walkers)
        print(N_w)

    histogram = numpy.zeros(N_w)
    for j in range(N_w):
        if walkers[j] > - domain_width / 2 and walkers[j] < domain_width / 2:
            histogram[round(walkers[j] / bins_width) + int(N_bins / 2)] += 1
        else:
            print('ocio!')
    histogram = histogram / (N_w * bins_width)
#    print(numpy.sum(histogram) * N_w * bins_width)
    output = open('tau' + str(tau) + '.dat', 'w')
    for i in range(N_bins):
        output.write(str((i - int(N_bins / 2)) * bins_width) + '\t' + str(histogram[i]) + '\n')
    output.close()
