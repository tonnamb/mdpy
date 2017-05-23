import numpy as np
import matplotlib.pyplot as plt

def plot(end_more=0):
    pmf = np.genfromtxt('out.pmf', comments='#')
    hist = np.genfromtxt('out.count', comments='#')

    first_pmf = pmf[0, 1]
    last_pmf = pmf[-1, 1]

    for i, val in enumerate(pmf):
        if first_pmf != val[1]:
            start_i = i
            break

    for i, val in enumerate(pmf[::-1]):
        if last_pmf != val[1]:
            end_i = -i - end_more
            break

    pmf[:, 1] = pmf[:, 1] - max(pmf[(end_i-100):end_i, 1])

    fig = plt.figure()
    plt.plot(pmf[start_i:end_i, 0], pmf[start_i:end_i, 1])
    plt.xlabel(r'z ($\AA$)')
    plt.ylabel('PMF (kcal/mol)')
    fig.savefig('pmf.png', dpi=fig.dpi)
    plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(hist[:,0], hist[:,1])

    plt.xlabel(r'z ($\AA$)')
    plt.ylabel('Histogram Count')
    fig.savefig('histogram.png', dpi=fig.dpi)
    plt.show()
    plt.close()

