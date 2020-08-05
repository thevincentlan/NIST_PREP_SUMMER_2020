import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import bruker_io
import os


def element_b(element, start=0, end=None):
    species = os.path.join("Trials", "Element B " + element + ".spx")
    dataset = bruker_io.FittingData(species)
    gauss_fit(dataset, element, start, end)


def glass(sequence, start=0, end=None):
    species = os.path.join("Glass", 'glass_chip_200s_12x12det_' + sequence + '.spx')
    dataset = bruker_io.FittingData(species)
    gauss_fit(dataset, sequence, start, end)


def gauss_fit(dataset, name, start, end):
    bruker_io.bruker_spx_import(dataset)

    x = dataset.energy_scale[start:end]
    y = dataset.channels[start:end]

    # weighted arithmetic mean (corrected - check the section below)
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))

    def Gauss(x, a, x0, sigma):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    popt, pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])

    plt.figure(num=None, figsize=(7, 3), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(x, y, 'k:', label='data')
    plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
    plt.legend()
    plt.title(name + ' XRF')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.show()
