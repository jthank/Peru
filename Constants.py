import numpy as np
from hapi import *
import math as Math
import matplotlib.pyplot as plt

FWHM = 10
wavelengthLow = 2278
wavelengthHigh = 2358
bands = 9
species = [['H20', 1, 1], ['N20', 4, 1], ['CH4', 6, 1]]
# Boundary is everything above the number ex: layer 2 (index 1) is 10->20
atmosLayersAltitude = [20, 10, 7, 6, 5, 4, 3, 2, 1, 0]
atmosLayerT = [300, 299, 298, 297, 296, 295, 294, 293, 292, 291]
atmosLayerP = [1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]


def retrieveA(altitude, SZA):
    SZAradians = np.deg2rad(SZA)
    arrA = np.zeros((len(atmosLayersAltitude) * len(species), 1))

    for iLayer in range(len(atmosLayersAltitude)):
        for iSpec in range(len(species)):
            # if the plane is entirely above the layer, it has the additional 1.
            arrA[iLayer + iSpec * len(atmosLayersAltitude)][0] = \
                1 / np.cos(SZAradians) if altitude < atmosLayersAltitude[iLayer] else 1 + 1 / np.cos(SZAradians)

    return np.matrix(arrA)


def retrieveI0():
    return 10


# def runIteration():
#     wavelengths = np.linspace(wavelengthLow, wavelengthHigh, bands)
#     TrefPerWavelength = np.zeros((len(atmosLayersAltitude) * len(species), bands))
#
#     # Take each value and multiply by A and X. Then take e ^ negative value. Multiply by I_0. That is F_hr
#
#     # Sample at each wavelength
#     for iWl in range(len(wavelengths)):
#         wn = swapWnWl(wavelengths[iWl])
#         index = int(round((wn - wnLow) / WnStep))
#         TrefPerWavelength[len(atmosLayersAltitude) * iSpec + iLayer][iWl] = 0 if index >= len(coef) else (
#         coef[index] * Math.log(10, Math.e))
#
#     # Graph
#     for j in range(len(coef)):
#         nu[j] = swapWnWl(nu[j])
#         # absorbance = optical density / ln(10)
#         coef[j] *= Math.log(10, Math.e)
#     plt.subplot(5, 2, iLayer + 1)
#     plt.plot(nu, coef)
#     plt.show()

def fitTrefToMatrix(TrefNusList, TrefCoefsList, wnStep):
    # First find the min and max over all Nus
    min = TrefNusList[0][0]
    max = TrefNusList[0][len(TrefNusList[0]) - 1]
    for nuList in TrefNusList:
        if(nuList[0] < min):
            min = nuList[0]
        if(nuList[len(nuList) - 1] > max):
            max = nuList[len(nuList) - 1]

    # Create the final matrix with correct dimensions (1.01 to account for double -> int precision errors)
    nuMat = np.empty((len(TrefCoefsList), int((max - min) / wnStep + 1.01)))
    coefMat = np.empty((len(TrefCoefsList), int((max - min) / wnStep + 1.01)))


    for iList in range(len(TrefNusList)):
        nuList = TrefNusList[iList]
        nu = np.asarray(nuList)
        coef = np.asarray(TrefCoefsList[iList])

        # Add linearly to nu and zeros to coef in order to equalize the dimensions for all lines
        if (nuList[0] > min):
            prependNu = np.arange(min, nuList[0] - wnStep / 2.0, wnStep)
            prependCoef = np.zeros(len(prependNu))
            nu = np.concatenate((prependNu, nu))
            coef = np.concatenate((prependCoef, coef))

        # Adds to the end rather than the beginning
        if (nuList[len(nuList) - 1] < max):
            appendNu = np.arange(nuList[len(nuList) - 1] + wnStep, max + wnStep / 2.0, wnStep)
            appendCoef = np.zeros(len(appendNu))
            nu = np.concatenate((nu, appendNu))
            coef = np.concatenate((coef, appendCoef))

        nuMat[iList] = nu
        coefMat[iList] = coef

    return nuMat, coefMat



def retrieveTref():
    wnLow = swapWnWl(wavelengthHigh + 2 * FWHM)
    wnHigh = swapWnWl(wavelengthLow - 2 * FWHM)
    wnStep = 0.01
    TrefNusList = []
    TrefCoefsList = []
    TrefNus = np.zeros((len(atmosLayersAltitude) * len(species), int((wnHigh - wnLow) / wnStep + 1)))
    TrefCoefs = np.zeros((len(atmosLayersAltitude) * len(species), int((wnHigh - wnLow) / wnStep + 1)))

    for iSpec in range(len(species)):
        fetch(species[iSpec][0], species[iSpec][1], species[iSpec][2], wnLow, wnHigh)

        for iLayer in range(len(atmosLayersAltitude)):
            nu, coef = absorptionCoefficient_Lorentz(SourceTables=species[iSpec][0], HITRAN_units=False,
                                                     Environment={'T': atmosLayerT[iLayer], 'p': atmosLayerP[iLayer]},
                                                     WavenumberStep=wnStep)
            nu, coef = absorptionSpectrum(nu, coef)

            for iNu in range(len(nu)):
                nu[iNu] = round(nu[iNu], 2)

            TrefNusList.append(nu)
            TrefCoefsList.append(coef)

    return fitTrefToMatrix(TrefNusList, TrefCoefsList, wnStep)

def swapWnWl(WlWn):
    return (1.0 / WlWn) * 10000000