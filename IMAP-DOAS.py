from hapi import *
import math as Math
from Constants import *
from DataReader import *
import matplotlib.pyplot as plt
import numpy as np


def gaussianConvolution(rgWl, Fhr, FWHM, wingLength):
    x = np.arange(-wingLength, wingLength + wlStep, wlStep)
    gaussian = SLIT_GAUSSIAN(x, FWHM)
    return rgWl, np.convolve(Fhr, gaussian, mode="same") * wlStep


# Don't forget to clone coefMat
def calcFhr(nu, I0, A, X, Tref):
    Fhr = np.zeros(len(Tref[0]))
    for wl in range(len(Tref[0])):
        exponent = 0
        for iLayer in range(len(X)):
            exponent -= A[iLayer] * X[iLayer] * Tref[iLayer][wl]

        Fhr[wl] = Math.exp(exponent) * 0.1367 * I0[wl]

    return nu, Fhr

def calcFlr():
    A = retrieveA(9, 12)
    X = np.ones((len(atmosLayersAltitude) * len(species), 1))
    rgWl, rgI0, matTref = retrieveLinearlyAlignedI0andTref()
    plot1Matrix(rgWl, matTref)
    rgWl, Fhr = calcFhr(rgWl, rgI0, A, X, matTref)
    rgWl, conv = gaussianConvolution(rgWl, Fhr, FWHM, 3 * FWHM)
    plt.subplot(2, 1, 1)
    plt.plot(rgWl, conv)
    plt.subplot(2, 1, 2)
    plt.plot(rgWl, Fhr)
    plt.show()

initializeConstants((-63, 5))
# verticalProfileCH4, verticalProfilesN2O, verticalProfilesH2O = dataReader.getProfilesAtCoordinate(())
