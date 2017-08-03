from hapi import *
import math as Math
import GlobalConstants as const
from Constants import *
from DataReader import *
import matplotlib.pyplot as plt
import numpy as np

def getValueAtWavelength(wl, rgWl, V):
    index = (wl - rgWl[0]) / const.wlStep
    return V[index]


def gaussianConvolution(Fhr, FWHM, wingLength):
    x = np.arange(-wingLength, wingLength + const.wlStep, const.wlStep)
    gaussian = SLIT_GAUSSIAN(x, FWHM)
    return np.convolve(Fhr, gaussian, mode="same") * const.wlStep


# Don't forget to clone coefMat
def calcFhr(X):
    Fhr = np.zeros(len(const.TrefMat[0]))
    for wl in range(len(const.TrefMat[0])):
        exponent = 0
        for iLayer in range(len(X)):
            exponent -= const.A[iLayer] * X[iLayer] * const.TrefMat[iLayer][wl]

        Fhr[wl] = Math.exp(exponent) * 0.1367 * const.I0[wl]

    return Fhr

def calcFlr(X):
    plot1Matrix(const.wlAxis, const.TrefMat)
    Fhr = calcFhr(X)
    conv = gaussianConvolution(Fhr, const.FWHM, 3 * const.FWHM)
    plt.subplot(2, 1, 1)
    plt.plot(const.wlAxis, conv)
    plt.subplot(2, 1, 2)
    plt.plot(const.wlAxis, Fhr)
    plt.show()


def analyzeImage():
    initializeConstants()

    print 'hi'

analyzeImage()
# initializeConstants(9000, 12, (-63, 5))
# X = const.verticalProfilesH2O.tolist() + const.verticalProfilesN2O.tolist() + const.verticalProfilesCH4.tolist()
# calcFlr(X)

# dataReader.readAvirisData()
#
# wlAxis, I0, TrefMat, X = initializeConstants((-63, 5))
# calcFlr(wlAxis, TrefMat, I0, X)