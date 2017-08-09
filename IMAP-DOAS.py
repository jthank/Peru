from hapi import *
import math as Math
import GlobalConstants as const
from Constants import *
from DataReader import *
import matplotlib.pyplot as plt
import numpy as np


def gaussianConvolution(Fhr, FWHM, wingLength):
    x = np.arange(-wingLength, wingLength + const.wlStep, const.wlStep)
    gaussian = SLIT_GAUSSIAN(x, FWHM)
    return np.convolve(Fhr, gaussian, mode="same") * const.wlStep


def calcFhr(A, T, X):
    Fhr = np.zeros(len(T[0]))
    for wl in range(len(T[0])):
        exponent = 0
        for iLayer in range(len(X)):
            exponent -= A[iLayer] * X[iLayer] * T[iLayer][wl]

        Fhr[wl] = Math.exp(exponent) * 0.1367 * const.I0[wl]

    return Fhr

def calcFlr(A, T, X):
    plot1Matrix(const.wlAxis, T)
    Fhr = calcFhr(A, T, X)

    conv = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv2 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv3 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv4 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv5 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv6 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # conv7 = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)

    # plt.subplot(2, 1, 1)
    plt.plot(const.wlAxis, conv, label='1')
    # plt.plot(const.wlAxis, conv2, label='2')
    # plt.plot(const.wlAxis, conv3, label='3')
    # plt.plot(const.wlAxis, conv4, label='4')
    # plt.plot(const.wlAxis, conv5, label='5')
    # plt.plot(const.wlAxis, conv6, label='6')
    # plt.plot(const.wlAxis, conv7, label='7')
    # conv2 = gaussianConvolution(np.multiply(const.I0, 0.1367), const.FWHM, 3 * const.FWHM)
    # plt.plot(const.wlAxis, conv2)
    #
    # for i in range(75):
    #     radiances = []
    #     for iC in range(len(const.AVIRIS_CenterIndexes)):
    #         radiances.append(0.04 * Math.pi * const.arrIMG[const.pixelX+i][const.pixelY+i][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC] / const.arrRFL[const.pixelX+i][const.pixelY+i][const.AVIRIS_CenterIndexes[iC]])
    #     plt.plot(const.bands, np.multiply(2,radiances))



    # plt.subplot(2, 1, 2)
    # plt.plot(const.wlAxis, Fhr)
    # plt.show()
    return conv


def isValidPixel():
    return True

def retrieveAMF(aircraftAlt, SZA, CZA, lowAlt, highAlt):
    cameraZenithFactor = 1 if aircraftAlt > highAlt else 0 if aircraftAlt < lowAlt else (aircraftAlt - lowAlt) / (highAlt - lowAlt)
    # replace 1 with (highAlt - lowAlt)
    return 1.0 / Math.cos(Math.radians(SZA)) + \
           cameraZenithFactor * 1.0 / Math.cos(Math.radians(CZA))

def retrieveHeightAdjustedATX(iWidth, iLength):
    pixelElevation = const.arrIGM[iWidth][iLength][2] #check this
    T = [row[:] for row in const.TrefMat]
    X = const.verticalProfiles[:]
    iAlt = len(const.levelAltitudes) - 1
    while iAlt > 0 and const.levelAltitudes[iAlt] < pixelElevation + 75:
        length = len(T)
        for iSpec in range(len(const.species), 0, -1):
            del T[length / len(const.species) * iSpec - 1]
            del X[length / len(const.species) * iSpec - 1]

        iAlt -= 1

    alts = const.levelAltitudes[:iAlt + 1]
    aircraftAlt = const.arrORT[iWidth][iLength][0] * Math.cos(Math.radians(const.arrORT[iWidth][iLength][2]))
    A = np.empty(len(X))
    CZA = const.arrORT[iWidth][iLength][2]
    SZA = const.arrORT[iWidth][iLength][4]

    for iAlt in range(len(alts)):
        low = 0.5 * alts[iAlt + 1] + 0.5 * alts[iAlt] if iAlt < len(alts) - 1 else pixelElevation
        high = 0.5 * alts[iAlt - 1] + 0.5 * alts[iAlt] if iAlt > 0 else 100000
        # times 100?
        AMF = retrieveAMF(aircraftAlt, SZA, CZA, low, high)
        for iSpec in range(0, len(const.species)):
            A[len(A) / len(const.species) * iSpec + iAlt] = AMF

    return A, T, X


def getError(radiances, Flr, currErrors):
    sum = 0
    for iBand in range(len(const.bands)):
        sum += abs(Flr[const.wlAxisBandIndexes[iBand]] - radiances[iBand])

    print sum
    currErrors.append(sum)

def applyGradient(A, T, X, radiances, alpha, currErrors):
    Flr = calcFlr(A, T, X)
    getError(radiances, Flr, currErrors)
    for iX in range(len(X)):
        for iBand in range(len(const.bands)):
            a = 2 * radiances[iBand] * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
            b = 2 * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
            c = Math.sqrt(radiances[iBand] * radiances[iBand] - 2 * radiances[iBand] * Flr[const.wlAxisBandIndexes[iBand]] + Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]])
            X[iX] -= alpha * (a - b) / c


def gradientDescent(A, T, X, radiances):
    base = Math.pow(10, 33)
    Xs = []
    alphas = [1 * base, 2 * base, 3 * base, 4 * base]
    allErrors = []
    for i in range(4):
        Xs.append(X[:])

    for iX in range(len(Xs)):
        currErrors = []
        for i in range(100):
            applyGradient(A, T, Xs[iX], radiances, alphas[iX], currErrors)

        allErrors.append(currErrors)

    for i in range(4):
        plt.subplot(2, 2, i + 1)
        plt.plot(np.arange(0, 100), allErrors[i])

    plt.show()

    for i in range(4):
        plt.subplot(2, 2, i + 1)
        plt.plot(np.arange(0, len(X)), np.divide(Xs[i], X))

    plt.show()

    for i in range(4):
        plt.subplot(2, 2, i + 1)
        conv = calcFlr(A, T, Xs[i])
        plt.plot(const.wlAxis, conv)
        plt.plot(const.bands, radiances)

    plt.show()

def adjustMeasuredRadiances():
    radiances = []
    for iC in range(len(const.AVIRIS_CenterIndexes)):
        radiances.append(0.04 * Math.pi * const.arrIMG[const.pixelX][const.pixelY][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC] / const.arrRFL[const.pixelX][const.pixelY][const.AVIRIS_CenterIndexes[iC]])

    return radiances

def adjustMeasuredRadiances2():
    radiances = []
    for iC in range(len(const.AVIRIS_CenterIndexes)):
        radiances.append(0.04 * Math.pi * const.arrIMG[const.pixelX][const.pixelY][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC])

    return radiances

def analyzeImage():
    initializeConstants()
    # for iWidth in range(const.arrIMG.shape[0]):
    #     for iLength in range(const.arrIMG.shape[1]):
    #         if not isValidPixel(iWidth, iLength):
    #             continue

    A, T, X = retrieveHeightAdjustedATX(const.pixelX, const.pixelY)
    calcFlr(A, T, X)

    radiances = adjustMeasuredRadiances()
    plt.plot(const.bands, radiances)
    radiances2 = adjustMeasuredRadiances2()
    plt.plot(const.bands, radiances2)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)
    plt.show()
    # gradientDescent(A, T, X, radiances)
    # calcFlr(A, T, X)



    print 'hi'

# dataReader.readAvirisData()
# name = 'blah'
# fetch(name, 1, 1, const.wnLow, const.wnHigh)
# nu, coef = absorptionCoefficient_Lorentz(SourceTables=name, HITRAN_units=True)
# plt.plot(nu, coef)
# plt.show()
analyzeImage()
# initializeConstants(9000, 12, (-63, 5))
# X = const.verticalProfilesH2O.tolist() + const.verticalProfilesN2O.tolist() + const.verticalProfilesCH4.tolist()
# calcFlr(X)

# dataReader.readAvirisData()
#
# wlAxis, I0, TrefMat, X = initializeConstants((-63, 5))
# calcFlr(wlAxis, TrefMat, I0, X)