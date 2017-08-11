from hapi import *
import math as Math
import GlobalConstants as const
from Constants import *
from DataReader import *
from multiprocessing.dummy import Pool as ThreadPool
import matplotlib.pyplot as plt
import numpy as np


def gaussianConvolution(Fhr, FWHM, wingLength):
    x = np.arange(-wingLength, wingLength + const.wlStep, const.wlStep)
    gaussian = SLIT_GAUSSIAN(x, FWHM)
    return np.convolve(Fhr, gaussian, mode="same") * const.wlStep



def calcFhr(A, T, X, pixel):
    Fhr = np.zeros(len(T[0]))
    for iWl in range(len(T[0])):
        exponent = 0
        for iLayer in range(len(X)):
            exponent -= A[iLayer] * X[iLayer] * T[iLayer][iWl]

        Fhr[iWl] = Math.exp(exponent) * const.I0[iWl] * Math.cos(Math.radians(const.arrORT[pixel[0]][pixel[1]][4]))

    return Fhr


def calcFlr(A, T, X, pixel):
    plot1Matrix(const.wlAxis, T)
    Fhr = calcFhr(A, T, X, pixel)
    conv = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    plt.plot(const.wlAxis, conv, label='Modeled Irradiance')
    return conv


def isValidPixel(pixelX, pixelY):
    if const.arrIMG[pixelX][pixelY][0] < 0:
        return False

    return True


def retrieveAMF(aircraftAlt, SZA, CZA, lowAlt, highAlt):
    cameraZenithFactor = 1 if aircraftAlt > highAlt else 0 if aircraftAlt < lowAlt else (aircraftAlt - lowAlt) / (highAlt - lowAlt)
    # replace 1 with (highAlt - lowAlt)
    return 1.0 / Math.cos(Math.radians(SZA)) + cameraZenithFactor * 1.0 / Math.cos(Math.radians(CZA))


def retrieveHeightAdjustedATX(pixelX, pixelY):
    pixelElevation = const.arrIGM[pixelX][pixelY][2]
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
    aircraftAlt = const.arrORT[pixelX][pixelY][0] * Math.cos(Math.radians(const.arrORT[pixelX][pixelY][2]))
    A = np.empty(len(X))
    CZA = const.arrORT[pixelX][pixelY][2]
    SZA = const.arrORT[pixelX][pixelY][4]

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


def applyGradient(A, T, X, radiances, alpha, currErrors, pixel):
    Flr = calcFlr(A, T, X, pixel)
    getError(radiances, Flr, currErrors)
    for iX in range(len(X)):
        for iBand in range(len(const.bands)):
            a = 2 * radiances[iBand] * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
            b = 2 * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
            c = Math.sqrt(radiances[iBand] * radiances[iBand] - 2 * radiances[iBand] * Flr[const.wlAxisBandIndexes[iBand]] + Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]])
            X[iX] -= alpha * (a - b) / c


def gradientDescent(A, T, X, radiances, pixel):
    base = Math.pow(10, 33)
    Xs = []
    alphas = [1 * base, 2 * base, 3 * base, 4 * base]
    allErrors = []
    for i in range(4):
        Xs.append(X[:])

    for iX in range(len(Xs)):
        currErrors = []
        for i in range(100):
            applyGradient(A, T, Xs[iX], radiances, alphas[iX], currErrors, pixel)

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
        conv = calcFlr(A, T, Xs[i], pixel)
        plt.plot(const.wlAxis, conv)
        plt.plot(const.bands, radiances)

    plt.show()


def adjustMeasuredRadiances3():
    radiances = []
    for iC in range(len(const.AVIRIS_CenterIndexes)):
        radiances.append(0.04 * Math.pi * const.arrIMG[const.pixelX][const.pixelY][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC] / const.arrRFL[const.pixelX][const.pixelY][const.AVIRIS_CenterIndexes[iC]])

    return radiances


def adjustMeasuredRadiances(pixelX, pixelY):
    radiances = []
    for iC in range(len(const.AVIRIS_CenterIndexes)):
        radiances.append(
            (0.01 * Math.pi * const.arrIMG[pixelX][pixelY][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC])
            / (const.arrRFL[pixelX][pixelY][const.AVIRIS_CenterIndexes[iC]] * Math.cos(Math.radians(const.arrORT[pixelX][pixelY][4]))))

    return radiances


xStartStop = (400, 420)
yStartStop = (400, 420)


def processPixel(pixels):
    print "Processing " + str(pixels[0]) + ", " + str(pixels[1])
    if not isValidPixel(pixels[0], pixels[1]):
        return 0

    A, T, X = retrieveHeightAdjustedATX(pixels[0], pixels[1])
    initialMethaneConcentraion = X[len(X) - 1]
    radiances = adjustMeasuredRadiances(pixels[0], pixels[1])
    gradientDescent(A, T, X, radiances, pixels)
    return X[len(X) - 1] / initialMethaneConcentraion


def analyzeImageMultiThreaded():
    initializeConstants()
    pixelTuples = []
    for pixelX in range(xStartStop[0], xStartStop[1]):
        for pixelY in range(yStartStop[0], yStartStop[1]):
            pixelTuples.append((pixelX, pixelY))


    pool = ThreadPool(8)
    print pool.map(processPixel, pixelTuples)

def analyzeImageSingleThread():
    pixel = (xStartStop[0], yStartStop[0])
    initializeConstants()
    A, T, X = retrieveHeightAdjustedATX(pixel[0], pixel[1])
    Flr = calcFlr(A, T, X, pixel)
    plt.plot(const.wlAxis, const.I0)
    print const.FWHM
    conv = gaussianConvolution(const.I0, const.FWHM, 1.5 * const.FWHM)
    plt.plot(const.wlAxis, conv)
    for pixelX in range(xStartStop[0], xStartStop[1]):
        for pixelY in range(yStartStop[0], yStartStop[1]):
            radiances = adjustMeasuredRadiances(pixelX, pixelY)
            plt.plot(const.bands, radiances)
    plt.show()

    # processPixel(pixel)


    # A, T, X = retrieveHeightAdjustedATX(const.pixelX, const.pixelY)
    # calcFlr(A, T, X)
    #
    # radiances = adjustMeasuredRadiances3()
    # radiances2 = adjustMeasuredRadiances(const.pixelX, const.pixelY)
    # plt.plot(const.bands, radiances2, label='Gathered Radiances')
    # # plt.plot(const.bands, radiances, label='Gathered Irradiance (adjusted)')
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=2, mode="expand", borderaxespad=0.)
    # plt.show()
    # gradientDescent(A, T, X, radiances2)
    # calcFlr(A, T, X)
    # plt.plot(const.bands, radiances2, label='Gathered Radiances')
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=2, mode="expand", borderaxespad=0.)
    # plt.show()

analyzeImageSingleThread()
