from hapi import *
import math as Math
import GlobalConstants as const
from Constants import *
from DataReader import *
from multiprocessing import Pool as ThreadPool
from PIL import Image
from matplotlib import colors
import matplotlib.pyplot as plt
import Queue
import numpy as np
import time




def plotBars(xAxis, bars):
    barWidth = 0.9 / len(bars)
    indexes = np.arange(len(bars[0]))
    for iBars in range(len(bars)):
        plt.bar(indexes + iBars * barWidth, bars[iBars], barWidth)

    plt.xticks(indexes, xAxis, rotation='vertical')

def gaussianConvolution(Fhr):
    return np.convolve(Fhr, const.gaussian, mode="same") * const.wlStep


def calcFhr(A, T, X, pixel):
    Fhr = np.zeros(len(T[0]))
    for iWl in range(len(T[0])):
        exponent = 0
        for iLayer in range(len(X)):
            exponent -= A[iLayer] * X[iLayer] * T[iLayer][iWl]

        Fhr[iWl] = Math.exp(exponent) * const.I0[iWl] * Math.cos(Math.radians(const.arrORT[pixel[0]][pixel[1]][4]))

    return Fhr


def calcFhr1(A1, A2, T, X, X0, pixel):
    Fhr = np.zeros(len(T[0]))
    X0H = np.zeros(len(X0))
    X0H[:len(X0) / 3] = X0[:len(X0) / 3]
    Xdiff = np.subtract(X, X0H)
    for iWl in range(len(T[0])):
        exponent = 0
        exponent2 = 0
        for iLayer in range(len(X)):
            exponent -= A1[iLayer] * X[iLayer] * T[iLayer][iWl]
            exponent2 -= A2[iLayer] * Xdiff[iLayer] * T[iLayer][iWl]

        Fhr[iWl] = Math.exp(exponent + exponent2) * const.I0[iWl] * Math.cos(Math.radians(const.arrORT[pixel[0]][pixel[1]][4]))

    return Fhr


def calcFhr2(A1, A2, T, X, X0, pixel):
    Fhr = np.zeros(len(T[0]))
    Xdiff = np.subtract(X, X0)

    for iWl in range(len(T[0])):
        exponent = 0
        exponent2 = 0
        for iLayer in range(len(X)):
            exponent -= A1[iLayer] * X[iLayer] * T[iLayer][iWl]
            exponent2 -= A2[iLayer] * Xdiff[iLayer] * T[iLayer][iWl]

        Fhr[iWl] = Math.exp(exponent + exponent2) * const.I0[iWl] * Math.cos(Math.radians(const.arrORT[pixel[0]][pixel[1]][4]))

    return Fhr


def calcFlr(A, T, X, pixel):
    plot1Matrix(const.wlAxis, T)
    Fhr = calcFhr(A, T, X, pixel)
    conv = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # plt.plot(const.wlAxis, conv, label='Modeled Irradiance')
    return conv

def calcFlr2(A1, A2, T, X, X0, pixel):
    plot1Matrix(const.wlAxis, T)
    Fhr = calcFhr2(A1, A2, T, X, X0, pixel)
    conv = gaussianConvolution(Fhr, const.FWHM, 1.5 * const.FWHM)
    # plt.plot(const.wlAxis, conv, label='Modeled Irradiance')
    return conv


def isValidPixel(pixelX, pixelY):
    if const.arrIMG[pixelX][pixelY][0] < 0:
        return False

    return True


def retrieveAMF(aircraftAlt, SZA, CZA, lowAlt, highAlt):
    cameraZenithFactor = 1 if aircraftAlt > highAlt else 0 if aircraftAlt < lowAlt else (aircraftAlt - lowAlt) / (highAlt - lowAlt)
    # replace 1 with (highAlt - lowAlt)
    return 1.0 / Math.cos(Math.radians(SZA)) + cameraZenithFactor * (1.0 / Math.cos(Math.radians(CZA)))

def retrieveAMF2(aircraftAlt, SZA, CZA, lowAlt, highAlt):
    cameraZenithFactor = 1 if aircraftAlt > highAlt else 0 if aircraftAlt < lowAlt else (aircraftAlt - lowAlt) / (highAlt - lowAlt)
    # replace 1 with (highAlt - lowAlt)
    return 1.0 / Math.cos(Math.radians(SZA)), cameraZenithFactor * (1.0 / Math.cos(Math.radians(CZA)))


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

def retrieveHeightAdjustedATX2(pixelX, pixelY):
    # pixelElevation = const.arrIGM[pixelX][pixelY][2]
    pixelElevation = const.arrIGM[300][300][2]
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
    A1 = np.empty(len(X))
    A2 = np.empty(len(X))
    CZA = const.arrORT[pixelX][pixelY][2]
    SZA = const.arrORT[pixelX][pixelY][4]

    for iAlt in range(len(alts)):
        low = 0.5 * alts[iAlt + 1] + 0.5 * alts[iAlt] if iAlt < len(alts) - 1 else pixelElevation
        high = 0.5 * alts[iAlt - 1] + 0.5 * alts[iAlt] if iAlt > 0 else 100000
        # times 100?
        AMF1, AMF2 = retrieveAMF2(aircraftAlt, SZA, CZA, low, high)
        for iSpec in range(0, len(const.species)):
            A1[len(A1) / len(const.species) * iSpec + iAlt] = AMF1
            A2[len(A2) / len(const.species) * iSpec + iAlt] = AMF2


    return A1, A2, T, X


def getError(radiances, Flr):
    sum = 0
    for iBand in range(len(const.bands)):
        sum += abs(Flr[const.wlAxisBandIndexes[iBand]] - radiances[iBand])

    print sum

    return sum


def applyGradient(A, T, X, radiances, alpha, currErrors, pixel, minX, minError):
    Flr = calcFlr(A, T, X, pixel)
    error = getError(radiances, Flr)
    currErrors.append(error)

    if error < minError:
        minError = error
        minX = X[:]

    layers = len(X) / len(const.species)
    for iSpecBoundaries in range(layers, len(const.species) * layers + 1, layers):
        for iX in range(iSpecBoundaries - 1, iSpecBoundaries - const.changeableLayers - 1, -1):
            for iBand in range(len(const.bands)):
                a = 2 * radiances[iBand] * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
                b = 2 * A[iX] * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
                c = Math.sqrt(radiances[iBand] * radiances[iBand] - 2 * radiances[iBand] * Flr[const.wlAxisBandIndexes[iBand]] + Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]])
                X[iX] -= alpha * (a - b) / c

    return minX, minError

def applyGradient2(A1, A2, T, X, X0, radiances, currErrors, pixel, minX, minError):
    Flr = calcFlr2(A1, A2, T, X, X0, pixel)
    error = getError(radiances, Flr)
    currErrors.append(error)

    if error < minError:
        minError = error
        minX = X[:]

    layers = len(X) / len(const.species)
    for iSpecBoundaries in range(layers, len(const.species) * layers + 1, layers):
        for iX in range(iSpecBoundaries - 1, iSpecBoundaries - const.changeableLayers - 1, -1):
            # springTopBound = 2.5
            # currentFactor = X[iX] / X0[iX]
            # springFactor = currentFactor if currentFactor < 1 else Math.pow((1 + 1.0 / (springTopBound - 1)) - currentFactor / (springTopBound - 1), 2)
            springFactor = 1

            for iBand in range(len(const.bands)):
                a = 2 * radiances[iBand] * (A1[iX] + A2[iX]) * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
                b = 2 * (A1[iX] + A2[iX]) * T[iX][const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]]
                c = Math.sqrt(radiances[iBand] * radiances[iBand] - 2 * radiances[iBand] * Flr[const.wlAxisBandIndexes[iBand]] + Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]])
                X[iX] -= springFactor * const.alpha * (a - b) / c

    return minX, minError


def gradientDescent(A, T, X, radiances, pixel):
    alpha = 4 * Math.pow(10, 32)
    errors = []
    minX = []
    minError = sys.float_info.max
    for i in range(100):
        minX, minError = applyGradient(A, T, X, radiances, alpha, errors, pixel, minX, minError)

    return minX

def gradientDescent2(A1, A2, T, X, X0, radiances, pixel):
    errors = []
    minX = []
    minError = sys.float_info.max
    q = Queue.Queue()
    for i in range(5):
        q.put(sys.float_info.max)

    for i in range(const.iterations):
        minX, minError = applyGradient2(A1, A2, T, X, X0, radiances, errors, pixel, minX, minError)
        if q.get() - minError < const.stopThreshold:
            return minX, minError, errors

        q.put(minError)

    return minX, minError, errors


def gradientDescent(A1, A2, T, X, X0, radiances, pixel):
    errors = []
    minX = []
    minError = sys.float_info.max
    q = Queue.Queue()
    for i in range(4):
        q.put(sys.float_info.max)

    for i in range(const.iterations):
        minX, minError = applyGradient2(A1, A2, T, X, X0, radiances, errors, pixel, minX, minError)

    return minX


def gradientDescentMultiAlphas(A, T, X, radiances, pixel):
    base = Math.pow(10, 32)
    Xs = []
    alphas = [base, 5 * base, 40 * base, 100 * base]
    #100 * base, 500 * base, 1000 * base, 5000 * base
    allErrors = []
    minXs = []
    for i in range(len(alphas)):
        Xs.append(X[:])

    for iX in range(len(Xs)):
        currErrors = []
        minX = []
        minError = sys.float_info.max
        for i in range(100):
            minX, minError = applyGradient(A, T, Xs[iX], radiances, alphas[iX], currErrors, pixel, minX, minError)

        minXs.append(minX)
        allErrors.append(currErrors)

    for i in range(len(Xs)):
        plt.subplot(2, 2, i + 1)
        plt.plot(np.arange(0, 100), allErrors[i])

    plt.show()

    for i in range(len(Xs)):
        plt.subplot(2, 2, i + 1)

        layers = len(X) / len(const.species)

        xAxis = []
        for iSpec in range(len(const.species)):
            xAxis.append(const.species[iSpec][0])


        scalingFactors = []
        for iScale in range(const.changeableLayers, 0, -1):
            scalingFactorsAtIndex = []
            for iSpecBoundaries in range(layers, len(const.species) * layers + 1, layers):
                scalingFactorsAtIndex.append(minXs[i][iSpecBoundaries - iScale] / X[iSpecBoundaries - iScale])

            scalingFactors.append(scalingFactorsAtIndex)

        plotBars(xAxis, scalingFactors)

    plt.show()



    for i in range(len(Xs)):
        plt.subplot(2, 2, i + 1)
        conv = calcFlr(A, T, minXs[i], pixel)
        plt.plot(const.wlAxis, conv)
        plt.plot(const.bands, radiances)

    plt.show()


def gradientDescentMultiAlphas2(A1, A2, T, X, X0, radiances, pixel):
    iterations = 100
    base = Math.pow(10, 33)
    Xs = []
    alphas = [const.alpha]
    allErrors = []
    minXs = []
    for i in range(len(alphas)):
        Xs.append(X[:])

    for iX in range(len(Xs)):
        minX, minError, currErrors = gradientDescent2(A1, A2, T, Xs[iX], X0, radiances, pixel)
        minXs.append(minX)
        allErrors.append(currErrors)

    for i in range(len(Xs)):
        plt.subplot(1, 1, i + 1)
        plt.plot(np.arange(0, len(allErrors[i])), allErrors[i])

    plt.show()

    for i in range(len(Xs)):
        plt.subplot(1, 1, i + 1)

        layers = len(X) / len(const.species)

        xAxis = []
        for iSpec in range(len(const.species)):
            xAxis.append(const.species[iSpec][0])

        scalingFactors = []
        for iScale in range(const.changeableLayers, 0, -1):
            scalingFactorsAtIndex = []
            for iSpecBoundaries in range(layers, len(const.species) * layers + 1, layers):
                scalingFactorsAtIndex.append(minXs[i][iSpecBoundaries - iScale] / X[iSpecBoundaries - iScale])

            scalingFactors.append(scalingFactorsAtIndex)

        plotBars(xAxis, scalingFactors)

    plt.show()

    for i in range(len(Xs)):
        plt.subplot(1, 1, i + 1)
        conv = calcFlr2(A1, A2, T, minXs[i], X0, pixel)
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

xStartStop = (5742, 5745)
yStartStop = (580, 583)


def processPixel(pixels):
    if not isValidPixel(pixels[0], pixels[1]):
        return 0

    A1, A2, T, X = retrieveHeightAdjustedATX2(pixels[0], pixels[1])
    X0 = X[:]
    radiances = adjustMeasuredRadiances(pixels[0], pixels[1])
    minX, minError, errors = gradientDescent2(A1, A2, T, X, X0, radiances, pixels)
    return minX[len(minX) - 1] / X0[len(X0) - 1]


def analyzeImage():
    initializeConstants()
    methaneFactors = []
    nPixels = 1.0 * (xStartStop[1] - xStartStop[0]) * (yStartStop[1] - yStartStop[0])
    nCompleted = 0
    for pixelX in range(xStartStop[0], xStartStop[1]):
        for pixelY in range(yStartStop[0], yStartStop[1]):
            methaneFactors.append(processPixel((pixelX, pixelY)))
            nCompleted += 1
            print str(100 * nCompleted / nPixels) + "%"

        currentTime = time.localtime()
        print "Writing File..."
        np.savetxt("TextOutputs/CH4Factors" + time.strftime('%a-%d-%b-%Y-%H-%M-%S-GMT-', currentTime) + str(pixelX) +
                   "-" + str(pixelY) + ".txt", np.reshape(methaneFactors, (1 + pixelX - xStartStop[0], yStartStop[1] - yStartStop[0])), delimiter=", ",newline=" |\n")
        print "Done"

    return np.reshape(methaneFactors, (xStartStop[1] - xStartStop[0], yStartStop[1] - yStartStop[0]))

def analyzeImageSingleThread():
    initializeConstants()
    allBeforeDiffs = []
    allAfterDiffs = []
    allRatios = []
    for pixelX in range(xStartStop[0], xStartStop[1]):
        for pixelY in range(yStartStop[0], yStartStop[1]):
            pixels = (pixelX, pixelY)
            A1, A2, T, X = retrieveHeightAdjustedATX2(pixels[0], pixels[1])
            X0 = X[:]
            # A2, T2, X2, = retrieveHeightAdjustedATX2(pixels[0], pixels[1])
            radiances = adjustMeasuredRadiances(pixels[0], pixels[1])
            # Flr = calcFlr2(A1, A2, T, X, pixels)
            # Flr2 = calcFlr(A2, T2, X2, pixels)
            # plt.plot(const.wlAxis, Flr)
            # plt.plot(const.wlAxis, Flr2)
            # plt.plot(const.bands, radiances)
            # for i in range(30):
            #     radiances = adjustMeasuredRadiances(pixels[0] + i, pixels[1])
            #     plt.plot(const.bands, radiances)

            # plt.show()
            gradientDescentMultiAlphas2(A1, A2, T, X, X0, radiances, pixels)


    #         pixel = (pixelX, pixelY)
    #         A, T, X = retrieveHeightAdjustedATX(pixel[0], pixel[1])
    #         Flr = calcFlr(A, T, X, pixel)
    #         radiances = adjustMeasuredRadiances(pixelX, pixelY)
    #         plt.plot(const.wlAxis, Flr)
    #         plt.plot(const.bands, radiances)
    #
    #         diffsBefore = []
    #         for iRad in range(len(radiances)):
    #             diffsBefore.append(radiances[iRad] - Flr[const.wlAxisBandIndexes[iRad]])
    #
    #         Xclone = X[:]
    #         gradientDescent(A, T, Xclone, radiances, pixel)
    #
    #         diffsAfter = []
    #         FlrAfter = calcFlr(A, T, Xclone, pixel)
    #         plt.plot(const.wlAxis, FlrAfter)
    #         plt.show()
    #
    #         for iRad in range(len(radiances)):
    #             diffsAfter.append(radiances[iRad] - FlrAfter[const.wlAxisBandIndexes[iRad]])
    #
    #         allBeforeDiffs.append(diffsBefore)
    #         allAfterDiffs.append(diffsAfter)
    #         allRatios.append(np.divide(Xclone, X))
    #
    # plt.subplot(2, 1, 1)
    # plotBars(const.bands, allBeforeDiffs)
    # plt.subplot(2, 1, 2)
    # plotBars(const.bands, allAfterDiffs)
    # plt.show()
    # xAxis2 = np.arange(len(allRatios[0]))
    # plotBars(xAxis2, allRatios)
    # plt.show()


def analyzeAndPlotImage():
    methaneFactors = analyzeImage()
    heatmap = plt.imshow(methaneFactors, cmap='plasma', interpolation='nearest', vmin=0.8, vmax=2)
    plt.colorbar(heatmap)
    plt.show()
    # im = Image.open("AerialImage.jpeg")
    # im2 = im.crop((yStartStop[0], xStartStop[0], yStartStop[1] - 1, xStartStop[1] - 1))
    # im2.show()

# analyzeAndPlotImage()
analyzeImageSingleThread()
# analyzeImageSingleThread()
# methaneFactors = analyzeImageMultiThreaded()
# heatmap = plt.imshow(methaneFactors, cmap='plasma', interpolation='nearest', vmin=0, vmax=3)
# plt.colorbar(heatmap)
# plt.show()
# print methaneFactors

# analyzeImageSingleThread()
