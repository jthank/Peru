from hapi import *
from os import sys
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
from scipy import stats


def plotBars(xAxis, bars):
    barWidth = 0.9 / len(bars)
    indexes = np.arange(len(bars[0]))
    labels = ["L_97", "L_98", "L_99"]
    colors = ["black", "gray", "lightblue"]
    for iBars in range(len(bars)):
        plt.bar(indexes + iBars * barWidth, bars[iBars], barWidth, label=labels[iBars], color=colors[iBars])

    plt.legend()
    plt.xticks(indexes, xAxis, rotation='vertical')

def gaussianConvolution(Fhr):

    return np.convolve(Fhr, const.gaussian, mode="same") * const.wlStep


def calcFhr(A1, A2, T, X, X0, pixel):
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

def calcFlr(A1, A2, T, X, X0, pixel):
    plot1Matrix(const.wlAxis, T)
    Fhr = calcFhr(A1, A2, T, X, X0, pixel)
    conv = gaussianConvolution(Fhr)
    return Fhr, conv


def isValidPixel(pixelX, pixelY):
    if const.arrIMG[pixelX][pixelY][0] < 0:
        return False

    return True


def retrieveAMF(aircraftAlt, SZA, CZA, lowAlt, highAlt):
    cameraZenithFactor = 1 if aircraftAlt > highAlt else 0 if aircraftAlt < lowAlt else (aircraftAlt - lowAlt) / (highAlt - lowAlt)
    # replace 1 with (highAlt - lowAlt)
    return 1.0 / Math.cos(Math.radians(SZA)), cameraZenithFactor * (1.0 / Math.cos(Math.radians(CZA)))

def retrieveHeightAdjustedATX(pixelX, pixelY):
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
        AMF1, AMF2 = retrieveAMF(aircraftAlt, SZA, CZA, low, high)
        for iSpec in range(0, len(const.species)):
            A1[len(A1) / len(const.species) * iSpec + iAlt] = AMF1
            A2[len(A2) / len(const.species) * iSpec + iAlt] = AMF2


    return A1, A2, T, X


def getError(radiances, Flr):
    sum = 0
    for iBand in range(len(const.bands)):
        sum += abs(Flr[const.wlAxisBandIndexes[iBand]] - radiances[iBand])

    # print sum

    return sum


def applyGradient(A1, A2, T, X, X0, radiances, currErrors, pixel, minX, minError):
    Fhr, Flr = calcFlr(A1, A2, T, X, X0, pixel)
    error = getError(radiances, Flr)
    currErrors.append(error)

    if error < minError:
        minError = error
        minX = X[:]

    layers = len(X) / len(const.species)
    for iSpecBoundaries in range(layers, len(const.species) * layers + 1, layers):
        for iX in range(iSpecBoundaries - 1, iSpecBoundaries - const.changeableLayers - 1, -1):
            # g is the partial derivate of FLR
            g = np.multiply(np.multiply(T[iX], Fhr), -1 * (A1[iX] + A2[iX]))
            for iBand in range(len(const.bands)):
                a = 2 * g[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]] - 2 * radiances[iBand] * g[const.wlAxisBandIndexes[iBand]]
                b = Math.sqrt(radiances[iBand] * radiances[iBand] - 2 * radiances[iBand] * Flr[const.wlAxisBandIndexes[iBand]] + Flr[const.wlAxisBandIndexes[iBand]] * Flr[const.wlAxisBandIndexes[iBand]])
                X[iX] -= const.alpha * a / b

    return minX, minError

def gradientDescent(A1, A2, T, X, X0, radiances, pixel):
    errors = []
    minX = []
    minError = sys.float_info.max
    q = Queue.Queue()
    for i in range(5):
        q.put(sys.float_info.max)

    for i in range(const.iterations):
        minX, minError = applyGradient(A1, A2, T, X, X0, radiances, errors, pixel, minX, minError)
        if q.get() - minError < const.stopThreshold:
            return minX, minError, errors

        q.put(minError)

    return minX, minError, errors


def gradientDescentMultiAlphas(A1, A2, T, X, X0, radiances, pixel):
    iterations = 100
    base = Math.pow(10, 34)
    Xs = []
    alphas = [0.02 * base, 0.2 * base, 2 * base, 20 * base]
    allErrors = []
    minXs = []
    for i in range(len(alphas)):
        Xs.append(X[:])

    for iX in range(len(Xs)):
        const.alpha = alphas[iX]
        minX, minError, currErrors = gradientDescent(A1, A2, T, Xs[iX], X0, radiances, pixel)
        minXs.append(minX)
        allErrors.append(currErrors)

    for i in range(len(Xs)):
        plt.subplot(2, 2, i + 1)
        plt.plot(np.arange(0, len(allErrors[i])), allErrors[i])

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
        Fhr, conv = calcFlr(A1, A2, T, minXs[i], X0, pixel)
        plt.plot(const.wlAxis, conv, label='Modeled Irradiance')
        plt.plot(const.bands, radiances, label='Gathered Irradiance')
        plt.legend()
    plt.show()


def adjustMeasuredRadiances(pixelX, pixelY):
    radiances = []
    for iC in range(len(const.AVIRIS_CenterIndexes)):
        radiances.append(
            (0.01 * Math.pi * const.arrIMG[pixelX][pixelY][const.AVIRIS_CenterIndexes[iC]] * const.bands[iC])
            / (const.arrRFL[pixelX][pixelY][const.AVIRIS_CenterIndexes[iC]] * Math.cos(Math.radians(const.arrORT[pixelX][pixelY][4]))))

    return radiances


def loadText(path):
    return np.loadtxt(path, delimiter=", ",)


def stitchImages(dirPath):
    files = os.listdir(dirPath)
    rgFactors = []
    rgErrors = []
    rgBounds = []
    for path in files:
        dashSplit = path.split("-")
        if dashSplit[1] == "Factors":
            rgFactors.append(loadText(dirPath + "/" + path))
            rgBounds.append(((int(dashSplit[2]), int(dashSplit[3])), (int(dashSplit[4]), int(dashSplit[5]))))

        else:
            rgErrors.append(loadText(dirPath + "/" + path))


    xMin = sys.maxint
    xMax = 0
    yMin = sys.maxint
    yMax = 0

    for bound in rgBounds:
        if bound[0][0] < xMin:
            xMin = bound[0][0]

        if bound[0][1] > xMax:
            xMax = bound[0][1]

        if bound[1][0] < yMin:
            yMin = bound[1][0]

        if bound[1][1] > yMax:
            yMax = bound[1][1]

    stitchedMap = np.zeros((xMax - xMin, yMax - yMin))
    stitchedErrors = np.zeros((xMax - xMin, yMax - yMin))

    for iImage in range(len(rgFactors)):
        xStart = rgBounds[iImage][0][0] - xMin
        yStart = rgBounds[iImage][1][0] - yMin
        stitchedMap[xStart: xStart + rgBounds[iImage][0][1] - rgBounds[iImage][0][0], yStart:yStart + rgBounds[iImage][1][1] - rgBounds[iImage][1][0]] = rgFactors[iImage]
        stitchedErrors[xStart: xStart + rgBounds[iImage][0][1] - rgBounds[iImage][0][0], yStart:yStart + rgBounds[iImage][1][1] - rgBounds[iImage][1][0]] = rgErrors[iImage]

    plt.subplot(1,3,3)
    errorHeatmap = plt.imshow(np.multiply(100,stitchedErrors), cmap='Purples', interpolation='nearest')
    plt.colorbar(errorHeatmap)

    plt.subplot(1,3,2)
    factorheatmap = plt.imshow(stitchedMap, cmap='plasma', interpolation='nearest', vmin=1, vmax=2.2)

    plt.colorbar(factorheatmap)
    plt.subplot(1,3,1)
    im = Image.open("AerialImage.jpeg")
    im2 = im.crop((yMin, xMin, yMax, xMax))
    plt.imshow(im2)
    plt.show()
    plotReflectionCorrelation(stitchedMap, (xMin, xMax), (yMin, yMax), 0)
    plotReflectionCorrelation(stitchedMap, (xMin, xMax), (yMin, yMax), 1.75)


xStartStop = (2070, 2100)
yStartStop = (315, 369)


def plotReflectionCorrelation(methaneScalingFactors, xSS, ySS, constrainValue):
    reflectionFactors = []
    flatFactors = []
    for pixelX in range(xSS[0], xSS[1]):
        for pixelY in range(ySS[0], ySS[1]):
            if methaneScalingFactors[pixelX - xSS[0]][pixelY - ySS[0]] <= constrainValue:
                continue

            sum = 0
            count = 0

            rng = const.AVIRIS_CenterIndexes
            if len(rng) == 0:
                rng = range(386, 397)
                rfl = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\imgHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\img')
                const.arrRFL = rfl.open_memmap()

            for iC in rng:
                sum += const.arrRFL[pixelX][pixelY][iC]
                count += 1

            flatFactors.append(methaneScalingFactors[pixelX - xSS[0]][pixelY - ySS[0]])
            reflectionFactors.append(sum / count)

    flatFactors = np.asarray(flatFactors)
    plt.plot(flatFactors, reflectionFactors, '.')
    slope, intercept, r_value, p_value, std_err = stats.linregress(flatFactors, reflectionFactors)
    plt.plot(flatFactors, intercept + slope * flatFactors, 'r')
    # plt.text(1, -0.02, "R^2 = " + str(r_value ** 2))
    print r_value ** 2
    plt.show()

def processPixel(pixels):
    if not isValidPixel(pixels[0], pixels[1]):
        return 1, 0

    A1, A2, T, X = retrieveHeightAdjustedATX(pixels[0], pixels[1])
    X0 = X[:]
    radiances = adjustMeasuredRadiances(pixels[0], pixels[1])
    minX, minError, errors = gradientDescent(A1, A2, T, X, X0, radiances, pixels)
    return minError / np.sum(radiances), minX[len(minX) - 1] / X0[len(X0) - 1]

def getMinimumScalingFactor():
    print "Finding minimum scaling factor"
    normalPixels = [(5360, 350), (5560, 650), (6040, 500), (6130, 680)]
    radius = 2
    minScalingFactor = sys.float_info.max

    for pixel in normalPixels:
        for pixelX in range(pixel[0], pixel[0] + radius):
            for pixelY in range(pixel[1], pixel[1] + radius):
                error, scalingFactor = processPixel((pixelX, pixelY))
                if scalingFactor < minScalingFactor:
                    minScalingFactor = scalingFactor

    print "Minimum scaling factor found: " + str(minScalingFactor)
    return minScalingFactor


def analyzeImage():
    initializeConstants()
    minScalingFactor = 1.62
    # minScalingFactor = getMinimumScalingFactor()
    methaneFactors = []
    averagePercentErrors = []
    nPixels = 1.0 * (xStartStop[1] - xStartStop[0]) * (yStartStop[1] - yStartStop[0])
    nCompleted = 0
    for pixelX in range(xStartStop[0], xStartStop[1]):
        for pixelY in range(yStartStop[0], yStartStop[1]):
            averagePercentError, methaneScalingFactor = processPixel((pixelX, pixelY))
            methaneFactors.append(methaneScalingFactor / minScalingFactor)
            averagePercentErrors.append(averagePercentError)
            nCompleted += 1
            print str(100 * nCompleted / nPixels) + "%"

        currentTime = time.localtime()
        print "Writing Files..."
        np.savetxt("TextOutputs/" + str(pixelX - xStartStop[0]) + "-Factors-" + str(xStartStop[0]) + "-" + str(xStartStop[1]) + "-" + str(yStartStop[0]) + "-" + str(yStartStop[1]) + "-"
                   + ".txt", np.reshape(methaneFactors, (1 + pixelX - xStartStop[0], yStartStop[1] - yStartStop[0])), delimiter=", ", newline="\n")

        np.savetxt("TextOutputs/" + str(pixelX - xStartStop[0]) + "-Errors-" + str(xStartStop[0]) + "-" + str(xStartStop[1]) + "-" + str(yStartStop[0]) + "-" + str(yStartStop[1]) + "-"
                   + ".txt", np.reshape(averagePercentErrors, (1 + pixelX - xStartStop[0], yStartStop[1] - yStartStop[0])), delimiter=", ", newline="\n")

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
            A1, A2, T, X = retrieveHeightAdjustedATX(pixels[0], pixels[1])
            X0 = X[:]
            radiances = adjustMeasuredRadiances(pixels[0], pixels[1])
            gradientDescentMultiAlphas(A1, A2, T, X, X0, radiances, pixels)


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
    heatmap = plt.imshow(methaneFactors, cmap='plasma', interpolation='nearest', vmin=1, vmax=2.5)
    plt.colorbar(heatmap)
    plt.show()


# stitchImages("./TextOutputs/Final")
analyzeAndPlotImage()
# analyzeImageSingleThread()
