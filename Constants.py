import numpy as np
from hapi import *
import math as Math
import matplotlib.pyplot as plt
import DataReader as dataReader
import GlobalConstants as const


def checkFlag(bit):
    return const.plotFlag >> bit & 0b1

def swapWnWl(WlWn):
    return (1.0 / WlWn) * 10000000


def plot1Matrix(rgWl, coef):
    if not checkFlag(1):
        return

    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[10 * i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[100 + 10 * i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[200 + 9 * i])

    plt.show()


def plot2Matrix(nu, coef):
    if not checkFlag(0):
        return

    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(nu[i],coef[i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(nu[10+i], coef[10 + i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(nu[20+i],coef[20 + i])

    plt.show()



def initializeConstants():
    dataReader.initializeDataReader()
    setLinearlyAlignedI0andTref()


# rounds both numbers to the nearest multiple of step, and then returns:
# -1 if one is 'rounded less than' two
# 0 if one is 'rounded equal to' two
# 1 if one is 'rounded greater than' two
def roundedEquals(one, two, step):
    oneRounded = round((1.0 * one) / step) * step
    twoRounded = round((1.0 * two) / step) * step
    if abs(twoRounded - oneRounded) < step / 2:
        return 0

    if oneRounded < twoRounded:
        return -1

    return 1


def fixXAxis(nu, coef, wlLow, wlHigh, step):
    nuFixed = np.arange(wlLow, wlHigh + step / 2, step)
    coefFixed = np.zeros(len(nuFixed))
    iFixed = 0
    iOrig = 0
    counter = 0
    sum = 0
    while iFixed < len(nuFixed) and iOrig < len(nu):
        compare = roundedEquals(nu[iOrig], nuFixed[iFixed], step)
        if compare < 0:
            iOrig += 1

        elif compare > 0:
            if counter == 0:
                coefFixed[iFixed] = 0

            else:
                coefFixed[iFixed] = (1.0 * sum) / counter

            sum = 0
            counter = 0
            iFixed += 1

        else:
            counter += 1
            sum += coef[iOrig]
            iOrig += 1

    if counter != 0 and iFixed < len(nuFixed):
        coefFixed[iFixed] = (1.0 * sum) / counter

    return nuFixed, coefFixed


def retrieveTref():
    TrefNusList = []
    TrefCoefsList = []

    for iSpec in range(len(const.species)):
        fetch(const.species[iSpec][0], const.species[iSpec][1], const.species[iSpec][2], const.wnLow, const.wnHigh)

        for iLayer in range(len(const.levelAltitudes)):
            nu, coef = absorptionCoefficient_Lorentz(SourceTables=const.species[iSpec][0], HITRAN_units=True,
                                                     Environment={'T': const.levelTemperatures[iLayer], 'p': const.levelPressures[iLayer]},
                                                     WavenumberStep=const.wnStep)

            # Outer layer of atmosphere is at 100 km
            # previousBoundary = 100000 if iLayer == 0 else atmosLayersAltitude[iLayer - 1]
            # distance = 100 * (previousBoundary - atmosLayersAltitude[iLayer])
            # coef *= distance

            TrefNusList.append(nu)
            TrefCoefsList.append(coef)

    return TrefNusList, TrefCoefsList

def getIncidentIntensity(iWl, pixel):
    wavelengths = [2262.5, 2282.5, 2302.5, 2322.5, 2342.5, 2362.5]
    I0Factor = [.0238207, .0185019, .0161892, .0174742, .0145124, .0163191]
    SZA = Math.radians(const.arrORT[pixel[0]][pixel[1]][4])
    FraunFactor = Math.cos(SZA) * const.I0[iWl]
    print FraunFactor, const.wlAxis[iWl]
    index = 0
    while index < len(wavelengths) and const.wlAxis[iWl] > wavelengths[index]:
        index += 1

    if index == len(wavelengths):
        return I0Factor[len(I0Factor) - 1] * FraunFactor

    if index == 0:
        return I0Factor[0] * FraunFactor

    return FraunFactor * ((const.wlAxis[iWl] - wavelengths[index - 1]) / (wavelengths[index] - wavelengths[index - 1]) * (I0Factor[index] - I0Factor[index - 1]) + I0Factor[index - 1])

def makeI0Absolute(yI0):
    # Added the 0 and 5000 term to avoid special casing an index past the end of the array and the start
    wavelengths = [0, 2262.5, 2282.5, 2302.5, 2322.5, 2342.5, 2362.5, 5000]
    I0Factor = [.0238207, .0238207, .0185019, .0161892, .0174742, .0145124, .0163191, .0163191]
    index = 1
    for iWl in range(len(yI0)):
        if const.wlAxis[iWl] >= wavelengths[index]:
            index += 1

        yI0[iWl] = 10000 * ((const.wlAxis[iWl] - wavelengths[index - 1]) / (wavelengths[index] - wavelengths[index - 1]) * (I0Factor[index] - I0Factor[index - 1]) + I0Factor[index - 1])


def convertAxisToWlAndFixIntervals(rgWn, rgCoefs):
    rgWl = list(reversed(rgWn))
    rgCoefs = list(reversed(rgCoefs))
    for iWl in range(len(rgWl)):
        rgWl[iWl] = swapWnWl(rgWl[iWl])

    return fixXAxis(rgWl, rgCoefs, const.wlLowPadded, const.wlHighPadded, const.wlStep)


def setLinearlyAlignedI0andTref():
    TrefCoefsFixedMatrix = []
    wnI0, valI0 = dataReader.readI0()
    xI0, yI0 = convertAxisToWlAndFixIntervals(wnI0, valI0)
    TrefNusList, TrefCoefsList = retrieveTref()
    plot2Matrix(TrefNusList, TrefCoefsList)

    for iList in range(len(TrefNusList)):
        rgWlFixed, rgCoefsFixed = convertAxisToWlAndFixIntervals(TrefNusList[iList], TrefCoefsList[iList])
        TrefCoefsFixedMatrix.append(rgCoefsFixed)

    const.wlAxis = xI0
    for wl in const.bands:
        const.wlAxisBandIndexes.append(int(round((wl - const.wlAxis[0]) / const.wlStep)))

    # makeI0Absolute(yI0)
    const.I0 = yI0
    const.TrefMat = TrefCoefsFixedMatrix


def retrieveA(altitude, SZA):
    SZAradians = np.deg2rad(SZA)
    const.A = np.empty(len(const.levelAltitudes) * len(const.species))

    for iLayer in range(len(const.levelAltitudes)):
        for iSpec in range(len(const.species)):
            # if the aircraft is entirely above the layer, it has the additional 1.
            const.A[iLayer + iSpec * len(const.levelAltitudes)] = 1 / np.cos(SZAradians) if altitude < const.levelAltitudes[iLayer] else 1 + 1 / np.cos(SZAradians)
