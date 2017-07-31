import numpy as np
from hapi import *
import math as Math
import matplotlib.pyplot as plt
import DataReader as dataReader

def checkFlag(bit):
    return plotFlag >> bit & 0b1

def swapWnWl(WlWn):
    return (1.0 / WlWn) * 10000000


def plot1Matrix(rgWl, coef):
    if not checkFlag(1):
        return

    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[10 + i])

    plt.show()
    for i in range(10):
        plt.subplot(5, 2, i + 1)
        plt.plot(rgWl, coef[20 + i])

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


plotFlag = 0b111
FWHM = 5
wnStep = 0.01
wlStep = 0.01
wavelengthLow = 2278 - 4 * FWHM
wavelengthHigh = 2358 + 4 * FWHM
wnLow = swapWnWl(wavelengthHigh)
wnHigh = swapWnWl(wavelengthLow)
bands = 9
species = [['H20', 1, 1], ['N20', 4, 1], ['CH4', 6, 1]]
# Boundary is everything above the number ex: layer 2 (index 1) is 10->20
atmosLayersAltitude = []
    # [16.5, 12.5, 10, 8, 6, 4.5, 3, 2, 1, 0]
atmosLayerT = []
    # = dataReader.createTempProfile()
atmosLayerP = []
    # = [0.07, 0.135, 0.21, 0.3, 0.41, 0.51, 0.625, 0.735, 0.83, 0.94]


def initializeConstants(coordinates):
    dataReader.initializeDataReader()
    coordinates, atmosLayerP, atmosLayersAltitude, atmosLayerT, verticalProfileCH4, verticalProfileN2O, verticalProfileH2O = dataReader.getProfilesAtCoordinate(coordinates)
    wlAxis, I0, TrefMat = retrieveLinearlyAlignedI0andTref()
    return wlAxis, I0, TrefMat, verticalProfileH2O + verticalProfileN2O + verticalProfileCH4


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

    for iSpec in range(len(species)):
        fetch(species[iSpec][0], species[iSpec][1], species[iSpec][2], wnLow, wnHigh)

        for iLayer in range(len(atmosLayersAltitude)):
            nu, coef = absorptionCoefficient_Lorentz(SourceTables=species[iSpec][0], HITRAN_units=True,
                                                     Environment={'T': atmosLayerT[iLayer], 'p': atmosLayerP[iLayer]},
                                                     WavenumberStep=wnStep)

            # Outer layer of atmosphere is at 100 km
            previousBoundary = 100 if iLayer == 0 else atmosLayersAltitude[iLayer - 1]
            distance = 100000 * previousBoundary - atmosLayersAltitude[iLayer]
            coef *= distance

            TrefNusList.append(nu)
            TrefCoefsList.append(coef)

    return TrefNusList, TrefCoefsList


def retrieveLinearlyAlignedI0andTref():
    TrefCoefsFixedMatrix = []
    wlI0, valI0 = dataReader.readI0()
    xI0, yI0 = fixXAxis(wlI0, valI0, wavelengthLow, wavelengthHigh, wlStep)
    TrefNusList, TrefCoefsList = retrieveTref()
    plot2Matrix(TrefNusList, TrefCoefsList)

    for iList in range(len(TrefNusList)):
        rgWl = list(reversed(TrefNusList[iList]))
        rgCoefs = list(reversed(TrefCoefsList[iList]))
        for iWl in range(len(rgWl)):
            rgWl[iWl] = swapWnWl(rgWl[iWl])

        rgWlFixed, rgCoefsFixed = fixXAxis(rgWl, rgCoefs, wavelengthLow, wavelengthHigh, wlStep)
        TrefCoefsFixedMatrix.append(rgCoefsFixed)

    return xI0, yI0, TrefCoefsFixedMatrix


def retrieveA(altitude, SZA):
    SZAradians = np.deg2rad(SZA)
    arrA = np.zeros((len(atmosLayersAltitude) * len(species), 1))

    for iLayer in range(len(atmosLayersAltitude)):
        for iSpec in range(len(species)):
            # if the aircraft is entirely above the layer, it has the additional 1.
            arrA[iLayer + iSpec * len(atmosLayersAltitude)][0] = \
                1 / np.cos(SZAradians) if altitude < atmosLayersAltitude[iLayer] else 1 + 1 / np.cos(SZAradians)

    return np.matrix(arrA)


# def alignI0andT():
#     nuI0, coefI0 = retrieveI0()
#     nuT, coefT = retrieveTref()
#     nuI0List = nuI0.tolist()
#     coefI0List = coefI0.tolist()
#     nuTList = nuT.tolist()
#     coefTList = []
#     for arr in coefT:
#         coefTList.append(arr.tolist())
#
#     #TODO: This works, but is super naive/inefficient. It is not a bottleneck though, and thus not a priority.
#     while nuI0List[0] - nuTList[0] > 0.00001:
#         del nuTList[0]
#         for list in coefTList:
#             del list[0]
#
#     while nuTList[0] - nuI0List[0] > 0.00001:
#         del nuI0List[0]
#         del coefI0List[0]
#
#     while nuI0List[len(nuI0List) - 1] - nuTList[len(nuTList) - 1] > 0.00001:
#         del nuI0List[len(nuI0List) - 1]
#         del coefI0List[len(coefI0List) - 1]
#
#     while nuTList[len(nuTList) - 1] - nuI0List[len(nuI0List) - 1] > 0.00001:
#         del nuTList[len(nuTList) - 1]
#         for list in coefTList:
#             del list[len(list) - 1]
#
#     return np.array(nuI0List), np.array(coefI0List), np.array(nuTList), np.array(coefTList)

# The Solar Spectrum information is actually too high resolution, so we lower it to match the absorption line
# resolution by average all the values that round to the lower resolution axis
# def averageIdenticalXAxisValues(nuI0file, coefI0file):
#     nuI0corrected = []
#     coefI0corrected = []
#     iNuFile = 0
#     while iNuFile < len(nuI0file):
#         cIdenticalNu = 1
#         sum = coefI0file[iNuFile]
#         while iNuFile + cIdenticalNu < len(nuI0file) and abs(nuI0file[iNuFile] - nuI0file[iNuFile + cIdenticalNu]) < wnStep / 1000.0:
#             sum += coefI0file[iNuFile + cIdenticalNu]
#             cIdenticalNu += 1
#
#         nuI0corrected.append(nuI0file[iNuFile])
#         coefI0corrected.append(1.0 * sum / cIdenticalNu)
#         iNuFile += cIdenticalNu
#
#     return np.flipud(nuI0corrected), np.flipud(coefI0corrected)
#
# def fitTrefToMatrix(TrefNusList, TrefCoefsList, wnStep):
#     # First find the min and max over all Nus
#     min = TrefNusList[0][0]
#     max = TrefNusList[0][len(TrefNusList[0]) - 1]
#     for nuList in TrefNusList:
#         if(nuList[0] < min):
#             min = nuList[0]
#         if(nuList[len(nuList) - 1] > max):
#             max = nuList[len(nuList) - 1]
#
#     # Create the final matrix with correct dimensions (1.01 to account for double -> int precision errors)
#     nuMat = np.empty((len(TrefCoefsList), int((max - min) / wnStep + 1.01)))
#     coefMat = np.empty((len(TrefCoefsList), int((max - min) / wnStep + 1.01)))
#
#
#     for iList in range(len(TrefNusList)):
#         nuList = TrefNusList[iList]
#         nu = np.asarray(nuList)
#         coef = np.asarray(TrefCoefsList[iList])
#
#         # Add linearly to nu and zeros to coef in order to equalize the dimensions for all lines
#         if (nuList[0] > min):
#             prependNu = np.arange(min, nuList[0] - wnStep / 2.0, wnStep)
#             prependCoef = np.zeros(len(prependNu))
#             nu = np.concatenate((prependNu, nu))
#             coef = np.concatenate((prependCoef, coef))
#
#         # Adds to the end rather than the beginning
#         if (nuList[len(nuList) - 1] < max):
#             appendNu = np.arange(nuList[len(nuList) - 1] + wnStep, max + wnStep / 2.0, wnStep)
#             appendCoef = np.zeros(len(appendNu))
#             nu = np.concatenate((nu, appendNu))
#             coef = np.concatenate((coef, appendCoef))
#
#         nuMat[iList] = nu
#         coefMat[iList] = coef
#
#     #TODO: Technically don't need to create entire matrix if all the rows are identical
#
#     return nuMat[0], coefMat
