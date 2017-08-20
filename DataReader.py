from netCDF4 import Dataset
import numpy as np
import math as Math
import pyproj as pyproj
import spectral.io.envi as envi
import GlobalConstants as const
from hapi import *



def barometricFormula(temperature, pressure):
    return -29.27175933 * temperature * Math.log(pressure)


def getProfilesAtCoordinate(inputCoordinate, coordinates, allVerticalProfiles, levelTemperaturesByCoordinate):
    minIndex = 0
    minValue = Math.sqrt((inputCoordinate[0] - coordinates[0][0])**2 + (inputCoordinate[1] - coordinates[0][1])**2)
    for iCoor in range(1, len(coordinates)):
        value = Math.sqrt((inputCoordinate[0] - coordinates[iCoor][0])**2 + (inputCoordinate[1] - coordinates[iCoor][1])**2)
        if value < minValue:
            minIndex = iCoor
            minValue = value

    for iProf in range(len(allVerticalProfiles)):
        const.verticalProfiles += allVerticalProfiles[iProf][minIndex].tolist()

    const.levelTemperatures = levelTemperaturesByCoordinate[minIndex]

    for iPress in range(len(const.levelPressures)):
        const.levelAltitudes.append(barometricFormula(const.levelTemperatures[iPress], const.levelPressures[iPress]))


def initializeDataReader():
    allVerticalProfiles = []
    for iSpec in range(len(const.species)):
        allVerticalProfiles.append([])

    levelTemperaturesByCoordinate = []

    files = os.listdir('./NUCAPS_SB')
    rootgrp = Dataset('NUCAPS_SB/' + files[0], "r", format="NETCDF4")

    for iPress in range(rootgrp.variables['Pressure'][0].size):
        const.levelPressures.append(rootgrp.variables['Pressure'][0][iPress] / 1013.25)

    for filename in files:
        rootgrp = Dataset('NUCAPS_SB/' + filename, "r", format="NETCDF4")
        for iCrIS in range(rootgrp.variables['Longitude'].size):
            const.coordinates.append(
                (
                    rootgrp.variables['Longitude'][iCrIS],
                    rootgrp.variables['Latitude'][iCrIS]
                )
            )
            levelTemperaturesByCoordinate.append(rootgrp.variables['Temperature'][iCrIS])
            for iSpec in range(len(const.species)):
                allVerticalProfiles[iSpec].append(rootgrp.variables[const.species[iSpec][0]][iCrIS])

    coordinate = readAvirisData()
    getProfilesAtCoordinate(coordinate, const.coordinates, allVerticalProfiles, levelTemperaturesByCoordinate)
    rootgrp.close()


def selectBands():
    bandsArr = np.asarray(const.bands)
    const.AVIRIS_CenterIndexes = np.where(np.logical_and(bandsArr >= const.wlLow, bandsArr <= const.wlHigh))[0]
    const.bands = const.bands[const.AVIRIS_CenterIndexes[0]:const.AVIRIS_CenterIndexes[len(const.AVIRIS_CenterIndexes) - 1] + 1]


def getLatLongCoordinates(rfl):
    # TODO: can read the region from relevant metadata and determine EPSG
    UTM11N = pyproj.Proj("+init=EPSG:32611")
    LonLat = pyproj.Proj("+init=EPSG:4326")
    x, y = rfl.metadata['map info'][3], rfl.metadata['map info'][4]
    lon, lat = pyproj.transform(UTM11N, LonLat, x, y)
    return (-119.863672, 34.423762)
    # return (lon, lat)

def readAvirisData():
    img = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\imgHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\img')
    ort = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\ortHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\ort')
    glt = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\gltHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\glt')
    igm = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\igmHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RDN\igm')
    rfl = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\imgHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\img')
    H2O = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\H2OHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG-COP\RFL\H2O')
    const.arrH2O = H2O.open_memmap()
    const.arrIMG = img.open_memmap()
    const.arrORT = ort.open_memmap()
    const.arrGLT = glt.open_memmap()
    const.arrIGM = igm.open_memmap()
    const.arrRFL = rfl.open_memmap()

    # view = imshow(img, (29, 19, 9))
    const.bands = img.bands.centers
    selectBands()
    const.FWHM = float(rfl.metadata['fwhm'][const.AVIRIS_CenterIndexes[len(const.AVIRIS_CenterIndexes) / 2]])
    const.gaussian = SLIT_GAUSSIAN(np.arange(-1.5 * const.FWHM, 1.5 * const.FWHM + const.wlStep, const.wlStep), const.FWHM)
    return getLatLongCoordinates(rfl)

def readI0():
    nuI0file = []
    coefI0file = []
    fraunSpec = open("ETSI.dat", "r")

    for line in fraunSpec:
        pair = line.split(",")
        nuI0file.append(float(pair[0]))
        coefI0file.append(float(pair[1]) * float(pair[0]))

    return nuI0file, coefI0file
