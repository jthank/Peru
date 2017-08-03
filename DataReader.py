import os
from netCDF4 import Dataset
import numpy as np
import math as Math
import pyproj as pyproj
from spectral import *
import spectral.io.envi as envi
import GlobalConstants as const
import matplotlib.pyplot as plt

import h5py as h5


def barometricFormula(temperature, pressure):
    return -29.27175933 * temperature * Math.log(pressure)

def getProfilesAtCoordinate(inputCoordinate, coordinates, verticalProfilesCH4ByCoordinate, verticalProfilesN2OByCoordinate, verticalProfilesH2OByCoordinate, levelTemperaturesByCoordinate):
    minIndex = 0
    minValue = Math.sqrt((inputCoordinate[0] - coordinates[0][0])**2 + (inputCoordinate[1] - coordinates[0][1])**2)
    for iCoor in range(1, len(coordinates)):
        value = Math.sqrt((inputCoordinate[0] - coordinates[iCoor][0])**2 + (inputCoordinate[1] - coordinates[iCoor][1])**2)
        if value < minValue:
            minIndex = iCoor
            minValue = value

    const.verticalProfiles = verticalProfilesCH4ByCoordinate[minIndex].tolist() + verticalProfilesN2OByCoordinate[minIndex].tolist() + verticalProfilesH2OByCoordinate[minIndex].tolist()
    const.levelTemperatures = levelTemperaturesByCoordinate[minIndex]

    for iPress in range(len(const.levelPressures)):
        const.levelAltitudes.append(barometricFormula(const.levelTemperatures[iPress], const.levelPressures[iPress]))

def initializeDataReader():
    verticalProfilesCH4ByCoordinate = []
    verticalProfilesN2OByCoordinate = []
    verticalProfilesH2OByCoordinate = []
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
            verticalProfilesCH4ByCoordinate.append(rootgrp.variables['CH4'][iCrIS])
            verticalProfilesN2OByCoordinate.append(rootgrp.variables['N2O'][iCrIS])
            verticalProfilesH2OByCoordinate.append(rootgrp.variables['H2O'][iCrIS])

    coordinate = readAvirisData()
    getProfilesAtCoordinate(coordinate, const.coordinates, verticalProfilesCH4ByCoordinate, verticalProfilesN2OByCoordinate, verticalProfilesH2OByCoordinate, levelTemperaturesByCoordinate)
    rootgrp.close()


def selectBands():
    bandsArr = np.asarray(const.bands)
    const.AVIRIS_CenterIndexes = np.where(np.logical_and(bandsArr >= const.wlLow, bandsArr <= const.wlHigh))[0]
    const.bands = const.bands[const.AVIRIS_CenterIndexes[0]:const.AVIRIS_CenterIndexes[len(const.AVIRIS_CenterIndexes) - 1]]

def getLatLongCoordinates(rfl):
    # TODO: can read the region from relevant metadata and determine EPSG
    UTM11N = pyproj.Proj("+init=EPSG:32611")
    LonLat = pyproj.Proj("+init=EPSG:4326")
    x, y = rfl.metadata['map info'][3], rfl.metadata['map info'][4]
    lon, lat = pyproj.transform(UTM11N, LonLat, x, y)
    return (lon, lat)


def readAvirisData():
    img = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\imgHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\img')
    ort = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\ortHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\ort')
    glt = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\gltHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\glt')
    igm = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\igmHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG2\RDN\igm')
    rfl = envi.open('C:\Users\Jordan\Desktop\AVIRIS-NG2\RFL\imgHDR.hdr', 'C:\Users\Jordan\Desktop\AVIRIS-NG2\RFL\img')

    const.arrIMG = img.open_memmap()
    const.arrORT = ort.open_memmap()
    arrGLT = glt.open_memmap()
    arrIGM = igm.open_memmap()
    const.arrRFL = rfl.open_memmap()

    const.bands = img.bands.centers
    const.FWHM = float(rfl.metadata['fwhm'][const.AVIRIS_CenterIndexes[len(const.AVIRIS_CenterIndexes) / 2]])
    selectBands()
    return getLatLongCoordinates(rfl)

def readI0():
    nuI0file = []
    coefI0file = []
    fraunSpec = open("FraunSpecHR_2258-2378.txt", "r")

    for line in fraunSpec:
        pair = line.split()
        # 10.0 because the input data is in Angstroms not nm
        nuI0file.append(float(pair[0]) / 10.0)
        coefI0file.append(float(pair[1]))

    return nuI0file, coefI0file


# def readFile(password, verbose):
#     # if (len(sys.argv) != 2):
#     #   print "usage: "+sys.argv[0]+" [-q] password_on_RDA_webserver"
#     #   print "-q suppresses the progress message for each file that is downloaded"
#     #   sys.exit(1)
#
#     # passwd_idx=1
#     # verbose=True
#     # if (len(sys.argv) == 3 and sys.argv[1] == "-q"):
#     #   passwd_idx=2
#     #   verbose=False
#
#     cj=cookielib.MozillaCookieJar()
#     opener=urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
#
#     # check for existing cookies file and authenticate if necessary
#     do_authentication=False
#     if (os.path.isfile("auth.rda.ucar.edu")):
#       cj.load("auth.rda.ucar.edu",False,True)
#       for cookie in cj:
#         if (cookie.name == "sess" and cookie.is_expired()):
#           do_authentication=True
#     else:
#       do_authentication=True
#     if (do_authentication):
#       login=opener.open("https://rda.ucar.edu/cgi-bin/login","email=jhank@stanford.edu&password="+password+"&action=login")
#
#     # save the authentication cookies for future downloads
#     # NOTE! - cookies are saved for future sessions because overly-frequent authentication to our server can cause your data access to be blocked
#       cj.clear_session_cookies()
#       cj.save("auth.rda.ucar.edu",True,True)
#
#     # download the data file(s)
#     listoffiles=["3HRLY/2008/NARRsfc_200806_2030.tar"]
#     for file in listoffiles:
#       idx=file.rfind("/")
#       if (idx > 0):
#         ofile=file[idx+1:]
#       else:
#         ofile=file
#       if (verbose):
#         sys.stdout.write("downloading "+ofile+"...")
#         sys.stdout.flush()
#       infile=opener.open("http://rda.ucar.edu/data/ds608.0/"+file)
#       outfile=open(ofile,"wb")
#       outfile.write(infile.read())
#       outfile.close()
#       if (verbose):
#         sys.stdout.write("done.\n")