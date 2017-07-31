import sys
import os
import urllib2
import cookielib
from netCDF4 import Dataset, num2date, date2num, MFDataset
import numpy as np
import math as Math
from datetime import datetime, timedelta

import h5py as h5

coordinates = []
pressures = []
altitudes = []
verticalProfilesTemperatures = []
verticalProfilesCH4 = []
verticalProfilesN2O = []
verticalProfilesH2O = []
# fileExtensions = ["70", "150", "200", "300", "400", "500", "600", "700", "850", "925"]

def barometricFormula(temperature, pressure):
    return -29.27175933 * temperature * Math.log(pressure)

def getProfilesAtCoordinate(inputCoordinate):
    minIndex = 0
    minValue = Math.sqrt((inputCoordinate[0] - coordinates[0][0])**2 + (inputCoordinate[1] - coordinates[0][1])**2)
    for iCoor in range(1, len(coordinates)):
        value = Math.sqrt((inputCoordinate[0] - coordinates[iCoor][0])**2 + (inputCoordinate[1] - coordinates[iCoor][1])**2)
        if value < minValue:
            minIndex = iCoor
            minValue = value

    temperatures = verticalProfilesTemperatures[minIndex]
    for iPress in range(len(pressures)):
        altitudes.append(barometricFormula(temperatures[iPress], pressures[iPress]))

    return coordinates, pressures, altitudes, temperatures, verticalProfilesCH4[minIndex], verticalProfilesN2O[minIndex], verticalProfilesH2O[minIndex]

def initializeDataReader():
    files = os.listdir('./Rename')
    rootgrp = Dataset('Rename/' + files[0], "r", format="NETCDF4")

    for iPress in range(rootgrp.variables['Pressure'][0].size):
        pressures.append(rootgrp.variables['Pressure'][0][iPress] / 1013.25)

    for filename in files:
        rootgrp = Dataset('Rename/' + filename, "r", format="NETCDF4")
        for iCrIS in range(rootgrp.variables['Longitude'].size):
            coordinates.append(
                (
                    rootgrp.variables['Longitude'][iCrIS],
                    rootgrp.variables['Latitude'][iCrIS]
                )
            )
            verticalProfilesTemperatures.append(rootgrp.variables['Temperature'][iCrIS])
            verticalProfilesCH4.append(rootgrp.variables['CH4'][iCrIS])
            verticalProfilesN2O.append(rootgrp.variables['N2O'][iCrIS])
            verticalProfilesH2O.append(rootgrp.variables['H2O'][iCrIS])

    rootgrp.close()

def createAtmosphericLayers():
    layers = np.empty(len(fileExtensions))
    for iExt in range(len(fileExtensions)):
        rootgrp = Dataset("GeopotentialProfilesPeru/" + fileExtensions[iExt] + ".nc", "r", format="NETCDF4")
        layers[iExt] = rootgrp.variables.values()[0][0][0][0]
        rootgrp.close()

    return layers

def createTempProfile():
    temperatures = np.empty(len(fileExtensions))
    for iExt in range(len(fileExtensions)):
        rootgrp = Dataset("TemperatureProfilesPeru/" + fileExtensions[iExt] + ".nc", "r", format="NETCDF4")
        temperatures[iExt] = rootgrp.variables.values()[0][0][0][0]
        rootgrp.close()

    return temperatures


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


def readFile(password, verbose):
    # if (len(sys.argv) != 2):
    #   print "usage: "+sys.argv[0]+" [-q] password_on_RDA_webserver"
    #   print "-q suppresses the progress message for each file that is downloaded"
    #   sys.exit(1)

    # passwd_idx=1
    # verbose=True
    # if (len(sys.argv) == 3 and sys.argv[1] == "-q"):
    #   passwd_idx=2
    #   verbose=False

    cj=cookielib.MozillaCookieJar()
    opener=urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))

    # check for existing cookies file and authenticate if necessary
    do_authentication=False
    if (os.path.isfile("auth.rda.ucar.edu")):
      cj.load("auth.rda.ucar.edu",False,True)
      for cookie in cj:
        if (cookie.name == "sess" and cookie.is_expired()):
          do_authentication=True
    else:
      do_authentication=True
    if (do_authentication):
      login=opener.open("https://rda.ucar.edu/cgi-bin/login","email=jhank@stanford.edu&password="+password+"&action=login")

    # save the authentication cookies for future downloads
    # NOTE! - cookies are saved for future sessions because overly-frequent authentication to our server can cause your data access to be blocked
      cj.clear_session_cookies()
      cj.save("auth.rda.ucar.edu",True,True)

    # download the data file(s)
    listoffiles=["3HRLY/2008/NARRsfc_200806_2030.tar"]
    for file in listoffiles:
      idx=file.rfind("/")
      if (idx > 0):
        ofile=file[idx+1:]
      else:
        ofile=file
      if (verbose):
        sys.stdout.write("downloading "+ofile+"...")
        sys.stdout.flush()
      infile=opener.open("http://rda.ucar.edu/data/ds608.0/"+file)
      outfile=open(ofile,"wb")
      outfile.write(infile.read())
      outfile.close()
      if (verbose):
        sys.stdout.write("done.\n")