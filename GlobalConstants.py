from Constants import swapWnWl
import math as Math


pixelX = 400
pixelY = 400

#2278
#2358
#2330
wlLow = 2278
wlHigh = 2330
FWHM = 5
wlLowPadded = wlLow - 4 * FWHM
wlHighPadded = wlHigh + 4 * FWHM

# All the same size
AVIRIS_CenterIndexes = []
bands = []
wlAxisBandIndexes = []

coordinates = []
levelPressures = []
levelAltitudes = []
levelTemperatures = []
verticalProfiles = []

wlAxis = []
TrefMat = []
I0 = []

A = []

plotFlag = 0b100

# Resolution up to 0.01
wnStep = 0.01
wlStep = 0.1

wnLow = swapWnWl(wlHighPadded)
wnHigh = swapWnWl(wlLowPadded)
species = [['H2O', 1, 1], ['N2O', 4, 1], ['CO', 5, 1], ['CH4', 6, 1]]

arrIMG = []
arrORT = []
arrRFL = []
arrIGM = []
arrH2O = []

changeableLayers = 8
alpha = 5 * Math.pow(10, 33)
iterations = 100
stopThreshold = 0.05