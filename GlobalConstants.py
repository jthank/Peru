from Constants import swapWnWl
import math as Math



FWHM = 5

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


species = [['H2O', 1, 1], ['N2O', 4, 1], ['CO', 5, 1], ['CH4', 6, 1]]

arrIMG = []
arrORT = []
arrRFL = []
arrIGM = []
arrH2O = []

#2278
#2358
#2330
wlLow = 2278
wlHigh = 2330

# Resolution up to 0.01
wnStep = 0.01
wlStep = 0.1
gaussian = []

changeableLayers = 3
alpha = 2 * Math.pow(10, 34)
iterations = 100
stopThreshold = 0.10

wlLowPadded = wlLow - 4 * FWHM
wlHighPadded = wlHigh + 4 * FWHM

wnLow = swapWnWl(wlHighPadded)
wnHigh = swapWnWl(wlLowPadded)
