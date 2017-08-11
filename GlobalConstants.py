from Constants import swapWnWl
pixelX = 400
pixelY = 400


wlLow = 2278
wlHigh = 2358
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
wnStep = 0.01
wlStep = 0.01

wnLow = swapWnWl(wlHighPadded)
wnHigh = swapWnWl(wlLowPadded)
#['CO', 5, 1],
species = [['H2O', 1, 1], ['N2O', 4, 1],  ['CH4', 6, 1]]

arrIMG = []
arrORT = []
arrRFL = []
arrIGM = []
arrH2O = []