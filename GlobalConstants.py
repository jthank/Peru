from Constants import swapWnWl

wlLow = 2278
wlHigh = 2358
FWHM = 5
wlLowPadded = wlLow - 4 * FWHM
wlHighPadded = wlHigh + 4 * FWHM
AVIRIS_CenterIndexes = []
bands = []
coordinates = []
levelPressures = []
levelAltitudes = []
levelTemperatures = []

wlAxis = []
TrefMat = []
I0 = []

A = []

plotFlag = 0b100
wnStep = 0.01
wlStep = 0.01

wnLow = swapWnWl(wlHighPadded)
wnHigh = swapWnWl(wlLowPadded)
species = [['H20', 1, 1], ['N20', 4, 1], ['CH4', 6, 1]]

arrIMG = []
arrORT = []
arrRFL = []
