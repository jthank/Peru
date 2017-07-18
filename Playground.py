from hapi import *
import math as Math
from Constants import *
import matplotlib.pyplot as plt
import numpy as np

fetch('CH4',6,1,4281,4292)

for i in range(len(atmosLayerT)):
    nu,coef = absorptionCoefficient_Lorentz(SourceTables='CH4',  HITRAN_units=False, Environment={'T':atmosLayerT[i],'p':atmosLayerP[i]})
    nu,coef = absorptionSpectrum(nu, coef)

    for j in range(len(coef)):
        nu[j] = Wn2Wl(nu[j])
        coef[j] *= Math.log(10, Math.e)

    plt.subplot(5,2,i+1)
    plt.plot(nu, coef)

plt.show()

