from hapi import *
import math as Math
from Constants import *
import matplotlib.pyplot as plt
import numpy as np

test, test2 = retrieveTref()

for iTest in range(30):
    plt.subplot(10, 3, iTest + 1)
    plt.plot(test[iTest], test2[iTest])

plt.show()

def calcFhr(A, I0, TrefNus, TrefCoef):
    print 'hi'