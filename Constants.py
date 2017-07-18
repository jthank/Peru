import numpy as np

# Boundary is everything above the number ex: layer 2 (index 1) is 10->20
atmosLayersAltitude = [20, 10, 7, 6, 5, 4, 3, 2, 1, 0]
atmosLayerT = [300, 299, 298, 297, 296, 295, 294, 293, 292, 291]
atmosLayerP = [1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]


def retrieveA(layers, species, altitude, SZA):
    SZAradians = np.deg2rad(SZA)
    arrA = np.zeros((layers * species, 1))

    for iLayer in range(len(atmosLayersAltitude)):
        for iSpec in range(species):
            # if the plane is entirely above the layer, it has the additional 1.
            arrA[iLayer + iSpec * layers][0] = 1 / np.cos(SZAradians) if altitude < atmosLayersAltitude[iLayer] else 1 + 1 / np.cos(SZAradians)

    return np.matrix(arrA)

def Wn2Wl(wavenumber):
    return (1 / wavenumber) * 10000000