import numpy as np

__all__ = [
    "flux2mag",
    "mag2flux"
]

def flux2mag(flux, f0):
    return -2.5*np.log10(flux / f0)

def mag2flux(mag, f0):
    return 10**(mag/(-2.5)) * f0