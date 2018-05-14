from astropy.units import Hz, MHz
from numpy import arange

def freqs(LO = 635):
	"""Frequency array in MHz for this LO setting."""
	return (arange(0,8192) * 1464.84375 * Hz + (LO * MHz * 2 + 144*  MHz)).to(MHz)