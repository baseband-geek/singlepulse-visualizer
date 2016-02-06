#!/usr/bin/python

"""
A set of (potentially) handy functions for quickly calculating common 
physical properties of pulsars. These include:

    -- dispersion delay:                           disp_delay
    -- rotational energy:                          E_rot
    -- rotational energy loss:                     E_loss
    -- minimum surface magnetic field strength:    min_mag_field
    -- magnetic energy density:                    mag_energy_density
    -- characterstic age:                          char_age

"""

from math import pi,sqrt



def disp_delay(freq1=0.184975, freq2=0.200335, DM=0):
	"""
	Accepts two frequencies as input (low then high in GHz) as well as a 
	DM (in pc/cm^3) and reports the delay in milliseconds.
	"""
	
	del_t = 4.148808 * DM * (freq1**(-2) - freq2**(-2))

	return del_t



def E_rot(p=1.0, *args):
	"""
	Accepts a spin period (in s) and returns the approximate roational
	energy in J. Also accepts an optional value for the moment of inertia I,
	else reverts to the standard estimate.
	"""
	
	if len(args) > 0:
		I = args[0]
	else:
		I = 1.0e46  # standard estiamte for pulsar moment of inertia
	
	e_rot = (2 * pi**2) * (I / p**2)

	return e_rot



def E_loss(p=1.0, pdot=1e-12, *args):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and
	returns the approximate roational energy loss rate in J/s. Also accepts
	optional value for moment of inertia I, else reverts to standard 
	estimate.
	"""
	
	if len(args) > 0:
		I = args[0]
	else:
		I = 1.0e46
	
	e_loss = (-4 * pi**2) * (I * pdot / p**3)

	return e_loss



def min_mag_field(p=1.0, pdot=1.0e-12):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and 
	returns the minimum surface magnetic field strength in Gauss. 
	
	Conversion to Tesla is: B(tesla) = 10^4 * B(Gauss)
	"""
		
	b_min = (3.2e19) * sqrt(p * pdot)

	return b_min



def mag_energy_density(b=1e12):
	"""
	Accepts a magnetic field strength in Gauss and returns the theoretical
	energy density in J/m^3. 

	Conversion to erg/cm^3 is: U(erg/cm^3) = U(J/m^3) / 10.
	"""

	u = 10 * b**2 / (8 * pi)

	return u



def char_age(p=1.0, pdot=1e-12):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and 
	returns the characteristic age (in years).
	"""

	tau = p / (2 * pdot)

	return tau / (60 * 60 * 24 * 365)




	
	
if __name__ == '__main__':
	print "Testing dispersion delay (ms) with freq1=100MHz, freq2=200MHz and DM=100pc/cm^3:",disp_delay(0.1, 0.2, 100)
	print "Testing E_rot (J), with P=0.0333s:",E_rot(0.0333)
	print "Testing E_loss (J/s), with P=0.0333s and Pdot=4.209e-13:",E_loss(0.0333, 4.209e-13)
	print "Testing B_surf (Gauss) with P=0.0333s and Pdot=4.209e-13:",min_mag_field(0.0333, 4.209e-13)
	print "Testing U_B with B=4e12 Gauss:",mag_energy_density(4e12)
	print "Testing char. age (years) with P=0.0333s and Pdot=4.209e-13:",char_age(0.0333, 4.209e-13)
