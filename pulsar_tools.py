#!/usr/bin/python

"""
A set of (potentially) handy functions for quickly calculating common 
physical properties of pulsars. These include:

    -- dispersion delay:                            disp_delay
    -- rotational energy:                           E_rot
    -- rotational energy loss:                      E_loss
    -- minimum surface magnetic field strength:     min_mag_field
    -- magnetic energy density:                     mag_energy_density
    -- characterstic age:                           char_age
    -- spin period/frequncy conversion:             freq_period_convert
    -- Goldreich-Julian polar cap plasma density:   polar_cap_plasma_density
    -- Goldreich-Julian polar cap radius:           polar_cap_radius
    -- light cylinder radius:                       light_cyl_radius

"""
from astropy.constants import c,m_e,e,pc,eps0
from math import pi,sqrt



def disp_delay(freq1=184.975, freq2=200.335, DM=0.0):
	"""
	Accepts two frequencies as input (low then high, in MHz) as well as a 
	DM (in pc/cm^3) and reports the delay in seconds.
	Default frequncies: 184.975 and 200.355 MHz, default DM=0.0 pc/cm^3
	"""

	# constant from analytic plasma dispersion
	A = e.value**2 / (8 * pi**2 * m_e.value * c.value * eps0.value)
	
	# convert into units: MHz^2 cm^3 s / pc
	A = A * (1.0/1.0e6)**2 * (1e6) * pc.to('m').value

	del_t = A * DM * ((freq1)**(-2) - (freq2)**(-2))

	return del_t



def E_rot(p=1.0, *args):
	"""
	Accepts a spin period (in s) and returns the approximate roational
	energy in J. 

	Also accepts an optional value for the moment of inertia I,
	else reverts to the standard estimate, I=1.0e46 kg/m^2.
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
	returns the approximate roational energy loss rate in J/s. 

	Also accepts optional value for moment of inertia I, else reverts to standard 
	estimate, I=1.0e46 kg/m^2.
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
	returns the characteristic age. Units: years.
	"""

	tau = p / (2 * pdot)

	return tau / (60 * 60 * 24 * 365)



def freq_period_convert(x=1.0, *args):
	"""
	Accepts a period or frequency and the corresponding derivatives (first and second).
	This function converts between period and frequncy (or vice versa). 
	If a derivative and/or second derivatives are provided, output is a list.

        Argument order is assumed to be: p, *p-dot, *p-dotdot OR f, *f-dot, *f-dotdot.
	Arguments marked with * are optional.
	"""

	if len(args) == 0:
		y = 1.0 / x

		return y

	elif len(args) == 1:
		y = 1.0 / x
		yd = -(args[0] / x**2)

		return [y, yd]

	elif len(args) == 2:
		y = 1.0 / x
		yd = -(args[0] / x**2)
		ydd = 2 * (args[0]**2 / x**3) - (args[1] / x**2)

		return [y, yd, ydd]

	else:
		print "freq_period_convert accepts 1, 2 or 3 args. You provided {0}.".format(len(args)+1)


def polar_cap_plasma_density(p=1.0, pdot=1.0e-12):
	"""
	Accepts the spin period and period derivative and returns the approximate polar cap plasma density in 
	the Goldreich-Julian model. Units: particles/cm^3.
	"""
	
	n = 7.0e10 * sqrt(1.0 / p) * sqrt(pdot / 1.0e-15)

	return n


def polar_cap_radius(p=1.0):
	"""
	Accepts the spin period and returns the approximate polar cap radius. Units: metres.
	"""
	
	r = 150 * sqrt(1 / p)

	return r


def light_cyl_radius(p=1.0):
	"""
	Accepts the spin period and returns the approximate light cylinder radius. Units: metres.
	"""

	rlc = c * p / (2 * pi)

	return rlc
	

	
	
if __name__ == '__main__':
	print "Testing dispersion delay (s) with freq1=100MHz, freq2=200MHz and DM=100pc/cm^3:",disp_delay(100.0, 200.0, 100.0)
	print "Testing E_rot (J), with P=0.0333s:",E_rot(0.0333)
	print "Testing E_loss (J/s), with P=0.0333s and Pdot=4.209e-13:",E_loss(0.0333, 4.209e-13)
	print "Testing B_surf (Gauss) with P=0.0333s and Pdot=4.209e-13:",min_mag_field(0.0333, 4.209e-13)
	print "Testing U_B with B=4e12 Gauss:",mag_energy_density(4e12)
	print "Testing char. age (years) with P=0.0333s and Pdot=4.209e-13:",char_age(0.0333, 4.209e-13)
	print "Testing frequency and period conversion..."
	print "f=29.946923 Hz, fdot=-3.77535E-10:",freq_period_convert(29.946923, -3.77535e-10)
	print "P=0.0333 s, Pdot=4.209e-13:",freq_period_convert(0.0333, 4.209e-13)
