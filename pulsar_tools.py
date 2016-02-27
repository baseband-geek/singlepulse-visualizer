#!/usr/bin/python

"""
A set of (potentially) handy functions for quickly calculating common 
physical properties of pulsars. These include:

    -- spin period/frequncy conversion:             freq_period_convert
    -- dispersion delay:                            disp_delay
    -- rotational energy:                           E_rot
    -- rotational energy loss:                      E_loss
    -- minimum surface magnetic field strength:     min_mag_field
    -- mean parallel B component (LOS):             mean_parallel_B
    -- magnetic energy density:                     mag_energy_density
    -- characterstic age:                           char_age
    -- braking index:                               braking_index
    -- Goldreich-Julian polar cap plasma density:   polar_cap_plasma_density
    -- Goldreich-Julian polar cap radius:           polar_cap_radius
    -- light cylinder (LC) radius:                  light_cyl_radius
    -- magnetic field strength at LC radius:        light_cyl_B

"""
from astropy.constants import c,m_e,e,pc,eps0
from math import pi,sqrt



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



def disp_delay(freq1=0.184975, freq2=0.200335, DM=0.0):
	"""
	Accepts two frequencies as input (low then high, in GHz) as well as a 
	DM (in pc/cm^3) and reports the delay in milliseconds (ms).
	Default frequncies: 184.975 and 200.355 MHz, default DM=0.0 pc/cm^3
	"""

	# constant from analytic cold plasma dispersion
	A = e.value**2 / (8 * pi**2 * m_e.value * c.value * eps0.value)
	
	# convert into units: GHz^2 cm^3 ms / pc
	A = A * (1.0/1.0e9)**2 * (1e6) * (1e3) *  pc.to('m').value

	del_t = A * DM * ((freq1)**(-2) - (freq2)**(-2))

	return del_t



def E_rot(p0=1.0, *args):
	"""
	Accepts a spin period (in s) and returns the approximate roational
	energy in J. 

	Also accepts an optional value for the moment of inertia I,
	else reverts to the standard estimate, I=1.0e38 kg m^2 (Lorimer & Kramer, 2004).
	"""
	
	if len(args) > 0:
		I = args[0]
	else:
		I = 1.0e38 
	
	e_rot = (2 * pi**2) * (I / p0**2)

	return e_rot



def E_loss(p0=1.0, p1=1e-12, *args):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and
	returns the approximate roational energy loss rate in J/s. 

	Also accepts optional value for moment of inertia I, else reverts to standard 
	estimate, I=1.0e38 kg m^2 (Lorimer & Kramer, 2004).
	"""
	
	if len(args) > 0:
		I = args[0]
	else:
		I = 1.0e38
	
	e_loss = (-4 * pi**2) * (I * p1 / p0**3)

	return e_loss



def min_mag_field(p0=1.0, p1=1.0e-12):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and 
	returns the minimum surface magnetic field strength in Gauss. 
	
	Conversion to Tesla is: B(tesla) = 10^4 * B(Gauss)
	"""
		
	b_min = (3.2e19) * sqrt(p0 * p1)

	return b_min



def mean_parallel_B(rm=0.0, dm=100.0):
	"""
	Accepts a rotation measure (in rad/m^2) and dispersion measure (in pc/cm^3) and 
	returns the approximate mean value for the paralelle magnetic field compoent along 
	the line of sight, in Gauss."
	"""

	b_los = 1.23 * (rm / dm)

	return b_los / 1e6



def mag_energy_density(b=1e12):
	"""
	Accepts a magnetic field strength in Gauss and returns the theoretical
	energy density in J/m^3. 

	Conversion to erg/cm^3 is: U(erg/cm^3) = U(J/m^3) / 10.
	"""

	u = 10 * b**2 / (8 * pi)

	return u



def char_age(p0=1.0, p1=1e-12):
	"""
	Accepts a spin period (in s) and spin period derivative (in s/s) and 
	returns the characteristic age. Units: years.
	"""

	tau = p0 / (2 * p1)

	return tau / (60 * 60 * 24 * 365)



def braking_index(p0=1.0, p1=1e-12, p2=1e-20):
	"""
	Accepts a spin period, pdot and pdotdot and returns the braking index, n.
	"""
	
	n = 2 - (p0 * p2) / p1**2

	return n



def polar_cap_plasma_density(p0=1.0, p1=1.0e-12):
	"""
	Accepts the spin period and period derivative and returns the approximate polar cap plasma density in 
	the Goldreich-Julian model. Units: particles/cm^3.
	"""
	
	n = 7.0e10 * sqrt(1.0 / p0) * sqrt(p1 / 1.0e-15)

	return n



def polar_cap_radius(p0=1.0):
	"""
	Accepts the spin period and returns the approximate polar cap radius. Units: metres.
	"""
	
	r = 150 * sqrt(1 / p0)

	return r



def light_cyl_radius(p0=1.0):
	"""
	Accepts the spin period and returns the approximate light cylinder radius. Units: metres.
	"""

	rlc = c.value * p0 / (2 * pi)

	return rlc


def light_cyl_B(p0=1.0, p1=1.0e-12):
	"""
	Accepts spin period and first time derivative and returns the approximate magnetic 
	field strength at the light cylinder radius. Units: Gauss.
	"""

	blc = 9.2 * p0**(-2.5) * sqrt(p1 / 1.0e-15)

	return blc

	
	
if __name__ == '__main__':
	print "Testing dispersion delay (s) with freq1=100MHz, freq2=200MHz and DM=100pc/cm^3:",disp_delay(0.1, 0.2, 100.0)
	print "Testing E_rot (J), with P=0.0333s:",E_rot(0.0333)
	print "Testing E_loss (J/s), with P=0.0333s and Pdot=4.209e-13:",E_loss(0.0333, 4.209e-13)
	print "Testing B_surf (Gauss) with P=0.0333s and Pdot=4.209e-13:",min_mag_field(0.0333, 4.209e-13)
	print "Testing U_B with B=4e12 Gauss:",mag_energy_density(4e12)
	print "Testing char. age (years) with P=0.0333s and Pdot=4.209e-13:",char_age(0.0333, 4.209e-13)
	print "Testing frequency and period conversion..."
	print "f=29.946923 Hz, fdot=-3.77535E-10:",freq_period_convert(29.946923, -3.77535e-10)
	print "P=0.0333 s, Pdot=4.209e-13:",freq_period_convert(0.0333, 4.209e-13)
