#!/usr/bin/python

"""
A set of (potentially) handy functions for quickly calculating common 
physical properties of pulsars. These include:

    -- retrive pulsar info from PSRCAT:             get_pulsar
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


def get_pulsar(name):
	"""
	Accepts pulsar name (J or B name), accesses pulsar parameters from PSRCAT and calculates
	all physical properties available. Property value defaults to None if required parameter 
	from PSRCAT is missing.
	Returns a dictionary of pulsar parameters and calculated quantities.

	Looks for PSRCAT directory by accessing the environment variable 'PSRCAT_DIR'. 
	This should point to the directory where the PSRCAT database and binary executable are located.
	"""
	from subprocess import Popen,PIPE
	from os import getenv,listdir
	import sys

	print "Searching PSRCAT for pulsar {} ...".format(name)
	
	# find where PSRCAT resides
	psrcat_dir = getenv("PSRCAT_DIR")
	if psrcat_dir is None:
		print "***WARNING: Environment variable PSRCAT_DIR is not defined." 
		print "***NOTE: Please set PSRCAT_DIR to point at the directory containing psrcat.db "\
		    "and the binary executable."
		sys.exit("unable to locate PSRCAT from env. variable")

	elif 'psrcat' and 'psrcat.db' not in listdir(psrcat_dir):
		print "PSRCAT_DIR does not point to the location of psrcat.db and psrcat. Exitting."
		sys.exit("PSRCAT not in location of PSRCAT_DIR")

	else:
		psrcat = psrcat_dir + '/psrcat'
		psrcat_db = psrcat_dir + '/psrcat.db'


	# call PSRCAT for the pulsar name
	args = [psrcat, "-db_file", psrcat_db, "-o", "short_csv", "-nohead", "-nonumber",\
			"-c", "PSRJ PSRB RAJ DECJ RAJD DECJD GL GB DM F0 F1 F2 W50 W10 RM", name]

	psr = Popen(args, stdout=PIPE).communicate()[0]
	
	# convert output string into list:
	#        [PSRJ, PSRB, RAJ, DECJ, RAJD, DECJD, GL, GB, DM, F0, F1, F2, W50, W10, RM]
	psr = psr.strip().split(';')[:-1]

	# replace and "*" with None
	psr = [None if p is '*' else p for p in psr]

	psrj, psrb, ra, dec, ra_deg, dec_deg, gl, gb, dm, f0, f1, f2, w50, w10, rm = psr
	print "!!NOTE: Pulsar parameters will be None if not available from PSRCAT!!"
	print 50*"="
	print "Right ascension (J2000):                         {0}".format(ra)
	print "Declination (J2000):                             {0}".format(dec)
	print "Dispersion measure, DM (pc/cm^3):                {0}".format(dm)
	print "Spin frequency, f0 (Hz):                         {0}".format(f0)
	print "Spin frequency 1st time-derivative, f1 (Hz/s):   {0}".format(f1)
	print "Spin frequency 2nd time-derivative, f2 (Hz/s^2): {0}".format(f2)
	print "Pulse width at 50%, W50 (ms):                    {0}".format(w50)
	print "Pulse width at 10%, W10 (ms):                    {0}".format(w10)
	print "Rotation measure, RM (rad/m^2):                  {0}".format(rm)
	print '\n'+50*"="

	# output all calculable quantities:
	print "Producing all calculable quantities for the given pulsar: {}".format(name)
	print "If quantity is not calculable, will be displayed as None.\n"

	if dm:
		ddelay = disp_delay(0.185, 0.21572, float(dm))
		ddelay_inf = disp_delay_inf(0.185, float(dm))
	else:
		ddelay = None
		ddelay_inf = None

	
	if dm and rm:
		b_los = mean_parallel_B(float(rm), float(dm))
	else:
		b_los = None


	if f0:
		p0 = freq_period_convert(float(f0))
		erot = E_rot(p0)
		pcr = polar_cap_radius(p0)
		lcr = light_cyl_radius(p0)
	else:
		erot = None
		pcr = None
		lcr = None
		

	if f0 and f1:
		p0, p1 = freq_period_convert(float(f0), float(f1))
		eloss = E_loss(p0, p1)
		bmin = min_mag_field(p0, p1)
		mag_E_density = mag_energy_density(bmin)
		tau = char_age(p0, p1)
		pcpd = polar_cap_plasma_density(p0, p1)
		lcrB = light_cyl_B(p0, p1)
	else:
		eloss = None
		bmin = None
		uB = None
		tau = None
		pcpd = None
		lcrB = None


	if f0 and f1 and f2:
		p0, p1, p2 = freq_period_convert(float(f0), float(f1), float(f2))
		n = braking_index(p0, p1, p2)
	else:
		n = None


	print "Dispersion delay (ms) between 185 and 215.72 MHz:          {0}".format(ddelay)
	print "Dispersion delay (ms) from 185 MHz to infinite freq.:      {0}".format(ddelay_inf)
	print "Rotational energy (J):                                     {0}".format(erot)	
	print "Rotational energy loss rate (J/s):                         {0}".format(eloss)	
	print "Minimum surface magnetic field strength, Bmin (G):         {0}".format(bmin)
	print "Mean parallel magnetic field component (G) along the LOS:  {0}".format(b_los)
	print "Magnetic energy density (J/m^3), assuming B=Bmin:          {0}".format(mag_E_density)
	print "Characteristic age (yr):                                   {0}".format(tau)
	print "Braking index:                                             {0}".format(n)
	print "GJ polar cap plasma density (particles/cm^3):              {0}".format(pcpd)
	print "Polar cap radius (m):                                      {0}".format(pcr)
	print "Light cylinder radius (m):                                 {0}".format(lcr)
	print "Magnetic field strength at light cylinder radius (G):      {0}".format(lcrB)
	print 50*"="


	psr_dict = {'psrj':psrj, 'psrb':psrb, 'ra_str':ra, 'dec_str':dec, 'ra':ra_deg, 'dec':dec_deg,\
			    'l':gl, 'b':gb, 'dm':dm, 'f0':f0, 'f1':f1, 'f2':f2, 'w50':w50, 'w10':w10,\
			    'rm':rm, 'Erot':erot, 'Eloss':eloss, 'bmin':bmin, 'b_parallel':b_los,\
			    'u_B':mag_E_density, 'char_age':tau, 'braking_ind':n,\
			    'gj_density':pcpd, 'pc_radius':pcr, 'lcyl_radius':lcr, 'lcyl_B':lcrB}

	return psr_dict



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



def disp_delay_inf(freq=0.184975, DM=0.0):
	"""
	Accepts lower band-edge frequency (in GHz) and a DM (in pc/cm^3), 
	returns the dispersion delay in milliseconds realtive to infinite frequncy.
	Default: frequency = 184.975 MHz, DM = 0 pc/cm^3
	"""
	A = e.value**2 / (8 * pi**2 * m_e.value * c.value * eps0.value)
	A = A * (1.0/1.0e9)**2 * (1e6) * (1e3) *  pc.to('m').value

	del_t = A * DM * freq**(-2)

	return del_t



def E_rot(p0=1.0, *args):
	"""
	Accepts a spin period (in s) and returns the approximate rotational
	energy in J. (1 J = 10^7 erg)

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
	returns the approximate roational energy loss rate in J/s. (1 J/s = 10^7 erg/s)

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
	
	Conversion to Tesla is: 1 T = 10^4 G
	"""
		
	b_min = (3.2e19) * sqrt(p0 * p1)

	return b_min



def mean_parallel_B(rm=0.0, dm=100.0):
	"""
	Accepts a rotation measure (in rad/m^2) and dispersion measure (in pc/cm^3) and 
	returns the approximate mean value for the paralelle magnetic field compoent along 
	the line of sight, in Gauss.

	Conversion to Tesla is: 1 T = 10^4 G
	"""

	b_los = 1.23 * (rm / dm)

	return b_los / 1.0e6



def mag_energy_density(b=1e12):
	"""
	Accepts a magnetic field strength in Gauss and returns the theoretical
	energy density in J/m^3. 

	Conversion to erg/cm^3 is: 1 J/m^3 = 10 erg/cm^3
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
	get_pulsar('B0531+21')
