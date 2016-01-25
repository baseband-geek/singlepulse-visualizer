#!/usr/bin/python






def disp_delay(freq1=0.184975, freq2=0.200335, DM=0):
	"""Accepts two frequencies as input (Low then High in GHz) as well as a DM and reports
	the delay in ms"""
	
	del_t = 4.148808 * DM * (freq1**(-2)-freq2**(-2))
	return del_t

	
	
if __name__ == '__main__':
	print "Testing with freq1=100MHz and freq2=200MHz with a DM of 100, dispersion delay is:"
	print disp_delay(0.1, 0.2, 100)