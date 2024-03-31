import scipy
from scipy.integrate import quad

def g(e):
	return e**(12/19)/(1-e**2)*(1+121/304*e**2)**(870/2299)

def Fe0(e): 
	y, abserr = quad(lambda x: g(x)**4*(1-x**2)**(5/2)/x/(1+121/304*x**2), 0, e)
	y = 48/19/g(e)**4*y
	return y

def a_to_p(a, m1, m2):
	"""
	a, m1, m2 are floats
	[a]: AU
	[m1]: solar masses
	[m2]: solar masses
	"""
	mtot = m1+m2 # in solar masses
	p = (a**3/mtot)**(1/2) # in years
	p = p*365*24 # in hours
	return p
	
def check_hubbletime(a, e, m1, m2):
	p = a_to_p(a, m1, m2)
	time = 0.009829*p**(8/3)*(m1+m2)**(-2/3)*(m1*m2/(m1+m2))**(-1)*Fe0(e)
	if time <= 13.8:
		print('Under Hubble time')
	else:
		print('Over Hubble time')
	return time

print(check_hubbletime(1, 0.84, 65, 47))
print(check_hubbletime(0.1, 0.01, 65, 47))