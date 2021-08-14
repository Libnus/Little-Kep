# UNUSED IMPORT 
# from matplotlib import pyplot as plt

import math
import numpy as np
import pandas as pd


pi = np.pi
cos = np.cos
sin = np.sin
sqrt = np.sqrt

#calculate leap years
def leap_year(d):
	start_date = "01/01/2000"
	end_date = str(pd.to_datetime(start_date) + pd.DateOffset(days=d))
	end_year = int(end_date.split('-')[0])
	
	current_year = 2000
	num_years = 0
	while current_year <= end_year:
		if current_year % 4 == 0:
			if current_year % 100 == 0:
				if current_year % 400 == 0:
					num_years += 1
				else:
					num_years += 0
			else:
				num_years += 1
		else:
			num_years += 0

		current_year += 1

	print(num_years)
	return num_years

t = (12400 - (leap_year(12400))) * 86400 #time in seconds
# t = 1071358080 + 2.3383944e+11

# arrays for stars and planets
# array[(body name, a, e, i, W, w, mean_anomaly, orbital_period)] | for barycenters
body_attributes = np.array([['Star A', '0.950357892', '0.361691535', '1.02346528', '0', '72.0673282', '-50.1187383', '1.48328971'], ["Star B", '1.24774234', '0.361691535', '1.02346528', '0', '252.067328', '-50.1187383', '1.48328971'], ["Barycenter 1-1.1", '9.38327752', '0.066734531', '-1.46038708', '2.68534913', '190.129327', '220.440033', '10.802328']])
satellite_attributes = np.array([['Vanir', '2.3737279e-05', '0.0327855125', '0.824298426', '-174.45128', '246.813784', '88.5398428', '0.0197265468'], ['Moon1.1', '0.000558072', '0.0327855125', '0.824298426', '-174.45128', '66.8137841', '88.5398428', '0.0197265468']])

# a = 0.950357892 # semi-major axis in AU
# e = 0.361691535 # eccentricity 0 <= e <= 1
# i = 1.02346528 / (180/pi) # inclination
# W =  0 # longitude of ascending node
# w = 72.0673282 / (180/pi) # argument of pericenter in radians
# gc = 6.67408e-11 # gravitational constant
# Mean_epoch = -50.1187383 / (180/pi) # mean anomaly at epoch converted to radians


def calculate_mean_anomaly(M0, t, ot):
	
	n = (2*pi)/ot
	print(n)
	M_final = M0 + n * (t - 0)
	print(M_final)

	return M_final

def calculate_kepler_equation(e, M_epoch, ot):
	M_final = calculate_mean_anomaly(M_epoch, t, ot) 

	# calculate eccentric anomaly using Newton's algorithm
	error = 1.e-8

	if M_final < pi:
		E = M_final + e/2
	else:
		E = M_final - e/2

	ratio = 1
	while abs(ratio) > error:
		ratio = (E - e*sin(E) - M_final)/(1 - e*cos(E))
		E = E - ratio

	return E

# def calculate_true_anomaly(E, e):
# 	E = E * (180/pi)
# 	print("Eccentric ANOMALY (in degrees): ", E)
# 	v = (cos(E)-e)/(1-e*cos(E))
# 	v = math.acos(v)
# 	print("\nDISTANCE:", calculate_distance(a, e, v))
# 	return v

# def calculate_distance(a, e, v):
# 	r0 = (a*(1-e**2))/1+e*cos(v)
# 	return r0/v

def convert_2d(a, e, E):
	global coordinates
	P = a * (cos(E) - e)
	Q = a * sin(E) * sqrt(1 - e**2)

	distance = sqrt(P**2 + Q**2)

	print("\n2D Coordinates: ({}, {})\nDistance: {}".format(P, Q, distance))
	return P, Q
	# coordinates.append([P, Q])

def convert_3d(P, Q, w, i, W):
	x = cos(w) * P - sin(w) * Q 
	y = sin(w) * P + cos(w) * Q 

	z = sin(i) * y 
	y = cos(i) * y 

	xtemp = x
	x = cos(W) * xtemp - sin(W) * y 
	y = sin(W) * xtemp + cos(W) * y

	distance = sqrt(x**2 + y**2 + z**2) # all values here will be in au

	print("\n3D Coordinates: ({}, {}, {})\n".format(x, y, z))
	return x, y, z

# print("Eccentric Anomaly:", calculate_kepler_equation(e, Mean_epoch))

def main():
	for body in body_attributes:
		body_name = body[0]
		e = float(body[2])
		mean_anomaly_at_epoch = float(body[6]) / (180/pi)
		a = float(body[1])
		i = float(body[3])
		W = float(body[4])
		w = float(body[5])
		orbital_period = float(body[7]) * 3.154e+7

		print(body_name + ": ")
		E = calculate_kepler_equation(e, mean_anomaly_at_epoch, orbital_period)
		P, Q = convert_2d(a, e, E) # Eccentric anomaly in radians
		x, y, z = convert_3d(P, Q, w, i, W)

	for sattelite in satellite_attributes:
		sattelite_name = sattelite[0]
		e = float(sattelite[2])
		mean_anomaly_at_epoch = float(sattelite[6]) / (180/pi)
		a = float(sattelite[1])
		i = float(sattelite[3])
		W = float(sattelite[4])
		w = float(sattelite[5])
		orbital_period = float(sattelite[7]) * 3.154e+7

		print(sattelite_name + ": ")
		E = calculate_kepler_equation(e, mean_anomaly_at_epoch, orbital_period)
		P, Q = convert_2d(a, e, E)
		x, y, z = convert_3d(P, Q, w, i, W)


if __name__ == '__main__':
	main()

#code to plot orbit on graph, used for testing :p
# while t < 1.07136e+9:
# 	calculate_kepler_equation(e, Mean_epoch)
# 	t += 86400

# x, y = zip(*coordinates)

# plt.scatter(x, y)
# plt.show()