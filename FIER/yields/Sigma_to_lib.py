#Sigma to FIER library
#Eric Matthews
#February 6, 2017

import sys 

file_in = sys.argv[1]
file_out = sys.argv[2]

file = open(file_in,'r')
lines = file.readlines()
file.close()

lines = lines[1:len(lines)]

file = open(file_out,'w')

TOT = 0.0

for line in lines:
	Z1 = line.split()[0]
	A1 = line.split()[1]
	i1 = line.split()[2]
	Y1 = float(line.split()[3])*100.0
	TOT += Y1
	dY1 = float(line.split()[4])*100.0
	try:
		Z2 = line.split()[5]
		A2 = line.split()[6]
		i2 = line.split()[7]
		Y2 = float(line.split()[8])*100.0
		TOT += Y2
		dY2 = float(line.split()[9])*100.0
	except IndexError as e:
		pass
	
	file.write(Z1 + ',' + A1 + ',' + i1 + ',' + "{:.6E}".format(Y1) + ',' + "{:.6E}".format(dY1) + '\n')
	file.write(Z2 + ',' + A2 + ',' + i2 + ',' + "{:.6E}".format(Y2) + ',' + "{:.6E}".format(dY2) + '\n')

file.close()
print(TOT)