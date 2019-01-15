#Eric Matthews
#June 5, 2018
#Testing suite for FIER

import os

#Run the first FIER test deck
if( os.name != 'nt' ):
	os.system( 'rm -f testing/*.csv' )
	os.system( './fier.exe testing/testdeck.txt > /dev/null' )
	print('Running deck 1...')

	#Run the second FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( './fier.exe testing/testdeck2.txt > /dev/null' )
	print('Running deck 2...')

	#Run the third FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( './fier.exe testing/testdeck3.txt > /dev/null' )
	print('Running deck 3...')

	#Run the fourth FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( './fier.exe testing/testdeck4.txt > /dev/null' )
	print('Running deck 4...')
else:
	os.system( 'rm -f testing/*.csv' )
	os.system( 'fier.exe testing/testdeck.txt > nul' )
	print('Running deck 1...')

	#Run the second FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( 'fier.exe testing/testdeck2.txt > nul' )
	print('Running deck 2...')

	#Run the third FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( 'fier.exe testing/testdeck3.txt > nul' )
	print('Running deck 3...')

	#Run the fourth FIER test deck
	os.system( 'rm -f testing/*.csv' )
	os.system( 'fier.exe testing/testdeck4.txt > nul' )
	print('Running deck 4...')



#Test that the populations of the first test match the standard
file = open( 'testing/reference/ref_populations.csv', 'r' )
ref_pops = file.readlines()
file.close()
ref_pops_vals = []
ref_pops_vals.append( [ int(i) for i in ref_pops[0].split(',')[1:] ] )
ref_pops_vals.append( [ int(i) for i in ref_pops[1].split(',')[1:] ] )
ref_pops_vals.append( [ int(i) for i in ref_pops[2].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[3].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[4].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[5].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[6].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[7].split(',')[1:] ] )
ref_pops_vals.append( [ float(i) for i in ref_pops[8].split(',')[1:] ] )

file1 = open( 'testing/output/populations.csv', 'r' )
res_pops1 = file1.readlines()
file1.close()
res_pops1_vals = []
res_pops1_vals.append( [ int(i) for i in res_pops1[0].split(',')[1:] ] )
res_pops1_vals.append( [ int(i) for i in res_pops1[1].split(',')[1:] ] )
res_pops1_vals.append( [ int(i) for i in res_pops1[2].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[3].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[4].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[5].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[6].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[7].split(',')[1:] ] )
res_pops1_vals.append( [ float(i) for i in res_pops1[8].split(',')[1:] ] )


test1_pass = True
if( ref_pops_vals[0] != res_pops1_vals[0] ):
	test1_pass = False
if( ref_pops_vals[1] != res_pops1_vals[1] ):
	test1_pass = False
if( ref_pops_vals[2] != res_pops1_vals[2] ):
	test1_pass = False
if( ref_pops_vals[3] != res_pops1_vals[3] ):
	test1_pass = False
if( ref_pops_vals[4] != res_pops1_vals[4] ):
	test1_pass = False

diff = 0.0
for i in range( 0,len(ref_pops_vals[5]) ):
	diff += abs(ref_pops_vals[5][i] - res_pops1_vals[5][i])
if( sum( ref_pops_vals[5] ) != 0.0 ):
	diff = diff / sum( ref_pops_vals[5] )
if( diff > 1e-12 ):
	test1_pass = False

diff = 0.0
for i in range( 0,len(ref_pops_vals[6]) ):
	diff += abs(ref_pops_vals[6][i] - res_pops1_vals[6][i])
if( sum( ref_pops_vals[6] ) != 0.0 ):
	diff = diff / sum( ref_pops_vals[6] )
if( diff > 1e-12 ):
	test1_pass = False

diff = 0.0
for i in range( 0,len(ref_pops_vals[7]) ):
	diff += abs(ref_pops_vals[7][i] - res_pops1_vals[7][i])
if( sum( ref_pops_vals[7] ) != 0.0 ):
	diff = diff / sum( ref_pops_vals[7] )
if( diff > 1e-12 ):
	test1_pass = False

diff = 0.0
for i in range( 0,len(ref_pops_vals[8]) ):
	diff += abs(ref_pops_vals[8][i] - res_pops1_vals[8][i])
if( sum( ref_pops_vals[8] ) != 0.0 ):
	diff = diff / sum( ref_pops_vals[8] )
if( diff > 1e-12 ):
	test1_pass = False


if( test1_pass ):
	print( 'Passed: Test 1 matches reference populations.' )
else:
	raise Exception('Test 1 failed. Populations calculated do not match reference.')



#Test that the gammas output of the first test match the standard
file = open( 'testing/reference/ref_gamma_output.csv', 'r' )
ref_gammas = file.readlines()
file.close()
ref_gammas_vals = []
ref_gammas_vals.append( [ int(i) for i in ref_gammas[0].split(',')[2:] ] )
ref_gammas_vals.append( [ int(i) for i in ref_gammas[1].split(',')[2:] ] )
ref_gammas_vals.append( [ int(i) for i in ref_gammas[2].split(',')[2:] ] )
ref_gammas_vals.append( [ float(i) for i in ref_gammas[3].split(',')[2:] ] )
ref_gammas_vals.append( [ float(i) for i in ref_gammas[4].split(',')[2:] ] )
ref_gammas_vals.append( [ float(i) for i in ref_gammas[5].split(',')[2:] ] )
ref_gammas_vals.append( [ float(i) for i in ref_gammas[6].split(',')[2:] ] )

file1 = open( 'testing/output/gamma_output.csv', 'r' )
res_gammas1 = file1.readlines()
file1.close()
res_gammas1_vals = []
res_gammas1_vals.append( [ int(i) for i in res_gammas1[0].split(',')[2:] ] )
res_gammas1_vals.append( [ int(i) for i in res_gammas1[1].split(',')[2:] ] )
res_gammas1_vals.append( [ int(i) for i in res_gammas1[2].split(',')[2:] ] )
res_gammas1_vals.append( [ float(i) for i in res_gammas1[3].split(',')[2:] ] )
res_gammas1_vals.append( [ float(i) for i in res_gammas1[4].split(',')[2:] ] )
res_gammas1_vals.append( [ float(i) for i in res_gammas1[5].split(',')[2:] ] )
res_gammas1_vals.append( [ float(i) for i in res_gammas1[6].split(',')[2:] ] )


test1_pass = True
if( ref_gammas_vals[0] != res_gammas1_vals[0] ):
	test1_pass = False
if( ref_gammas_vals[1] != res_gammas1_vals[1] ):
	test1_pass = False
if( ref_gammas_vals[2] != res_gammas1_vals[2] ):
	test1_pass = False
if( ref_gammas_vals[3] != res_gammas1_vals[3] ):
	test1_pass = False
if( ref_gammas_vals[4] != res_gammas1_vals[4] ):
	test1_pass = False
if( ref_gammas_vals[5] != res_gammas1_vals[5] ):
	test1_pass = False

diff = 0.0
for i in range( 0,len(ref_gammas_vals[6]) ):
	diff += abs(ref_gammas_vals[6][i] - res_gammas1_vals[6][i])
if( sum( ref_gammas_vals[6] ) != 0.0 ):
	diff = diff / sum( ref_gammas_vals[6] )
if( diff > 1e-12 ):
	test1_pass = False


if( test1_pass ):
	print( 'Passed: Test 1 matches reference gamma output.' )
else:
	raise Exception('Test 1 failed. Gamma output calculated does not match reference.')



#Test that the population of the second test, which is strictly decay, are never negative
file2 = open( 'testing/output/populations2.csv', 'r' )
res_pops2 = file2.readlines()
file2.close()

res_pops2 = res_pops2[6] #We're only interested in the one time step
failed = False
for item in res_pops2.split(','):
	if( float(item) < 0.0 ):
		print(item)
		failed = True
if(failed):
	raise Exception('Test 2 failed. Negative populations from batch decay solution found.')
print('Passed: Test 2 produces no negative populations from the batch decay solution.')



#Test that the population of the third test, which is strictly continuous production, are never negative
file3 = open( 'testing/output/populations3.csv', 'r' )
res_pops3 = file3.readlines()
file3.close()
res_pops3 = res_pops3[6] #We're only interested in the one time step
failed = False
for item in res_pops3.split(','):
	if( float(item) < 0.0 ):
		print(item)
		failed = True
if(failed):
	raise Exception('Test 3 failed. Negative populations from continuous production solution found.')
print('Passed: Test 3 produces no negative populations from the continuous production solution.')
print('Note: negative populations can be produced to due bit resolution issues when the length of the time period is very short. The continuous production solution is more sensitive to this.')



#Test that the population of the fourth test, which is the same as the first test but with production broken into two neighboring and identical steps,match the first test
#This tests that the tracking of populations between time steps works
file4 = open( 'testing/output/populations4.csv', 'r' )
res_pops4 = file4.readlines()
file4.close()

res_pops4 = [ float(i) for i in res_pops4[7].split(',') ]
ref_pops_t = [ float(i) for i in ref_pops[6].split(',') ]
#Due to bit precision the results will not match exactly, however they should match to within reasonable error, much less than one part in a billion
diff_tot = 0.0
for i in range(1,len(res_pops4)):
	diff_tot += abs( ref_pops_t[i] - res_pops4[i] )
rel_err = diff_tot/sum(ref_pops_t)
if( rel_err < 1e-9 ):
	print( 'Passed: Test 4 matches reference populations to within numerical error. Gamma-ray emission tracking works.' )
	print( '   Numerical error: ' + str(rel_err) )
else:
	raise Exception('Test 4 failed. Populations calculated do not match reference to within numerical error.')

#And test the tracking of populations between time steps is correct and is reflected in the gammas output
file4 = open( 'testing/output/gamma_output4.csv', 'r' )
res_gammas4 = file4.readlines()
file4.close()

res_gammas4 = [ float(i) for i in res_gammas4[6].split(',') ]
ref_gammas_t = [ float(i) for i in ref_gammas[6].split(',') ]
#Due to bit precision the results will not match exactly, however they should match to within reasonable error, much less than one part in a billion
diff_tot = 0.0
for i in range(1,len(res_gammas4)):
	diff_tot += abs( ref_gammas_t[i] - res_gammas4[i] )
rel_err = diff_tot/sum(ref_gammas_t)
if( rel_err < 1e-9 ):
	print( 'Passed: Test 4 matches reference gammas to within numerical error. Gamma-ray emission tracking works.' )
	print( '   Numerical error: ' + str(rel_err) )
else:
	raise Exception('Test 4 failed. Gamma outputs do not match reference to within numerical error.')