"""
Test the new covariance matrix procedure.

"""

import numpy as np
import sys
import matplotlib.pyplot as plt

def proba(dtmodels, dtobs, covmat, n):

	dtmodels = np.matrix([dtmodels])
	dtobs = np.matrix([dtobs])
	covmat = np.matrix(covmat)

	diffdt = dtmodels-dtobs
	covmatdet = np.linalg.det(covmat)

	#print diffdt
	#print covmat.I
	#sys.exit()
	#blah = 1.0 / ((2*np.pi)**(n/2.0) * covmatdet**(1./2.0))
	#argexp = (-0.5 * (diffdt * covmat.I * diffdt.T)[0, 0])
	#print blah
	#print argexp
	#print int(argexp) ** float(blah)
	#print np.exp(-0.5 * 1e16 * (diffdt * covmat.I * diffdt.T)[0,0])
	#sys.exit()

	return (1./((2*np.pi)**(n/2.0) * covmatdet**(n/2.0)) * np.exp(-0.5*(diffdt * covmat.I * diffdt.T)[0,0]))

### Let's test this on a handmade three delay system
if 1:
	# We have AB, AC, BC

	# we measure two time shifts, B and C (assume A is 0)
	mean = (-8.8, -1.1)
	cov = [[0.4, 0.3],
		   [0.2, 0.4]]
	samples = []
	nsamples = 500
	delays = []

	for i in np.arange(nsamples):
		#shifts = np.random.multivariate_normal(mean, cov)
		shiftA = np.random.normal(loc=0, scale=0.5)
		shiftB = np.random.normal(loc=-8.8, scale=0.3)
		shiftC = np.random.normal(loc=-1.1, scale=0.6)
		shifts = (shiftA, shiftB, shiftC)

		samples.append(shifts)
		delays.append([0 - shifts[0], 0 - shifts[1], shifts[0] - shifts[1]])

	covsamples = np.cov(samples, rowvar=False)
	cov = np.cov(np.transpose(delays))
	# TODO: there is something I don't get here. If I use the covariance matrix of delays (cov) calculated from above, that matrix is singular, i.e. cannot be inverted in the proba function. (!?)

	# I use covsamples instead of cov in the following. That's not at all a matrix of delays, but for testing purposes it works...
	cov = cov

# This is the result I get from above
AB = -8.8
AC = -1.1
BC = 7.7

#ABvAB = cov[0][0]
#ABvAC = cov[0][1]
#ABvBC = cov[0][2]
#ACvAC = cov[1][1]
#ACvBC = cov[1][2]
#BCvBC = cov[2][2]

#print ABvAB, ABvAC, ABvBC, ACvAC, ACvBC, BCvBC

ABvAB = 0.74
ABvAC = 0.43
ABvBC = -0.08
ACvAC = 0.75
ACvBC = 0.08
BCvBC = 0.39


dtobs = [AB, AC, BC]
covmat = [
	[ABvAB, ABvAC, ABvBC],
	[ABvAC, ACvAC, ACvBC],
	[ABvBC, ACvBC, BCvBC]
]

# AB, AC
covmatreducA = [
	[ABvAB, ABvAC],
	[ABvAC, ACvAC]
]
dtobsreducA = [AB, AC]

# BA, BC
covmatreducB = [
	[ABvAB, ABvBC],
	[ABvBC, BCvBC]
]
dtobsreducB = [-AB, BC]

# CA, CB
covmatreducC = [
	[ACvAC, ACvBC],
	[ACvBC, BCvBC]
]
dtobsreducC = [-AC, -BC]

### AB
plt.figure()
dtAB = np.linspace(AB-3, AB+3, 100)
probAB = []
probABreduc = []
for dt in dtAB:
	dtmodels = [dt, AC, BC]
	probAB.append(proba(dtmodels, dtobs, covmat, n=3))
	dtmodelsreduc = [dt, AC]
	probABreduc.append(proba(dtmodelsreduc, dtobsreducA, covmatreducA, n=2))

plt.plot(dtAB, probAB/max(probAB), 'b')
plt.plot(dtAB, probABreduc/max(probABreduc), '--r')


dtBA = np.linspace(-AB-3, -AB+3, 100)
probBAreduc = []
for dt in dtBA:
	dtmodelsreduc = [dt, BC]
	probBAreduc.append(proba(dtmodelsreduc, dtobsreducB, covmatreducA, n=2))

plt.plot(dtAB, probBAreduc/max(probBAreduc), '--k')
plt.show()


### AC
plt.figure()
dtAC = np.linspace(AC-3, AC+3, 100)
probAC = []
probACreduc = []
for dt in dtAC:
	dtmodels = [AB, dt, BC]
	probAC.append(proba(dtmodels, dtobs, covmat, n=3))
	dtmodelsreduc = [-dt, -BC]
	probACreduc.append(proba(dtmodelsreduc, dtobsreducC, covmatreducC, n=2))

plt.plot(dtAC, probAC/max(probAC), 'b')
#plt.plot(dtAC, probACreduc/max(probACreduc), '--r')
plt.show()


### BC
plt.figure()
dtBC = np.linspace(BC-3, BC+3, 100)
probBC = []
probBCreduc = []
for dt in dtBC:
	dtmodels = [AB, AC, dt]
	probBC.append(proba(dtmodels, dtobs, covmat, n=3))
	dtmodelsreduc = [-AB, dt]
	probBCreduc.append(proba(dtmodelsreduc, dtobsreducB, covmatreducB, n=2))

plt.plot(dtBC, probBC/max(probBC), 'b')
plt.plot(dtBC, probBCreduc/max(probBCreduc), '--r')
plt.show()




sys.exit()
### let's define the obs delays and covcoeffs

# obs delays
AB = -8.8
AC = -1.1
AD = -13.8
BC = 7.7
BD = -5.1
CD = -12.7

# mean cov
ABvAB = 0.74
ABvAC = 0.43
ABvAD = 0.43
ABvBC = -0.08
ABvBD = -0.08
ABvCD = 0.00

ACvAC = 0.75
ACvAD = 0.44
ACvBC = 0.08
ACvBD = -0.01
ACvCD = -0.09

ADvAD = 1.0
ADvBC = -0.01
ADvBD = 0.40
ADvCD = 0.41

BCvBC = 0.39
BCvBD = 0.07
BCvCD = -0.09

BDvBD = 0.73
BDvCD = 0.41

CDvCD = 0.74



# full set
dtobs = [AB, AC, AD, BC, BD, CD]
covmat = [
		  [ABvAB, ABvAC, ABvAD, ABvBC, ABvBD, ABvCD],
		  [ABvAC, ACvAC, ACvAD, ACvBC, ACvBD, ACvCD],
		  [ABvAD, ACvAD, ADvAD, ADvBC, ADvBD, ADvCD],
		  [ABvBC, ACvBC, ADvBC, BCvBC, BCvBD, BCvCD],
		  [ABvBD, ACvBD, ADvBD, BCvBD, BDvBD, BDvCD],
		  [ABvCD, ACvCD, ADvCD, BCvCD, BDvCD, CDvCD]
		  ]

# AB, AC, AD
covmatreducA = [
	[ABvAB, ABvAC, ABvAD],
	[ABvAC, ACvAC, ACvAD],
	[ABvAD, ACvAD, ADvAD]
]
dtobsreducA = [AB, AC, AD]


# BA, BC, BD
covmatreducB = [
	[ABvAB, ABvBC, ABvBD],
	[ABvBC, BCvBC, BCvBD],
	[ABvBD, BCvBD, BDvBD]
]
dtobsreducB = [-AB, BC, BD]


# CA, CB, CD
covmatreducC = [
	[ACvAC, ACvBC, ACvCD],
	[ACvBC, BCvBC, BCvCD],
	[ACvCD, BCvCD, CDvCD]
]
dtobsreducC = [-AC, -BC, CD]


# DA, DB, DC
covmatreducD = [
	[ADvAD, ADvBD, ADvCD],
	[ADvBD, BDvBD, BDvCD],
	[ADvCD, BDvCD, CDvCD]
]
dtobsreducD = [-AD, -BD, -CD]




### AB
plt.figure()
dtAB = np.linspace(AB-3, AB+3, 100)
probAB = []
probABreduc = []
for dt in dtAB:
	dtmodels = [dt, AC, AD, BC, BD, CD]
	probAB.append(proba(dtmodels, dtobs, covmat, n=6))
	dtmodelsreduc = [dt, AC, AD]
	probABreduc.append(proba(dtmodelsreduc, dtobsreducA, covmatreducA, n=3))

dtBA = np.linspace(-AB-3, -AB+3, 100)
probBAreduc = []
for dt in dtBA:
	dtmodelsreduc = [dt, BC, BD]
	probBAreduc.append(proba(dtmodelsreduc, dtobsreducB, covmatreducB, n=3))

plt.plot(dtAB, probAB/max(probAB), 'b')
plt.plot(dtAB, probABreduc/max(probABreduc), 'r')
plt.plot(-dtBA, probBAreduc/max(probBAreduc), 'g')


### BC
plt.figure()
dtBC = np.linspace(BC-3, BC+3, 100)
probBC = []
probBCreduc = []
for dt in dtBC:
	dtmodels = [AB, AC, AD, dt, BD, CD]
	probBC.append(proba(dtmodels, dtobs, covmat, n=6))
	dtmodelsreduc = [-AB, dt, BD]
	probBCreduc.append(proba(dtmodelsreduc, dtobsreducB, covmatreducB, n=3))

dtCB = np.linspace(-BC-3, -BC+3, 100)
probCBreduc = []
for dt in dtCB:
	dtmodelsreduc = [-AC, dt, CD]
	probCBreduc.append(proba(dtmodelsreduc, dtobsreducC, covmatreducC, n=3))


plt.plot(dtBC, probBC/max(probBC), 'b')
plt.plot(dtBC, probBCreduc/max(probBCreduc), 'r')
plt.plot(-dtCB, probCBreduc/max(probCBreduc), 'g')


### CD
plt.figure()
dtCD = np.linspace(CD-3, CD+3, 100)
probCD = []
probCDreduc = []
for dt in dtCD:
	dtmodels = [AB, AC, AD, BC, BD, dt]
	probCD.append(proba(dtmodels, dtobs, covmat, n=6))
	dtmodelsreduc = [-AC, -BC, dt]
	probCDreduc.append(proba(dtmodelsreduc, dtobsreducC, covmatreducC, n=3))

dtDC = np.linspace(-CD-3, -CD+3, 100)
probDCreduc = []
for dt in dtDC:
	dtmodelsreduc = [-AD, -BD, dt]
	probDCreduc.append(proba(dtmodelsreduc, dtobsreducD, covmatreducD, n=3))

plt.plot(dtCD, probCD/max(probCD), 'b')
plt.plot(dtCD, probCDreduc/max(probCDreduc), 'r')
plt.plot(-dtDC, probCDreduc/max(probCDreduc), 'g')


### AD
plt.figure()
dtAD = np.linspace(AD-3, AD+3, 100)
probAD = []
probADreduc = []
for dt in dtAD:
	dtmodels = [AB, AC, dt, BC, BD, CD]
	probAD.append(proba(dtmodels, dtobs, covmat, n=6))
	dtmodelsreduc = [AB, AC, dt]
	probADreduc.append(proba(dtmodelsreduc, dtobsreducA, covmatreducA, n=3))

dtDA = np.linspace(-AD-3, -AD+3, 100)
probDAreduc = []
for dt in dtDA:
	dtmodelsreduc = [dt, -BD, -CD]
	probDAreduc.append(proba(dtmodelsreduc, dtobsreducD, covmatreducD, n=3))

plt.plot(dtAD, probAD/max(probAD), 'b')
plt.plot(dtAD, probADreduc/max(probADreduc), 'r')
plt.plot(-dtDA, probDAreduc/max(probDAreduc), 'g')


plt.show()
