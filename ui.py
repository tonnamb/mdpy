import numpy as np
import matplotlib.pyplot as plt
import shutil
import os

def plot(end_more = 0, num_windows = None):
	pmf = np.genfromtxt('fe_ui.xy', comments='#')
	hist = np.genfromtxt('global_histogram.xy', comments='#')

	first_pmf = pmf[0, 1]
	last_pmf = pmf[-1, 1]

	for i, val in enumerate(pmf):
		if first_pmf != val[1]:
			start_i = i
			break

	for i, val in enumerate(pmf[::-1]):
		if last_pmf != val[1]:
			end_i = -i - end_more
			break

	pmf[:, 1] = pmf[:, 1] - max(pmf[(end_i-100):end_i, 1])

	fig = plt.figure()
	plt.plot(pmf[start_i:end_i, 0], pmf[start_i:end_i, 1])
	plt.xlabel(r'z ($\AA$)')
	plt.ylabel('PMF (kcal/mol)')
	fig.savefig('pmf.png', dpi=fig.dpi)
	plt.show()
	plt.close()

	fig = plt.figure()
	plt.plot(hist[:,0], hist[:,1])

	if type(num_windows) == type(50):
		histograms = []
		for i in range(1, num_windows+1):
			histograms.append(np.genfromtxt('histograms/%i' % i, comments='#'))
			plot_len = int(len(histograms[i-1])/2)
			plt.plot(histograms[i-1][:plot_len,0], histograms[i-1][:plot_len,1], ':')

	plt.xlabel(r'z ($\AA$)')
	plt.ylabel('Histogram Count')
	fig.savefig('histogram.png', dpi=fig.dpi)
	plt.show()
	plt.close()

	if type(num_windows) == type(50):
		return (pmf, hist, histograms)

def read_write_colvar(runrange, k_kcal, prevfolder):
	for i in runrange:
		with open(str(i) + "_window/out.colvars.traj") as file:
			line = file.readline() # skip first line
			print(i)
			if prevfolder:
				# skip one line since it duplicates with the last line from previously
				line = file.readline()
			line = file.readline()
			with open("ui/data/" + str(i),"a") as out:
				while line:
					rows = line.split()
					if rows[0] != '#':
						out.write(rows[1] + ' {0:.2f} '.format(k_kcal) + rows[2] + '\n')
					line = file.readline()

def extract_colvars(
	start, stop, k_eV, prevfolder=None,
	T=413, min_rc=0, max_rc=30, n=2000, ss=0):

	# example prevfolder = "2ns/k_1.7"
	if prevfolder:
		shutil.copytree("../../%s/ui/data" % prevfolder, "ui/data")
	else:
		os.makedirs("ui/data")
	
	eVtoKcal = 23.06
	k_kcal = k_eV*eVtoKcal

	runrange = range(start, stop+1)

	read_write_colvar(runrange, k_kcal, prevfolder)

	os.chdir("ui")
	os.system("ui.out -ui -T {0:d} -min {1:d} -max {2:d} -n {3:d} -u kcal -ss {4:d} -r -1 -v 2 > log.txt".format(T, min_rc, max_rc, n, ss))

def extract_2k(
	range1, k1, range2, k2, prevfolder=None,
	T=900, min_rc=0, max_rc=70, n=200, ss=1000):

	# example prevfolder = "2ns/k_1.7"
	if prevfolder:
		shutil.copytree("../../%s/ui/data" % prevfolder, "ui/data")
	else:
		os.makedirs("ui/data")

	eVtoKcal = 23.06
	k1_kcal = k1*eVtoKcal
	k2_kcal = k2*eVtoKcal

	read_write_colvar(range1, k1_kcal, prevfolder)
	read_write_colvar(range2, k2_kcal, prevfolder)

	os.chdir("ui")
	os.system("ui.out -ui -T {0:d} -min {1:d} -max {2:d} -n {3:d} -u kcal -ss {4:d} -r -1 -v 2 > log.txt".format(T, min_rc, max_rc, n, ss))
