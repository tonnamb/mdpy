import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem

def read(fname, skip=0):

	data = []

	with open(fname) as file:
		count = 0
		for line in file:
			if count > skip: # skip header at count=0
				i = list(map(float, line.split()))
				data.append(i)
			count = count + 1

	return data

def read_str_array(fname, skip=0):

	if (type(fname) == type('')):
		head = header(fname)
		data = read(fname, skip)
	elif (type(fname) == type([])):
		head = header(fname[0])
		data = read(fname[0])
		for f in fname[1:]:
			next_data = read(f)
			if (data[-1][0] == next_data[0][0]): # remove overlapping frame
				next_data = next_data[1:]
			data = data + next_data
		data = data[skip:]
	return (data, head)

def header(fname):

	with open(fname) as file:
		h = file.readline()

	return h.strip()

def plot_graph(data, time_interval, name, unit):
	fig = plt.figure()
	avg = np.mean(data)
	time = np.linspace(0, len(data)*time_interval/1000000., len(data))
	plt.plot(time, data, '.', markersize=1)
	plt.axhline(y=avg, color='k')
	plt.xlabel('Time (ns)')
	plt.ylabel('%s (%s)' % (name, unit))
	plt.title('Mean = %.2f %s' % (avg, unit))
	fig.savefig('%s.png' % name, dpi=fig.dpi)
	plt.show()
	plt.close()

def plot(fname, time_interval=2000, skip_ns=0, poteng=False, temp=False, lx=False, lz=False):

	print("")

	skip_timestep = int(skip_ns*1000000/time_interval)
	data, head = read_str_array(fname, skip_timestep)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
		print("")
	else:
		print("WARNING: Header assumption is FALSE!")
		raise NameError('WARNING: Header assumption is FALSE!')

	data = np.array(data)

	time_interval = float(time_interval)

	if (skip_ns == 0):
		if (poteng):
			plot_graph(data[:, 1], time_interval, "potential_energy", "eV")
	
		if (temp):
			plot_graph(data[:, 4], time_interval, "temperature", "K")
	
		if (lx):
			plot_graph(data[:, 7], time_interval, "lx", "A")
	
		if (lz):
			plot_graph(data[:, 8], time_interval, "lz", "A")
	else:
		if (poteng):
			plot_graph(data[:, 1], time_interval, "potential_energy_skip%sns" % skip_ns, "eV")
	
		if (temp):
			plot_graph(data[:, 4], time_interval, "temperature_skip%sns" % skip_ns, "K")
	
		if (lx):
			plot_graph(data[:, 7], time_interval, "lx_skip%sns" % skip_ns, "A")
	
		if (lz):
			plot_graph(data[:, 8], time_interval, "lz_skip%sns" % skip_ns, "A")

	return data

def segment(fname, skip_ns=0, time_interval=2000, n_seg=20, conv_eV_to_kcal=False):

	skip_timestep = int(skip_ns*1000000/time_interval)
	data, head = read_str_array(fname, skip_timestep)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
		print("")

	data = np.array(data)

	if conv_eV_to_kcal:
		eV_to_kcal = 23.0605
		data[:, 1] = data[:, 1] * eV_to_kcal
		segment_plot(data[:, 1], skip_ns=skip_ns, time_interval=time_interval, n_seg=n_seg, unit='kcal/mol')
	else:
		segment_plot(data[:, 1], skip_ns=skip_ns, time_interval=time_interval, n_seg=n_seg)

	return data

def segment_plot(data, skip_ns=0, time_interval=2000, n_seg=20, name='', unit='eV'):

	seg_pe = np.array_split(data, n_seg)
	seg_len = np.mean(list(map(len, seg_pe)))*time_interval/1000000. # unit = ns
	seg_mean = list(map(np.mean, seg_pe))

	print("Skip: {0:.2f} ns".format(skip_ns))
	print("Production run: {0:.2f} ns".format(len(data)*2000./1000000.))
	    # 2000 = frequency of thermo output
	    # 1000000 = timesteps per ns
	avg_pe = np.mean(data)
	sem_pe = sem(data)
	print("Mean potential energy = {0:.2f} +/- {1:.2f} {2}".format(avg_pe, sem_pe, unit))
	print("Number of segments: {0}".format(n_seg))
	print("Segment length: {0:.2f} ns".format(seg_len))

	fig = plt.figure()
	plt.plot(seg_mean)
	plt.axhline(y=avg_pe, color='k', linestyle='--')
	ax = plt.gca()
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	plt.xlabel("Segment")
	plt.title("Potential energy (%s)" % unit)
	fig.savefig('poteng_segment_%s.png' % name, dpi=fig.dpi)
	plt.show()
	plt.close()


