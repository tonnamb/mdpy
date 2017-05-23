import mdtraj
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import chemfiles
import time
from mdpy.entropy import pbc_wrap
from scipy.stats import sem
import pdb

def avg_eg_peak_width(f_traj, eg_types=[7, 8, 9, 10], num_bins=5000):

	atoms_in_eg = 10.0

	if type(f_traj) == type([]):
		traj = chemfiles.Trajectory(f_traj[0])
	elif type(f_traj) == type(''):
		traj = chemfiles.Trajectory(f_traj)

	n_frame = traj.nsteps()
	frame = traj.read()
	eg_atoms = chemfiles.Selection("atoms: name {0:d} or name {1:d} or name {2:d} or name {3:d}".format(*eg_types)).evaluate(frame)

	box_x, box_y, box_z = frame.cell().lengths()
	box_area = box_x*box_y
	eg_z = []

	positions = frame.positions()
	eg_z.append( list(positions[eg_atoms, 2]) )
	atoms_per_frame = len(eg_z[0])

	for i in range(traj.nsteps()-1):
		frame = traj.read()
		positions = frame.positions()
		eg_z.append( list(positions[eg_atoms, 2]) )

	if type(f_traj) == type([]) and len(f_traj) > 1:
		for ft in f_traj[1:]:
			traj = chemfiles.Trajectory(ft)
			n_frame += traj.nsteps()
			for i in range(traj.nsteps()-1):
				frame = traj.read()
				positions = frame.positions()
				eg_z.append( list(positions[eg_atoms, 2]) )

	eg_z = [pbc_wrap(x, box_z) for li in eg_z for x in li]

	count, binEdges = np.histogram(eg_z, num_bins, density=False)
	dist_bin = binEdges[1]-binEdges[0] # Unit = A
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

	density = count/(box_area*dist_bin*n_frame*atoms_in_eg) # Convert to density [=] EG / A^3

	for index, val in enumerate(density):
		print('{0:d} {1:.5f}'.format(index, val))

	fig = plt.figure()
	plt.plot(bincenters,density, 'r-')
	plt.xlabel('z (nm)')
	plt.ylabel(r'$\rho_{EG}$ ($EG/A^3$)')
	fig.savefig('eg_density.png', dpi=fig.dpi)
	plt.show()
	plt.close()

	print('Bin with density=0.0001 contribution: %.3f' % (dist_bin*box_area*0.0001))
	bin_to_cover = 4.5/dist_bin
	print('Number of bin to cover 4.5 A: %.1f' % bin_to_cover)

	bottom_2 = int(input('Enter index of bottom peak end (> 0.00010): '))
	print('Recommended bottom peak start: %.0d' % (bottom_2 - round(bin_to_cover)))
	bottom_1 = int(input('Enter index of bottom peak start: '))
	top_1 = int(input('Enter index of top peak start (> 0.00010): '))
	print('Recommended top peak end: %.0d' % (top_1 + round(bin_to_cover)))
	top_2 = int(input('Enter index of top peak end: '))

	width_bottom = (bottom_2-bottom_1)*dist_bin
	width_top = (top_2-top_1)*dist_bin
	print('Width of bottom peak = %.2f A' % width_bottom)
	print('Width of top peak = %.2f A' % width_top)

	with open('eg_peak_width.txt', 'a') as f:
		f.write(time.strftime("%c"))
		f.write('\n')
		f.write('Width of bottom peak = %.2f A\n' % width_bottom)
		f.write('Width of top peak = %.2f A\n' % width_top)

def count_eg_displaced(
	f_traj, eg_types=[7, 8, 9, 10], num_bins=5000, n_seg=5):

	atoms_in_eg = 10.0

	if type(f_traj) == type([]):
		traj = chemfiles.Trajectory(f_traj[0])
	elif type(f_traj) == type(''):
		traj = chemfiles.Trajectory(f_traj)

	n_frame = traj.nsteps()
	frame = traj.read()
	eg_atoms = chemfiles.Selection("atoms: name {0:d} or name {1:d} or name {2:d} or name {3:d}".format(*eg_types)).evaluate(frame)

	box_x, box_y, box_z = frame.cell().lengths()
	box_area = box_x*box_y
	eg_z = []

	positions = frame.positions()
	eg_z.append( list(positions[eg_atoms, 2]) )
	atoms_per_frame = len(eg_z[0])

	for i in range(traj.nsteps()-1):
		frame = traj.read()
		positions = frame.positions()
		eg_z.append( list(positions[eg_atoms, 2]) )

	if type(f_traj) == type([]) and len(f_traj) > 1:
		for ft in f_traj[1:]:
			traj = chemfiles.Trajectory(ft)
			n_frame += traj.nsteps()
			for i in range(traj.nsteps()-1):
				frame = traj.read()
				positions = frame.positions()
				eg_z.append( list(positions[eg_atoms, 2]) )

	eg_z = [pbc_wrap(x, box_z) for li in eg_z for x in li]

	seg_egz = np.array_split(eg_z, n_seg)
	seg_len = np.mean(list(map(len, seg_egz)))

	displaced_array = []

	with open('displaced_eg.txt', 'a') as f:
		f.write(time.strftime("%c"))
		f.write('\nNumber of segments: %d\n' % n_seg)

	for i in range(n_seg):
		count, binEdges = np.histogram(seg_egz[i], num_bins, density=False)
		dist_bin = binEdges[1]-binEdges[0] # Unit = A
		bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
		n_frame_seg = len(seg_egz[i])/atoms_per_frame

		density = count/(box_area*dist_bin*n_frame_seg*atoms_in_eg) # Convert to density [=] EG / A^3

		for index, val in enumerate(density):
			print('{0:d} {1:.5f}'.format(index, val))

		fig = plt.figure()
		plt.plot(bincenters,density, 'r-')
		plt.xlabel('z (nm)')
		plt.ylabel(r'$\rho_{EG}$ ($EG/A^3$)')
		fig.savefig('eg_density_%d.png' % (i+1), dpi=fig.dpi)
		plt.show()
		plt.close()

		print('Bin with density=0.0002 contribution: %.3f' % (dist_bin*box_area*0.0002))
		bin_to_cover = 4.5/dist_bin
		print('Number of bin to cover 4.5 A: %.1f' % bin_to_cover)

		bottom_2 = int(input('Enter index of bottom peak end (> 0.00020): '))
		print('Recommended bottom peak start: %.0d' % (bottom_2 - round(bin_to_cover)))
		# bottom_1 = int(input('Enter index of bottom peak start: '))
		bottom_1 = int(bottom_2 - round(bin_to_cover))
		top_1 = int(input('Enter index of top peak start (> 0.00020): '))
		print('Recommended top peak end: %.0d' % (top_1 + round(bin_to_cover)))
		# top_2 = int(input('Enter index of top peak end: '))
		top_2 = int(top_1 + round(bin_to_cover))

		displaced_array.append(get_displaced(density, box_area, dist_bin, f_traj, bottom_1, bottom_2, top_1, top_2, i, n_frame_seg))

	d = {}
	d['displaced_array'] = displaced_array
	d['displaced_mean'] = np.mean(displaced_array)
	d['displaced_sem'] = sem(displaced_array)

	print('Displaced array:')
	print('%.2f '*len(d['displaced_array']) % tuple(d['displaced_array']))
	print('Displaced: %.2f +/- %.2f' % (d['displaced_mean'], d['displaced_sem']))

	with open('displaced_eg.txt', 'a') as f:
		f.write('\nSummary\n')
		f.write('Displaced array: ')
		f.write('%.2f '*len(d['displaced_array']) % tuple(d['displaced_array']))
		f.write('\n')
		f.write('Displaced: %.2f +/- %.2f\n\n' % (d['displaced_mean'], d['displaced_sem']))

	# d['density'] = density
	# d['box_area'] = box_area
	# d['dist_bin'] = dist_bin
	# d['f_traj'] = f_traj
	# d['bottom_1'] = bottom_1
	# d['bottom_2'] = bottom_2
	# d['top_1'] = top_1
	# d['top_2'] = top_2

	return d

def get_displaced(density, box_area, dist_bin, f_traj,
	bottom_1, bottom_2, top_1, top_2, seg_i, n_frame_seg):

	total_count = np.sum(density)*box_area*dist_bin
	bottom_count = np.sum(density[bottom_1:bottom_2])*box_area*dist_bin
	top_count = np.sum(density[top_1:top_2])*box_area*dist_bin
	displaced = top_count - bottom_count

	print(seg_i+1)
	print(f_traj)
	print('Number of frames in segment: %d' % n_frame_seg)
	print('Bin range bottom: %d' % (bottom_2-bottom_1))
	print('Bin range top: %d' % (top_2-top_1))
	print('Total count: %.2f' % total_count)
	print('Top count: %.2f' % top_count)
	print('Bottom count: %.2f' % bottom_count)
	print('Displaced: %.2f' % displaced)

	with open('displaced_eg.txt', 'a') as f:
		f.write('\n%d\n' % (seg_i+1))
		f.write('Number of frames in segment: %d\n' % n_frame_seg)
		f.write('Bin with density=0.00020 contribution: %.3f\n' % (dist_bin*box_area*0.0002))
		f.write('Number of bin to cover 4.5 A: %.1f\n' % (4.5/dist_bin))
		if type(f_traj) == type([]):
			f.write(', '.join(f_traj))
		elif type(f_traj) == type(''):
			f.write(f_traj)
		f.write('\n')
		f.write('Bottom peak 1: %d\n' % bottom_1)
		f.write('Bottom peak 2: %d\n' % bottom_2)
		f.write('Top peak 1: %d\n' % top_1)
		f.write('Top peak 2: %d\n' % top_2)
		f.write('Bin range bottom: %d\n' % (bottom_2-bottom_1))
		f.write('Bin range top: %d\n' % (top_2-top_1))
		f.write('Total count: %.2f\n' % total_count)
		f.write('Top count: %.2f\n' % top_count)
		f.write('Bottom count: %.2f\n' % bottom_count)
		f.write('EG Displaced: %.2f\n' % displaced)
	
	return displaced

def get_density(F_TRAJ, ATOM_TYPE, SURF_ATOMS, NUMBINS=100, ATOM_NAME='pvp', ON_TOP=False):

	if type(F_TRAJ) == type([]):
		traj = chemfiles.Trajectory(F_TRAJ[0])
	elif type(F_TRAJ) == type(''):
		traj = chemfiles.Trajectory(F_TRAJ)

	n_frame = traj.nsteps()
	frame = traj.read()
	ATOMS = chemfiles.Selection("atoms: name {0:d}".format(ATOM_TYPE)).evaluate(frame)
	ADJ_SURF_ATOMS = [ x - 1 for x in SURF_ATOMS ]

	BOX_X, BOX_Y, BOX_Z = frame.cell().lengths()
	BOX_AREA = BOX_X * BOX_Y

	atom_z = []
	surf_z = []
	positions = frame.positions()
	atom_z.append( list(positions[ATOMS, 2]) )
	surf_z.append( list(positions[ADJ_SURF_ATOMS, 2]) )
	ATOMS_PER_FRAME = len(atom_z[0])

	for i in range(traj.nsteps()-1):
		frame = traj.read()
		positions = frame.positions()
		atom_z.append( list(positions[ATOMS, 2]) )
		surf_z.append( list(positions[ADJ_SURF_ATOMS, 2]) )

	if type(F_TRAJ) == type([]) and len(F_TRAJ) > 1:
		for ft in F_TRAJ[1:]:
			traj = chemfiles.Trajectory(ft)
			n_frame += traj.nsteps()
			for i in range(traj.nsteps()-1):
				frame = traj.read()
				positions = frame.positions()
				atom_z.append( list(positions[ATOMS, 2]) )
				surf_z.append( list(positions[ADJ_SURF_ATOMS, 2]) )

	atom_z = [ [pbc_wrap(x, BOX_Z) for x in frame] for frame in atom_z ]
	surf_z = [ [pbc_wrap(x, BOX_Z) for x in frame] for frame in surf_z ]
	avg_surf_z = [ np.mean(frame) for frame in surf_z ]

	fig = plt.figure()
	plt.hist([val for frame in surf_z for val in frame], 100)
	plt.xlabel('z of surface atoms (A)')
	plt.ylabel('counts')
	plt.title(ATOM_NAME)
	fig.savefig('ag_%s.png' % ATOM_NAME.replace(' ','_'), dpi=fig.dpi)
	plt.show()
	plt.close()

	adj_atom_z = []
	for atom_frame, avg_surf in zip(atom_z, avg_surf_z):
		if ON_TOP:
			adj_atom_z.extend([atom - avg_surf for atom in atom_frame])
		else:
			adj_atom_z.extend([avg_surf - atom for atom in atom_frame])

	count, binEdges = np.histogram(adj_atom_z, NUMBINS, density=False)
	dist_bin = binEdges[1]-binEdges[0] # Unit = A
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

	density = count/(BOX_AREA*dist_bin*n_frame) # Convert to density [=] atom / A^3

	data = {}
	data['name'] = ATOM_NAME
	data['density'] = density
	data['bincenters'] = bincenters

	return data

def plot_density(data_array, plot_name, plot_title='', plot_range=[0, 15]):

	fig = plt.figure()
	for data in data_array:
		plt.plot(data['bincenters'], data['density'], label=data['name'])
	plt.legend()
	plt.xlabel('z (nm)')
	plt.ylabel(r'$\rho$ (atom / $A^3$)')
	plt.xlim(plot_range)
	plt.title(plot_title)
	fig.savefig('density_%s.png' % plot_name, dpi=fig.dpi)
	plt.show()
	plt.close()

def calc_conc(read, pdb, ag_type, num_pvp, name):
	# Load trajectory
	t = mdtraj.load(read, top=pdb)
	unitcell_x = t.unitcell_lengths[-1][0]
	unitcell_y = t.unitcell_lengths[-1][1]
	unitcell_z = t.unitcell_lengths[-1][2]
	topology = t.topology

	# Select Ag atoms
	agatom = [atom.index for atom in topology.atoms if (atom.name == str(ag_type))]
	agpos_z = t.xyz[-1, agatom, 2].flatten()

	# Assumes that Ag layer is not in between periodic boundary
	agslab_thickness = np.max(agpos_z) - np.min(agpos_z) # nm

	# Effective z for liquid occupation
	free_z = unitcell_z - agslab_thickness # nm

	# Liquid volume
	liq_vol = free_z * unitcell_x * unitcell_y # nm^3

	# Avogadro's constant
	avo = 6.022e23 # num/mol

	# Concentration of PVP
	pvp_conc = (num_pvp/avo) / (liq_vol*1e-24) # mol/L

	# Output to screen
	print('Lx = {0:.2f} nm'.format(unitcell_x))
	print('Ly = {0:.2f} nm'.format(unitcell_y))
	print('Lz = {0:.2f} nm'.format(unitcell_z))
	print('Bounds of Ag atoms = [{0:.2f}, {1:.2f}]'.format(np.min(agpos_z), np.max(agpos_z)))
	print('Thickness of Ag slab = {0:.2f} nm'.format(agslab_thickness))
	print('Effective z for liquid occupation = {0:.2f} nm'.format(free_z))
	print('Liquid volume = {0:.1f} nm^3'.format(liq_vol))
	print('Number of {0} = {1}'.format(name, num_pvp))
	print('Concentration of {0} = {1:.4f} mol/L = {2:.1f} mM'.format(name, pvp_conc, (pvp_conc * 1000)))

	with open('calc_conc_%s.txt' % name, 'a') as f:
		f.write('\nLx = {0:.2f} nm\n'.format(unitcell_x))
		f.write('Ly = {0:.2f} nm\n'.format(unitcell_y))
		f.write('Lz = {0:.2f} nm\n'.format(unitcell_z))
		f.write('Bounds of Ag atoms = [{0:.2f}, {1:.2f}]\n'.format(np.min(agpos_z), np.max(agpos_z)))
		f.write('Thickness of Ag slab = {0:.2f} nm\n'.format(agslab_thickness))
		f.write('Effective z for liquid occupation = {0:.2f} nm\n'.format(free_z))
		f.write('Liquid volume = {0:.1f} nm^3\n'.format(liq_vol))
		f.write('Number of {0} = {1}\n'.format(name, num_pvp))
		f.write('Concentration of {0} = {1:.4f} mol/L = {2:.1f} mM\n'.format(name, pvp_conc, (pvp_conc * 1000)))


def plot_eg_density(traj_f, pdb_f, atom_type):

	# Load trajectory file
	t = mdtraj.load(traj_f, top=pdb_f)
	topology = t.topology

	# Get unitcell lengths and number of frames
	unitcell_x = t.unitcell_lengths[0,0] # nm
	unitcell_y = t.unitcell_lengths[0,1] # nm
	unitcell_area = unitcell_x * unitcell_y # nm^2
	n_frame = t.n_frames

	# Get atom index
	str_atom_type = list(map(str, atom_type))
	egatom = [atom.index for atom in topology.atoms if (atom.name in str_atom_type)]

	# Get atom position
	opos = t.xyz[:, egatom]
	# All frames
	# Returns numpy array with dimensions [timeframe, atom, xyz]
	z = opos[:, :, 2].flatten()

	num_bins = 200
	count, binEdges = np.histogram(np.ndarray.flatten(z), num_bins, density=False)
	dist_bin = binEdges[1]-binEdges[0] # Unit = nm
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

	atoms_in_eg = 10.0
	av = 6.022e23 # #/mol
	mw = 62.07 # g/mol
	cm3_nm3 = 1e-21 # cm^3/nm^3

	density = (count*mw/av)/(unitcell_area*dist_bin*n_frame*atoms_in_eg*cm3_nm3) # Convert to density

	fig = plt.figure()
	plt.plot(bincenters,density, 'r-')
	plt.xlabel('z (nm)')
	plt.ylabel(r'$\rho_{EG}$ ($g/cm^3$)')
	fig.savefig('eg_density.png', dpi=fig.dpi)
	plt.show()
	plt.close()

	return (bincenters, density)

def setup_dict(path, traj_f, pdb_f, plot_type, surf_id, name='', bottom=True):

	dictObj = {}
	dictObj['name'] = name
	dictObj['path'] = path

	if type(traj_f) == type([]):
		dictObj['traj_f'] = list(map(lambda x: path+x, traj_f))
	elif type(traj_f) == type(""):
		dictObj['traj_f'] = path+traj_f

	dictObj['pdb_f'] = path+pdb_f
	dictObj['traj'] = mdtraj.load(dictObj['traj_f'], discard_overlapping_frames=True, top=dictObj['pdb_f'])
	dictObj['topo'] = dictObj['traj'].topology

	dictObj['x'] = dictObj['traj'].unitcell_lengths[0,0] # nm
	dictObj['y'] = dictObj['traj'].unitcell_lengths[0,1] # nm
	dictObj['area'] = dictObj['x'] * dictObj['y']
	# print("Area = %.2f nm^2" % dictObj['area'])
	dictObj['n_frames'] = dictObj['traj'].n_frames

	# Get atom index
	if type(plot_type) == type([]):
		dictObj['str_plot_type'] = list(map(str, plot_type))
	elif type(plot_type) == type(1):
		dictObj['str_plot_type'] = str(plot_type)
	
	dictObj['plot_atom'] = [atom.index for atom in dictObj['topo'].atoms if (atom.name in dictObj['str_plot_type'])]
	
	if type(surf_id) == type(range(1)):
		surf_atom = list(surf_id)
	elif type(surf_id) == type([]):
		surf_atom = surf_id
	elif type(surf_id) == type((1, 2)):
		surf_atom = [ i for j in surf_id for i in j ]

	dictObj['surf_atom'] = list(map(lambda x: x-1, surf_atom))

	# Get surface position
	dictObj['surf_z'] = dictObj['traj'].xyz[:, dictObj['surf_atom'], 2]
	dictObj['surf_avg_z'] = np.mean(dictObj['surf_z'].flatten())

	# Check if surface index is correct by histogram
	fig = plt.figure()
	plt.hist(dictObj['surf_z'].flatten(), 50)
	plt.axvline(x=dictObj['surf_avg_z'], linewidth=2, color='r')
	plt.xlabel('z (nm)')
	plt.ylabel('Histogram count')
	plt.title('Z-position of surface layer for %s' % name)
	fig.savefig('hist_surf_%s.png' % name, dpi=fig.dpi)
	plt.show()
	plt.close()

	# Get plot atom position
	dictObj['plot_z'] = dictObj['traj'].xyz[:, dictObj['plot_atom'], 2]
	dictObj['plot_xyz'] = dictObj['traj'].xyz[:, dictObj['plot_atom']]
	dictObj['bottom'] = bottom
	if bottom:
		dictObj['plot_z_adj'] = dictObj['surf_avg_z'] - dictObj['plot_z']
	else:
		dictObj['plot_z_adj'] = dictObj['plot_z'] - dictObj['surf_avg_z']

	return dictObj

def extend_dict(dictObj, extend_type, extend_name='2p'):

	d = deepcopy(dictObj)
	# Get atom index
	if type(extend_type) == type([]):
		d['str_%s_type' % extend_name] = list(map(str, extend_type))
	elif type(extend_type) == type(1):
		d['str_%s_type' % extend_name] = str(extend_type)

	d['%s_atom' % extend_name] = [atom.index for atom in d['topo'].atoms if (atom.name in d['str_%s_type' % extend_name])]

	# Get plot atom position
	d['%s_z' % extend_name] = d['traj'].xyz[:, d['%s_atom' % extend_name], 2]
	d['%s_xyz' % extend_name] = d['traj'].xyz[:, d['%s_atom' % extend_name]]

	if d['bottom']:
		d['%s_z_adj' % extend_name] = d['surf_avg_z'] - d['%s_z' % extend_name]
	else:
		d['%s_z_adj' % extend_name] = d['%s_z' % extend_name] - d['surf_avg_z']

	return d

def near_pvp_dict(dictObj, under_z=2.0, search_radius=0.5):

	if not dictObj['2p_xyz'].any():
		raise NameError('Extend the dictionary object to 2P/MP first')

	d = deepcopy(dictObj)
	d['2p_below_xyz'] = list(map( lambda x: filter_below_z(x, d['surf_avg_z'], under_z) , d['2p_xyz'] ))
	d['2p_near_xyz'] = list(map( lambda pair: filter_near_o(pair[0], pair[1], search_radius) , zip(d['plot_xyz'], d['2p_below_xyz']) ))
	d['2p_near_z_adj'] = flatten_list(list(map( lambda x: get_z_adj(x, d['surf_avg_z']), d['2p_near_xyz'])))

	return d

def filter_below_z(tp_xyz, surf_avg_z, under_z):
    return list(filter(lambda xyz: ((xyz[2] - surf_avg_z < under_z) and (xyz[2] - surf_avg_z > 0.0)) , tp_xyz))

def near_o(pvp_xyz, tp_xyz, search_radius):
    diff = list(map(np.linalg.norm, np.subtract(pvp_xyz, tp_xyz)))
    return len([x for x in diff if x < search_radius]) > 0

def filter_near_o(pvp_xyz, tp_xyz, search_radius):
    return list(filter(lambda p_xyz: near_o(pvp_xyz, p_xyz, search_radius), tp_xyz))

def get_z_adj(xyz, surf_avg_z):
    return list(map(lambda x: x[2] - surf_avg_z, xyz))

def flatten_list(the_list):
    return [y for x in the_list for y in x]

def plot_dict(dictObjs, labels, bins=100, f_name='hist_plot', title='', plot_name='plot', xlim=''):
	fig = plt.figure()

	if (type(dictObjs) == type([]) and type(labels) == type([])):
		for dictObj, label in zip(dictObjs, labels):
			count, binEdges = np.histogram(dictObj['%s_z_adj' % plot_name], bins, density=False)
			dist_bin = binEdges[1]-binEdges[0]
			bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
			density = (count)/(dictObj['area']*dist_bin*dictObj['n_frames']) # count per nm^3
			plt.plot(bincenters, density, label=label)
		plt.legend()

	elif (type(dictObjs) == type({}) and type(labels) == type("")):
		count, binEdges = np.histogram(dictObjs['%s_z_adj' % plot_name], bins, density=False)
		dist_bin = binEdges[1]-binEdges[0]
		bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
		density = (count)/(dictObjs['area']*dist_bin*dictObjs['n_frames']) # count per nm^3
		plt.plot(bincenters, density, label=labels)

	if type(xlim) == type([]):
		plt.xlim(xlim[0], xlim[1])

	plt.xlabel('z (nm)')
	plt.ylabel('density (count / nm^3)')
	plt.title(title)
	fig.savefig('density_%s.png' % f_name, dpi=fig.dpi)
	plt.show()
	plt.close()

def filter_z_between(dictObj, start, stop):
	newDictObj = deepcopy(dictObj)
	newDictObj['plot_z_adj'] = list(map( lambda x: list(filter( lambda y: (y >= start) and (y <= stop) , x)), dictObj['plot_z_adj']))
	newDictObj['plot_z_adj'] = [ i for j in newDictObj['plot_z_adj'] for i in j ]
	return newDictObj

def count_switch(dictObj, thres_top, time_interval, time_step=1.0):
	# time_step [=] femtosecond

	d = deepcopy(dictObj)
	d['count_switch_array'] = [mapper_switch(data, thres_top, time_interval, time_step, index) for (index, data) in enumerate(d['plot_z_adj'].transpose())]
	d['count_switch'] = np.sum(d['count_switch_array'])
	d['len_ns'] = d['plot_z_adj'].shape[0]*time_interval*time_step/1000000.
	d['switch_every'] = np.divide(d['len_ns'],d['count_switch_array']) # ns per swith
	d['avg_switch_every'] = d['len_ns']/d['count_switch'] # ns per switch

	print("Run time: %d ns" % d['len_ns'])
	print("Count switch:")
	print(d['count_switch_array'])
	print("Total switch: %d" % d['count_switch'])
	print("")
	print("Switch every (ns/switch):")
	np.set_printoptions(precision=1)
	print(d['switch_every'])
	print("Average switch every: %.1f ns/switch" % d['avg_switch_every'])
	print("")
	print("Atom ID:")
	print(np.add(d['plot_atom'],1))
	print("")
	return d

def mapper_switch(data, thres_top, time_interval, time_step, index):
	ns_step = time_interval*time_step/1000000.
	time = np.linspace(0, len(data)*ns_step, len(data))
	count = 0
	switch_pos = []
	top = data[0] > thres_top
	for (i, val) in enumerate(data):
		if top != (val > thres_top):
			count += 1
			switch_pos.append(i*ns_step)
		top = (val > thres_top)

	fig = plt.figure()
	plt.plot(time, data)
	plt.axhline(y=thres_top, color='k')
	for s in switch_pos:
		plt.axvline(x=s, color='r', linestyle=":")
	plt.xlabel('Time (ns)')
	plt.ylabel('z (nm)')
	plt.title('Atom: %d, Switch Count: %d' % (index+1, count))
	fig.savefig('switch_%d.png' % (index+1), dpi=fig.dpi)
	plt.show()
	plt.close()

	return count
