import numpy as np
import mdtraj
from scipy.stats import sem
from numpy.linalg import norm
from mdpy import thermo
import sys
import time
import chemfiles

def read_and_check(fname, skipns=0, time_interval=2000):

	skipts = int(skipns*1000000./time_interval)
	data, head = thermo.read_str_array(fname, skipts)
	print(fname)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
		print("")
	else:
		print("WARNING: Header assumption is FALSE!")
		raise NameError('WARNING: Header assumption is FALSE!')

	data = np.array(data)

	out_dict = {}
	out_dict['pe'] = data[:, 1]
	out_dict['pe_mean'] = np.mean(data[:, 1])
	out_dict['pe_sem'] = sem(data[:, 1])
	out_dict['len_ns'] = calc_ns(data, time_interval)

	return out_dict

def calc_ns(data, time_interval=2000):
	return len(data)*time_interval/1000000.

def vacuum(
	bound_fname, ag_fname, free_fname, 
	bound_skipns=0, ag_skipns=0, free_skipns=0, 
	time_interval=2000, temp=413,
	binding_free_kcal=0, binding_free_sem_kcal=0):

	print("Reading bound")
	bound = read_and_check(bound_fname, bound_skipns, time_interval)

	print("Reading ag")
	ag = read_and_check(ag_fname, ag_skipns, time_interval)

	print("Reading free")
	free = read_and_check(free_fname, free_skipns, time_interval)

	bind_pot = bound['pe_mean'] - (free['pe_mean'] + ag['pe_mean'])
	bind_pot_sem = norm([bound['pe_sem'], free['pe_sem'], ag['pe_sem']])

	eV_to_kcal = 23.0605
	bind_pot_kcal = bind_pot*eV_to_kcal
	bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

	entropy = (bind_pot_kcal - binding_free_kcal)/temp
	entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

	entro_kcal = bind_pot_kcal - binding_free_kcal
	entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

	print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
	print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
	print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
	print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
	print('')
	print('Bound = %.2f ns (skipped %i ns)' % (bound['len_ns'], bound_skipns))
	print('Ag = %.2f ns (skipped %i ns)' % (ag['len_ns'], ag_skipns))
	print('Free = %.2f ns (skipped %i ns)' % (free['len_ns'], free_skipns))
	print('')

def vacuum_eg(
	bound_fname, ag_fname, free_fname, 
	bound_traj, bound_pdb,
	safe_dist=0.5, oxy_type=4, ag_type=5, 
	bool_polymer=True, extra_skip_traj=0,
	bound_skipns=0, ag_skipns=0, free_skipns=0, 
	time_interval=2000, temp=413,
	binding_free_kcal=0, binding_free_sem_kcal=0):

	print("Reading free")
	free = read_and_check(free_fname, free_skipns, time_interval)

	print("Reading ag")
	ag = read_and_check(ag_fname, ag_skipns, time_interval)

	print("Reading bound")
	bound = read_bound_traj(bound_fname, bound_traj, bound_pdb, 
												  safe_dist, oxy_type, ag_type,
												  bound_skipns, time_interval, 
												  bool_polymer, extra_skip_traj)

	bind_pot = bound['pe_mean'] - (free['pe_mean'] + ag['pe_mean'])
	bind_pot_sem = norm([bound['pe_sem'], free['pe_sem'], ag['pe_sem']])

	eV_to_kcal = 23.0605
	bind_pot_kcal = bind_pot*eV_to_kcal
	bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

	entropy = (bind_pot_kcal - binding_free_kcal)/temp
	entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

	entro_kcal = bind_pot_kcal - binding_free_kcal
	entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

	print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
	print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
	print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
	print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
	print('')
	print('Bound = %.2f ns (skipped %i ns)' % (bound['len_ns'], bound_skipns))
	print('Ag = %.2f ns (skipped %i ns)' % (ag['len_ns'], ag_skipns))
	print('Free = %.2f ns (skipped %i ns)' % (free['len_ns'], free_skipns))
	print('')

def read_bound_traj(
	bound_fname, bound_traj, bound_pdb,
	safe_dist=0.5, oxy_type=4, ag_type=5,
	bound_skipns=0, time_interval=2000, 
	bool_polymer=True, extra_skip_traj=0):

	skipts = int(bound_skipns*1000000./time_interval)
	data, head = thermo.read_str_array(bound_fname, skipts)
	print(bound_fname)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
	else:
		print("WARNING: Header assumption is FALSE!")
		raise NameError('WARNING: Header assumption is FALSE!')

	data = np.array(data)
	print("Frames in thermo: %i" % len(data))
	print("")

	# Load trajectory file
	print(bound_traj)
	print("")
	t = mdtraj.load(bound_traj, discard_overlapping_frames=True, top=bound_pdb)
	t = t[extra_skip_traj:]
	oxy_atom = [atom.index for atom in t.topology.atoms if (atom.name == str(oxy_type))]
	ag_atom = [atom.index for atom in t.topology.atoms if (atom.name == str(ag_type))]
	if bool_polymer:
		oxy_z = t.xyz[skipts:, oxy_atom, 2]
	else:
		oxy_z = t.xyz[skipts:, oxy_atom, 2].flatten()
	ag_z = t.xyz[skipts:, ag_atom, 2].flatten()
	print("Frames in traj: %i" % len(oxy_z))

	# Find safe zone
	safe_low = ag_z.min() - safe_dist
	safe_high = ag_z.max() + safe_dist
	print('safe_dist = %.1f nm' % safe_dist)
	print('safe_low = %.2f nm' % safe_low)
	print('safe_high = %.2f nm' % safe_high)
	print("")

	# Get safe index
	if bool_polymer:
		safe_bottom = (oxy_z.min(axis=1) > safe_low) * (oxy_z.max(axis=1) < ag_z.min())
		safe_top = (oxy_z.max(axis=1) < safe_high) * (oxy_z.min(axis=1) > ag_z.max())
	else:
		safe_bottom = (oxy_z > safe_low) * (oxy_z < ag_z.min())
		safe_top = (oxy_z < safe_high) * (oxy_z > ag_z.max())
	safe_traj = np.where(safe_top + safe_bottom)[0]

	# Get safe potential energy
	safe_pot = list(data[:, 1][safe_traj])

	# Output dictionary object
	out_dict = {}
	out_dict['pe_mean'] = np.mean(safe_pot)
	out_dict['pe_sem'] = sem(safe_pot)
	out_dict['len_ns'] = calc_ns(safe_pot, time_interval)

	return out_dict

def load_bound_free(
	f_thermo, f_traj, f_pdb,
	oxy_type, ag_type, bound_dist=0.5,
	skipns=0, time_interval=2000, 
	bool_polymer=True, extra_skip_traj=0):

	skipts = int(skipns*1000000./time_interval)
	data, head = thermo.read_str_array(f_thermo, skipts)
	print(f_thermo)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
	else:
		print("WARNING: Header assumption is FALSE!")
		raise NameError('WARNING: Header assumption is FALSE!')

	data = np.array(data)
	print("Frames in thermo: %i" % len(data))
	print("")

	# Load trajectory file
	print(f_traj)
	print("")
	t = mdtraj.load(f_traj, discard_overlapping_frames=True, top=f_pdb)
	t = t[extra_skip_traj:]
	oxy_atom = [atom.index for atom in t.topology.atoms if (atom.name == str(oxy_type))]
	ag_atom = [atom.index for atom in t.topology.atoms if (atom.name == str(ag_type))]
	if bool_polymer:
		oxy_z = t.xyz[skipts:, oxy_atom, 2]
	else:
		oxy_z = t.xyz[skipts:, oxy_atom, 2].flatten()
	ag_z = t.xyz[skipts:, ag_atom, 2].flatten()
	print("Frames in traj: %i" % len(oxy_z))

	# Find bound zone
	bound_low = ag_z.min() - bound_dist
	bound_high = ag_z.max() + bound_dist
	print('bound_dist = %.1f nm' % bound_dist)
	print('bound_low = %.2f nm' % bound_low)
	print('bound_high = %.2f nm' % bound_high)
	print("")

	# Get bound index
	if bool_polymer:
		bound_bottom = (oxy_z.min(axis=1) > bound_low) * (oxy_z.max(axis=1) < ag_z.min())
		bound_top = (oxy_z.max(axis=1) < bound_high) * (oxy_z.min(axis=1) > ag_z.max())
	else:
		bound_bottom = (oxy_z > bound_low) * (oxy_z < ag_z.min())
		bound_top = (oxy_z < bound_high) * (oxy_z > ag_z.max())

	bound_combine = bound_top + bound_bottom
	free_combine = [not x for x in bound_combine]
	bound_index = np.where(bound_combine)[0]
	free_index = np.where(free_combine)[0]

	# Get potential energy
	bound_pot = list(data[:, 1][bound_index])
	free_pot = list(data[:, 1][free_index])

	# Output dictionary object
	out_dict = {}
	out_dict['bound_pe'] = bound_pot
	out_dict['bound_pe_mean'] = np.mean(bound_pot)
	out_dict['bound_pe_sem'] = sem(bound_pot)
	out_dict['bound_len'] = calc_ns(bound_pot, time_interval)
	out_dict['free_pe'] = free_pot
	out_dict['free_pe_mean'] = np.mean(free_pot)
	out_dict['free_pe_sem'] = sem(free_pot)
	out_dict['free_len'] = calc_ns(free_pot, time_interval)

	return out_dict

def solvent(
	thermo_traj, f_pdb,
	oxy_type, ag_type, bound_dist=0.5,
	time_interval=2000, bool_polymer=True,
	temp=413, binding_free_kcal=0, binding_free_sem_kcal=0):

	d = {}
	d['bound_pe'] = []
	d['free_pe'] = []
	for f_thermo, f_traj in thermo_traj:
		loaded = load_bound_free(f_thermo, f_traj, f_pdb,
														 oxy_type, ag_type, bound_dist=bound_dist,
														 time_interval=time_interval, bool_polymer=bool_polymer)
		d['bound_pe'] = d['bound_pe'] + loaded['bound_pe']
		d['free_pe'] = d['free_pe'] + loaded['free_pe']

	d['bound_pe_mean'] = np.mean(d['bound_pe'])
	d['bound_pe_sem'] = sem(d['bound_pe'])
	d['bound_len'] = calc_ns(d['bound_pe'], time_interval)
	d['free_pe_mean'] = np.mean(d['free_pe'])
	d['free_pe_sem'] = sem(d['free_pe'])
	d['free_len'] = calc_ns(d['free_pe'], time_interval)

	bind_pot = d['bound_pe_mean'] - d['free_pe_mean']
	bind_pot_sem = norm([d['bound_pe_sem'], d['free_pe_sem']])

	eV_to_kcal = 23.0605
	bind_pot_kcal = bind_pot*eV_to_kcal
	bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

	entropy = (bind_pot_kcal - binding_free_kcal)/temp
	entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

	entro_kcal = bind_pot_kcal - binding_free_kcal
	entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

	print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
	print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
	print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
	print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
	print('')
	print('Bound = %.2f ns' % d['bound_len'])
	print('Free = %.2f ns' % d['free_len'])
	print('')

	d['bound_pe_kcal'] = np.array(d['bound_pe'])*eV_to_kcal
	d['free_pe_kcal'] = np.array(d['free_pe'])*eV_to_kcal

	print('Bound potential energy segments')
	thermo.segment_plot(d['bound_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='bound_pe', unit='kcal/mol')
	print('')

	print('Free potential energy segments')
	thermo.segment_plot(d['free_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='free_pe', unit='kcal/mol')
	print('')

	return d

def minus_one(val):
  return val-1

def solvent_traj(
	thermo_traj,
	oxy_type, agbot_id, agtop_id,
	adjust_ag_index_please=False,
	bound_dist=7, time_interval=2000, bool_polymer=False,
	temp=413, binding_free_kcal=0, binding_free_sem_kcal=0):

	if bool_polymer:
		raise NameError('ERROR: Do not support polymer yet!')

	if adjust_ag_index_please:
		agbot_id = list(map(minus_one, agbot_id))
		agtop_id = list(map(minus_one, agtop_id))

	d = {}
	d['bound_pe'] = []
	d['free_pe'] = []
	for f_thermo, f_traj in thermo_traj:
		loaded = chemfiles_bound_free(f_thermo, f_traj,
																	oxy_type, agbot_id, agtop_id, 
																	bound_dist=bound_dist,
																	time_interval=time_interval, bool_polymer=bool_polymer)
		d['bound_pe'] = d['bound_pe'] + loaded['bound_pe']
		d['free_pe'] = d['free_pe'] + loaded['free_pe']
	
	d['bound_pe_mean'] = np.mean(d['bound_pe'])
	d['bound_pe_sem'] = sem(d['bound_pe'])
	d['bound_len'] = calc_ns(d['bound_pe'], time_interval)
	d['free_pe_mean'] = np.mean(d['free_pe'])
	d['free_pe_sem'] = sem(d['free_pe'])
	d['free_len'] = calc_ns(d['free_pe'], time_interval)

	bind_pot = d['bound_pe_mean'] - d['free_pe_mean']
	bind_pot_sem = norm([d['bound_pe_sem'], d['free_pe_sem']])

	eV_to_kcal = 23.0605
	bind_pot_kcal = bind_pot*eV_to_kcal
	bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

	entropy = (bind_pot_kcal - binding_free_kcal)/temp
	entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

	entro_kcal = bind_pot_kcal - binding_free_kcal
	entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

	print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
	print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
	print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
	print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
	print('')
	print('Bound = %.2f ns' % d['bound_len'])
	print('Free = %.2f ns' % d['free_len'])
	print('')

	d['bound_pe_kcal'] = np.array(d['bound_pe'])*eV_to_kcal
	d['free_pe_kcal'] = np.array(d['free_pe'])*eV_to_kcal

	print('Bound potential energy segments')
	thermo.segment_plot(d['bound_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='bound_pe', unit='kcal/mol')
	print('')

	print('Free potential energy segments')
	thermo.segment_plot(d['free_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='free_pe', unit='kcal/mol')
	print('')

	return d

def pbc_wrap(val, box_size):
  if (val >= 0 and val <= box_size):
    return val

  elif val < 0:
    newval = val
    while newval < 0:
      newval = newval + box_size
    return newval

  elif val > box_size:
    newval = val
    while newval > box_size:
      newval = newval - box_size
    return newval

def chemfiles_bound_free(
	f_thermo, f_traj,
	oxy_type, agbot_id, agtop_id, 
	bound_dist=7,
	skipns=0, time_interval=2000, bool_polymer=False):

	if bool_polymer:
		raise NameError('ERROR: Do not support polymer yet!')

	skipts = int(skipns*1000000./time_interval)
	data, head = thermo.read_str_array(f_thermo, skipts)
	print(f_thermo)

	if (head == "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz"):
		print("Header assumption is TRUE")
		print("# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz")
	else:
		print("WARNING: Header assumption is FALSE!")
		raise NameError('WARNING: Header assumption is FALSE!')

	data = np.array(data)
	print("Frames in thermo: %i" % len(data))
	print("")

	# Load trajectory file
	print(f_traj)
	print("")
	traj = chemfiles.Trajectory(f_traj)
	frame = traj.read()
	oxy_atom = chemfiles.Selection("atoms: name == %s" % oxy_type).evaluate(frame)
	box_z = frame.cell().lengths()[2]

	oxy_z = []
	agbot_z = []
	agtop_z = []

	positions = frame.positions()

	if bool_polymer:
		oxy_z.append( list(positions[oxy_atom, 2]) )
	else:
		oxy_z.append( positions[oxy_atom, 2][0] )
	
	agbot_z.append( np.min(positions[agbot_id, 2]) )
	agtop_z.append( np.max(positions[agtop_id, 2]) )

	for i in range(traj.nsteps()-1):
		frame = traj.read()
		positions = frame.positions()
		if bool_polymer:
			oxy_z.append( list(positions[oxy_atom, 2]) )
		else:
			oxy_z.append( positions[oxy_atom, 2][0] )
		agbot_z.append( np.min(positions[agbot_id, 2]) )
		agtop_z.append( np.max(positions[agtop_id, 2]) )

	oxy_z = oxy_z[skipts:]
	agbot_z = agbot_z[skipts:]
	agtop_z = agtop_z[skipts:]
	
	oxy_z = list(map(lambda x:pbc_wrap(x, box_z), oxy_z))

	print("Frames in traj: %i" % len(oxy_z))
	print("")

	# Bound zone
	agbot_z = np.array(agbot_z)
	agtop_z = np.array(agtop_z)
	bound_low = agbot_z - bound_dist
	bound_high = agtop_z + bound_dist
	
	# Get bound index
	if bool_polymer:
		raise NameError('ERROR: Do not support polymer yet!')
		# bound_bottom = (oxy_z.min(axis=1) > bound_low) * (oxy_z.max(axis=1) < ag_z.min())
		# bound_top = (oxy_z.max(axis=1) < bound_high) * (oxy_z.min(axis=1) > ag_z.max())
	else:
		bound_bottom = (oxy_z > bound_low) * (oxy_z < agbot_z)
		bound_top = (oxy_z < bound_high) * (oxy_z > agtop_z)

	bound_combine = bound_top + bound_bottom
	free_combine = [not x for x in bound_combine]
	bound_index = np.where(bound_combine)[0]
	free_index = np.where(free_combine)[0]

	# Get potential energy
	bound_pot = list(data[:, 1][bound_index])
	free_pot = list(data[:, 1][free_index])

	# Output dictionary object
	out_dict = {}
	out_dict['bound_pe'] = bound_pot
	out_dict['bound_pe_mean'] = np.mean(bound_pot)
	out_dict['bound_pe_sem'] = sem(bound_pot)
	out_dict['bound_len'] = calc_ns(bound_pot, time_interval)
	out_dict['free_pe'] = free_pot
	out_dict['free_pe_mean'] = np.mean(free_pot)
	out_dict['free_pe_sem'] = sem(free_pot)
	out_dict['free_len'] = calc_ns(free_pot, time_interval)

	return out_dict

def solvent_no_traj(
	bound_fname, free_fname, 
	bound_skipns=0, free_skipns=0, 
	time_interval=2000, temp=413,
	binding_free_kcal=0, binding_free_sem_kcal=0):

	print("Reading bound")
	bound = read_and_check(bound_fname, bound_skipns, time_interval)

	print("Reading free")
	free = read_and_check(free_fname, free_skipns, time_interval)

	bind_pot = bound['pe_mean'] - free['pe_mean']
	bind_pot_sem = norm([bound['pe_sem'], free['pe_sem']])

	eV_to_kcal = 23.0605
	bind_pot_kcal = bind_pot*eV_to_kcal
	bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

	entropy = (bind_pot_kcal - binding_free_kcal)/temp
	entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

	entro_kcal = bind_pot_kcal - binding_free_kcal
	entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

	print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
	print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
	print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
	print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
	print('')
	print('Bound = %.2f ns (skipped %i ns)' % (bound['len_ns'], bound_skipns))
	print('Free = %.2f ns (skipped %i ns)' % (free['len_ns'], free_skipns))
	print('')

	bound['pe_kcal'] = bound['pe']*eV_to_kcal
	free['pe_kcal'] = free['pe']*eV_to_kcal

	print('Bound potential energy segments')
	thermo.segment_plot(bound['pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='bound_pe', unit='kcal/mol')
	print('')

	print('Free potential energy segments')
	thermo.segment_plot(free['pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='free_pe', unit='kcal/mol')
	print('')

class Logger:
    def __init__(self, fname):
        self.terminal = sys.stdout
        self.log = open(fname, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

def start_logger(fname):
	sys.stdout = Logger(fname)
	print(time.strftime("%c"))
	print('')
