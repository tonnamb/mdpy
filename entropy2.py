"""Calculates the entropy using center of mass of PVP chain as bound/free criteria"""

import chemfiles
from mdpy import thermo
import numpy as np
from mdpy.entropy import pbc_wrap
from mdpy.entropy import in_between
from mdpy.entropy import calc_ns
from scipy.stats import sem
from numpy.linalg import norm
import linecache
import os

def calc_entropy(
        f_thermo, f_traj,
        pvp_types, pvp_masses, ag_type,
        bound_range=[3, 6], free_dist=10, time_interval=2000,
        temp=413, binding_free_kcal=0, binding_free_sem_kcal=0):

    data = {}
    data['bound_pe'] = []
    data['free_pe'] = []

    for thermo_file, traj_file in zip(f_thermo, f_traj):

        loaded = get_bound_free_pe(thermo_file, traj_file,
                                   pvp_types, pvp_masses, ag_type,
                                   bound_range, free_dist, time_interval)

        # loaded['bound_pe'], loaded['free_pe'] assumed to be python list
        data['bound_pe'] = data['bound_pe'] + loaded['bound_pe']
        data['free_pe'] = data['free_pe'] + loaded['free_pe']

    data['bound_pe_mean'] = np.mean(data['bound_pe'])
    data['bound_pe_sem'] = sem(data['bound_pe'])
    data['bound_len'] = calc_ns(data['bound_pe'], time_interval)
    data['free_pe_mean'] = np.mean(data['free_pe'])
    data['free_pe_sem'] = sem(data['free_pe'])
    data['free_len'] = calc_ns(data['free_pe'], time_interval)

    bind_pot = data['bound_pe_mean'] - data['free_pe_mean']
    bind_pot_sem = norm([data['bound_pe_sem'], data['free_pe_sem']])

    eV_to_kcal = 23.0605
    bind_pot_kcal = bind_pot*eV_to_kcal
    bind_pot_sem_kcal = bind_pot_sem*eV_to_kcal

    entropy = (bind_pot_kcal - binding_free_kcal)/temp
    entropy_sem = norm([bind_pot_sem_kcal, binding_free_sem_kcal])/temp

    entro_kcal = bind_pot_kcal - binding_free_kcal
    entro_sem_kcal = norm([bind_pot_sem_kcal, binding_free_sem_kcal])

    print('')
    print('Bound range = [{0}, {1}]'.format(*bound_range))
    print('Free dist = {0}'.format(free_dist))
    print('Binding potential energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(bind_pot_kcal, bind_pot_sem_kcal))
    print('Binding free energy = {0:.2f} +/- {1:.2f} kcal/mol'.format(binding_free_kcal, binding_free_sem_kcal))
    print('Entropic contribution = {0:.2f} +/- {1:.2f} kcal/mol'.format(entro_kcal, entro_sem_kcal))
    print('Entropy at {0} K = {1:.4f} +/- {2:.4f} kcal/mol/T'.format(temp, entropy, entropy_sem))
    print('')
    print('Bound = %.2f ns' % data['bound_len'])
    print('Free = %.2f ns' % data['free_len'])
    print('')

    data['bound_pe_kcal'] = np.array(data['bound_pe'])*eV_to_kcal
    data['free_pe_kcal'] = np.array(data['free_pe'])*eV_to_kcal

    print('Bound potential energy segments')
    thermo.segment_plot(data['bound_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='bound_pe', unit='kcal/mol')
    print('')

    print('Free potential energy segments')
    thermo.segment_plot(data['free_pe_kcal'], skip_ns=0, time_interval=time_interval, n_seg=20, name='free_pe', unit='kcal/mol')
    print('')

    return data

def get_bound_free_pe(
        f_thermo, f_traj,
        pvp_types, pvp_masses, ag_type,
        bound_range, free_dist, time_interval):

    data, head = thermo.read_str_array(f_thermo)

    if head != "# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz":
        print("WARNING: Header assumption is FALSE!")
        raise NameError('WARNING: Header assumption is FALSE!')

    # assumes thermo.lammps is
    # Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz format
    data = np.array(data)
    pot_eng = data[:, 1]

    # Load trajectory file
    print(f_traj)

    pvp_atoms = []
    pvp_com_z = []
    agbot_z = []
    agtop_z = []

    traj = chemfiles.Trajectory(f_traj)
    frame = traj.read()

    for pvp_type in pvp_types:
        pvp_atoms.append(get_atoms_id(pvp_type, frame))
    ag_atom = get_atoms_id(ag_type, frame)

    box_z = frame.cell().lengths()[2]

    positions = frame.positions()

    pvp_com_z.append(get_pvp_com_z(positions, pvp_atoms, pvp_masses))
    agbot_z.append(np.min(positions[ag_atom, 2]))
    agtop_z.append(np.max(positions[ag_atom, 2]))

    for i in range(traj.nsteps()-1):
        frame = traj.read()
        positions = frame.positions()

        pvp_com_z.append(get_pvp_com_z(positions, pvp_atoms, pvp_masses))
        agbot_z.append(np.min(positions[ag_atom, 2]))
        agtop_z.append(np.max(positions[ag_atom, 2]))

    pvp_com_z = [pbc_wrap(x, box_z) for x in pvp_com_z]

    agbot_z = np.array(agbot_z)
    agtop_z = np.array(agtop_z)

    # Bound zone
    bound_bottom_low = agbot_z - bound_range[1]
    bound_bottom_high = agbot_z - bound_range[0]
    bound_top_low = agtop_z + bound_range[0]
    bound_top_high = agtop_z + bound_range[1]

    # Get bound index
    bound_bottom = in_between(pvp_com_z, bound_bottom_low, bound_bottom_high)
    bound_top = in_between(pvp_com_z, bound_top_low, bound_top_high)
    bound_combine = bound_top + bound_bottom

    # Get free index
    free_bottom_low = agbot_z - free_dist
    free_top_high = agtop_z + free_dist
    free_bound = in_between(pvp_com_z, free_bottom_low, free_top_high)
    free_combine = get_reverse_bool(free_bound)

    bound_index = get_index_from_bool(bound_combine)
    free_index = get_index_from_bool(free_combine)

    # Get potential energy
    bound_pot = list(pot_eng[bound_index])
    free_pot = list(pot_eng[free_index])

    # Output dictionary object
    out_dict = {}
    out_dict['bound_pe'] = bound_pot
    out_dict['free_pe'] = free_pot

    return out_dict

def get_pvp_com_z(positions, pvp_atoms, pvp_masses):

    com_numerator = 0
    com_denom = 0

    for index, pvp_atom in enumerate(pvp_atoms):
        pvp_pos = positions[pvp_atom, 2]
        atom_mass = pvp_masses[index]
        for pos in pvp_pos:
            com_numerator += pos * atom_mass
            com_denom += atom_mass

    return com_numerator/com_denom

def get_atoms_id(atom_type, frame):
    selection = "atoms: name {0:d}".format(atom_type)
    return chemfiles.Selection(selection).evaluate(frame)

def get_index_from_bool(bool_array):
    return np.where(bool_array)[0]

def get_reverse_bool(bool_array):
    return [not x for x in bool_array]

def check_lammpstrj():
    thermo_timestep_start = int(linecache.getline('thermo.lammps', 2).strip().split()[0])
    thermo_timestep_end = int(os.popen('tail -n 1 thermo.lammps').read().strip().split()[0])
    linecache.clearcache()
    
    pvpag_timestep_start = int(linecache.getline('pvpag.lammpstrj', 2).strip())
    n_atoms = int(linecache.getline('pvpag.lammpstrj', 4).strip())
    pvpag_timestep_end = int(os.popen('tail -n %d pvpag.lammpstrj | sed -n 2p' % (n_atoms + 9)).read().strip())
    linecache.clearcache()
    
    if pvpag_timestep_start == thermo_timestep_start:
        print('start: same')
    elif (pvpag_timestep_start + 2000) == thermo_timestep_start:
        print('start: pvpag is one step behind')
        os.system("sed -i.backup -e '1,{0:d}d' pvpag.lammpstrj".format(n_atoms + 9))
        print('start: removed first step in pvpag')
        # re-check
        pvpag_timestep_start = int(linecache.getline('pvpag.lammpstrj', 2).strip())
        linecache.clearcache()
        if pvpag_timestep_start == thermo_timestep_start:
            print('start: same')
        else:
            print('start: different')
            print('start: thermo = %d' % thermo_timestep_start)
            print('start: pvpag = %d' % pvpag_timestep_start)
            raise Exception('Starts are different')
    else:
        print('start: different')
        print('start: thermo = %d' % thermo_timestep_start)
        print('start: pvpag = %d' % pvpag_timestep_start)
        raise Exception('Starts are different')
    
    if pvpag_timestep_end == thermo_timestep_end:
        print('end: same')
    else:
        print('end: different')
        print('end: thermo = %d' % thermo_timestep_end)
        print('end: pvpag = %d' % pvpag_timestep_end)
        raise Exception('Ends are different')

def check_pvpag_recursive():
    pvpag_paths = [os.path.abspath(os.path.dirname(path)) for path in os.popen('find . -name \pvpag.lammpstrj').read().strip('/pvpag.lammpstrj').split()]
    for path in pvpag_paths:
        os.chdir(path)
        print(path)
        check_lammpstrj()

def get_pvpag_list(postpend_path=''):
    pvpag_paths = [os.path.abspath(os.path.dirname(path)) for path in os.popen('find . -name \pvpag.lammpstrj').read().strip('/pvpag.lammpstrj').split()]
    return ['{0}/{1}'.format(path,postpend_path) for path in pvpag_paths]

