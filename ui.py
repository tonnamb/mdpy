import numpy as np
import matplotlib.pyplot as plt
import shutil
import os

def plot(end_more = 0, num_windows = None, verticals = []):
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
    for vertical in verticals:
        plt.axvline(x=vertical, color='r', linestyle='--')
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

def plot_traj_ui(num_windows):
    if not os.path.exists('traj_ui'):
        os.makedirs('traj_ui')
    for i in range(1, num_windows+1):
        ui_traj = np.genfromtxt('data/%d' % i)
        fig = plt.figure()
        plt.plot(ui_traj[:, 0])
        plt.plot(ui_traj[:, 2])
        plt.ylabel(r'z ($\AA$)')
        plt.xlabel('Samples')
        plt.title(i)
        fig.savefig('traj_ui/%d.png' % i, dpi=fig.dpi)
        plt.close()

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

def read_write_colvar_auto_k(runrange, prevfolder, row_sample=1, row_restrain=2):
    eVtoKcal = 23.06
    for i in runrange:
        with open(str(i) + "_window/pull_unix.in") as file:
            for line in file:
                line_split = line.split()
                if len(line_split) > 0 and line_split[0] == 'forceConstant':
                    k_eV = float(line_split[1])
        k_kcal = k_eV*eVtoKcal
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
                        out.write(rows[row_sample] + ' {0:.2f} '.format(k_kcal) + rows[row_restrain] + '\n')
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

def extract_colvars_auto_k(
        start, stop, prevfolder=None,
        T=413, min_rc=0, max_rc=30, n=2000, ss=0,
        row_sample=1, row_restrain=2):

    # example prevfolder = "2ns/k_1.7"
    if prevfolder:
        shutil.copytree("../../%s/ui/data" % prevfolder, "ui/data")
    else:
        os.makedirs("ui/data")

    runrange = range(start, stop+1)

    read_write_colvar_auto_k(runrange, prevfolder, row_sample, row_restrain)

    os.chdir("ui")
    # ui.out -ui -T 900 -min 0 -max 70 -n 500 -u kcal -ss 1000 -r -1 -v 2 > log.txt
    os.system("ui.out -ui -T {0:d} -min {1:d} -max {2:d} -n {3:d} -u kcal -ss {4:d} -r -1 -v 2 > log.txt".format(T, min_rc, max_rc, n, ss))

def read_out_colvars_abf(
        fname, nwin, start, step, windowstart, windowsize, fout='abf_range'):

    data = np.genfromtxt(fname)
    boundaries = [windowstart - windowsize * (i - 1) for i in range(1, nwin + 2)]
    arg_start = np.argwhere(data[:, 0] == start)[0, 0]
    arg_step = np.argwhere(data[:, 0] == start + step)[0, 0] - arg_start
    arg_stop = arg_start + arg_step * nwin
    x = data[arg_start:arg_stop:arg_step][:, 0]
    y = data[arg_start:arg_stop:arg_step][:, 1]

    fig = plt.figure()
    plt.scatter(x, y)
    for boundary in boundaries:
        plt.axhline(y=boundary)
    plt.xlabel('Restart timesteps')
    plt.ylabel('Colvars position')
    plt.title('nwin:%d start:%d step:%d windowstart:%d windowsize:%d' % (nwin, start, step, windowstart, windowsize))
    fig.savefig('%s.png' % fout, dpi=fig.dpi)
    plt.show()
    plt.close()

def find_nearest_idx(array, value):
    return (np.abs(array-value)).argmin()

def read_out_colvars_us(
        fname, nwin, windowstart, windowstep, fout='us_range'):

    data = np.genfromtxt(fname)
    centers = [windowstart + windowstep * (i - 1) for i in range(1, nwin + 1)]
    data_idx = [find_nearest_idx(data[:, 1], center) for center in centers]
    x = data[data_idx][:, 0]
    print(x)
    y = data[data_idx][:, 1]

    fig = plt.figure()
    plt.scatter(x, y)
    for center in centers:
        plt.axhline(y=center)
    plt.xlabel('Restart timesteps')
    plt.ylabel('Colvars position')
    plt.title('nwin:%d windowstart:%d windowstep:%d' % (nwin, windowstart, windowstep))
    fig.savefig('%s.png' % fout, dpi=fig.dpi)
    plt.show()
    plt.close()

def read_out_colvars_us_center(fname, centers, fout='us_range'):

    data = np.genfromtxt(fname)
    data_idx = [find_nearest_idx(data[:, 1], center) for center in centers]
    plot_out_colvers(data, data_idx, centers, fout)

def plot_out_colvers(data, data_idx, centers, fout):

    x = data[data_idx][:, 0]
    print(x)
    y = data[data_idx][:, 1]

    fig = plt.figure()
    plt.scatter(x, y, marker='v')
    plt.scatter(x, centers, s=3)
    plt.xlabel('Restart timesteps')
    plt.ylabel('Colvars position')
    fig.savefig('%s.png' % fout, dpi=fig.dpi)
    plt.show()
    plt.close()
