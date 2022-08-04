import sys
import matplotlib.pyplot as plt
import numpy as np


def thermal_eq_stats(filename,plot=False):
    with open(filename) as f:
        lines = f.readlines()

    thermal_eq_flag = False
    dump_interval = 1.1
    for i in range(len(lines)):
        spl = lines[i].split()
        if spl == ['#', '----------', 'Thermal', 'Equilibriation', '---------------------']:
            thermal_eq_flag = True
        if thermal_eq_flag:
            if spl == ['Step', 'PotEng', 'Temp']:
                start = i+1
        
        if len(spl) > 0:
            if spl[0] == 'run':
                if spl[1][0] != '$':
                    length = int(spl[1])
            if spl[0] == 'dump':
                dump_interval = int(spl[4])
                
    dt = int(lines[start+1].split()[0])


    end = start + length//dt
    ts = []
    temps = []
    dump_t = []
    dump_temp = []
    for line in lines[start:end+1]:
        l = line.split()
        t = float(l[0])
        temp = float(l[-1])
        ts.append(t)
        temps.append(temp)
        if t % dump_interval == 0:
            dump_t.append(t)
            dump_temp.append(temp)

    last = temps[-101:]
    av = np.average(last)
    std = np.std(last)

    print('Av: {}; std: {}'.format(round(av,2),round(std,2)))

    if plot == True:
        plt.plot(ts,temps)
        plt.scatter(dump_t,dump_temp)
        plt.xlabel('Timestep')
        plt.ylabel('Temperature K')
        plt.title('{} - last datapoints average: {}; std: {}'.format(filename,round(av,2),round(std,2)))
        plt.show()
    

if len(sys.argv) != 3:
    print('Incorrect arguments: [filename,plot]')
else:
    filename = sys.argv[1]
    plot = bool(sys.argv[2])
    thermal_eq_stats(filename,plot)
