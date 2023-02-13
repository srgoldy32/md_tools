import atomman as am
import numpy as np
import time
import sys

if len(sys.argv) != 4:
    print('Incorrect arguments: [ref filename, end filename, timestep]')
else:
    ref_filename = sys.argv[1]
    end_filename = sys.argv[2]
    timestep = int(sys.argv[3])
    split = end_filename.split(".")
    descriptor = split[1]
    end = int(split[-1])
    start = int(ref_filename.split(".")[-1])
    times = range(timestep,end+timestep,timestep)

t_0 = time.time()
d0 = am.load('atom_dump', ref_filename, symbols='Al')

ne11 = am.NeighborList(system=d0, cutoff=4.07)
ne12 = am.NeighborList(system=d0, cutoff=3.52)
ne23 = am.NeighborList(system=d0, cutoff=2.86)
tol = 2.86/np.sqrt(3)/2

load = time.time()
load_time = round(load-t_0,3)
print("Load time: ",load_time)
for t in times:
    st = time.time()
    def_filename = 'out.{}.{}'.format(descriptor,t)

    dd = am.load('atom_dump', def_filename, symbols='Al')
    dd.atoms.slip_vec = np.zeros([1,3])
    dd.atoms.slip_mag = 0.0
    for ind in range(len(d0)):
        atype = d0.atoms[ind].atype
        l1 = ne11[ind]
        l2 = ne12[ind]
        l3 = ne23[ind]
        if atype == [1]:
            types_1 = d0.atoms.prop(key='atype',index=l1)
            types_2 = d0.atoms.prop(key='atype',index=l2)
            ones = np.nonzero(types_1 == 1)
            twos = np.nonzero(types_2 == 2)
            one_ne = l1[ones]
            two_ne = l2[twos]
            nes = np.concatenate((one_ne,two_ne))
        if atype == [2]:
            types_2 = d0.atoms.prop(key='atype',index=l2)
            types_3 = d0.atoms.prop(key='atype',index=l3)
            twos = np.nonzero(types_2 == 1)
            thrs = np.nonzero(types_3 == 3)
            one_ne = l2[twos]
            thr_ne = l3[thrs]
            nes = np.concatenate((one_ne,thr_ne))
        if atype == [3]:
            types_2 = d0.atoms.prop(key='atype',index=l3)
            twos = np.nonzero(types_2 == 2)
            nes = l3[twos]
        
        me_r = d0.atoms.prop(key='pos',index=ind)
        me_d = dd.atoms.prop(key='pos',index=ind)
        r = d0.atoms.prop(key='pos',index=nes)
        d = dd.atoms.prop(key='pos',index=nes)

        del_r = am.dvect(me_r,r,d0.box,d0.pbc)
        del_d = am.dvect(me_d,d,dd.box,dd.pbc)

        full_del = del_d - del_r

        slipmag = np.linalg.norm(full_del,axis=1)

        check = slipmag > tol
        inds = np.nonzero(check)
        ns = len(inds)


        slip_vec = full_del[inds].sum(axis=0)/-ns
        slip_mag = np.linalg.norm(slip_vec)
        
        dd.atoms.prop('slip_vec',ind,slip_vec)
        dd.atoms.prop('slip_mag',ind,slip_mag)


    atom_dump = dd.dump('atom_dump')
    slip_filename = 'out.{}.slip.{}_me'.format(descriptor,timestep)
    f = open(slip_filename, "w")
    f.write(atom_dump)
    done = time.time()
    run_time = round(done-st,3)
    print("Wrote: {} in {} s".format(slip_filename,run_time))


        
        
