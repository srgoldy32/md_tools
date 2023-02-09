import sys
import atomman as am
import numpy as np
import time

tol=2.98/np.sqrt(3)/2;

def get_slip_1s(i):
    nlist = ne11[i]
    ref1 = d0.atoms[i].pos
    def1 = dd.atoms[i].pos

    ns = 0
    slip_vec = np.array([[0.0,0.0,0.0]])
    
    if nlist.size > 0:
        for j in nlist:
            on = False
            ref2 = d0.atoms[j].pos
            def2 = dd.atoms[j].pos
            dmag = d0.dmag(ref1,ref2)
            if j in twos:
                if dmag < 3.52:
                    on = True
            elif j in ones:
                on = True
            if on:
                del_ref = am.dvect(ref1,ref2,d0.box,d0.pbc)
                del_def = am.dvect(def1,def2,dd.box,dd.pbc)
                full_del = del_def - del_ref
                slipmag = np.linalg.norm(full_del)
                if slipmag > tol:
                    ns += 1
                    slip_vec += full_del
        if ns == 0:
            slip_vec = np.array([[0.0,0.0,0.0]])
        else:
            slip_vec = slip_vec/-ns
    slip_mag = np.linalg.norm(slip_vec)
    dd.atoms.prop('slip_vec',i,slip_vec)
    dd.atoms.prop('slip_mag',i,slip_mag)
    return slip_vec, slip_mag

def get_slip_2s(i):
    nlist = ne12[i]
    ref1 = d0.atoms[i].pos
    def1 = dd.atoms[i].pos

    atype1 = d0.atoms[i].atype
    ns = 0
    slip_vec = np.array([[0.0,0.0,0.0]])
    
    if nlist.size > 0:
        for j in nlist:
            on = False
            ref2 = d0.atoms[j].pos
            def2 = dd.atoms[j].pos
            dmag = d0.dmag(ref1,ref2)
            if j in threes:
                if dmag < 2.86:
                    on = True
            if j in ones:
                on = True

            if on:
                del_ref = am.dvect(ref1,ref2,d0.box,d0.pbc)
                del_def = am.dvect(def1,def2,dd.box,dd.pbc)
                full_del = del_def - del_ref
                
                
                slipmag = np.linalg.norm(full_del)
                if slipmag > tol:
                    ns += 1
                    slip_vec += full_del
        if ns == 0:
            slip_vec = np.array([[0.0,0.0,0.0]])
        else:
            slip_vec = slip_vec/-ns
    slip_mag = np.linalg.norm(slip_vec)
    dd.atoms.prop('slip_vec',i,slip_vec)
    dd.atoms.prop('slip_mag',i,slip_mag)
    return slip_vec, slip_mag

def get_slip_3s(i):
    nlist = ne23[i]
    ref1 = d0.atoms[i].pos
    def1 = dd.atoms[i].pos
    ns = 0
    slip_vec = np.array([[0.0,0.0,0.0]])
    
    if nlist.size > 0:
        for j in nlist:
            if j in twos:
                ref2 = d0.atoms[j].pos
                def2 = dd.atoms[j].pos

                del_ref = am.dvect(ref1,ref2,d0.box,d0.pbc)
                del_def = am.dvect(def1,def2,dd.box,dd.pbc)
                full_del = del_def - del_ref


                slipmag = np.linalg.norm(full_del)
                if slipmag > tol:
                    ns += 1
                    slip_vec += full_del
        if ns == 0:
            slip_vec = np.array([[0.0,0.0,0.0]])
        else:
            slip_vec = slip_vec/-ns
    slip_mag = np.linalg.norm(slip_vec)
    dd.atoms.prop('slip_vec',i,slip_vec)
    dd.atoms.prop('slip_mag',i,slip_mag)
    return slip_vec, slip_mag

def parse_ref(ref_filename):
    d0 = am.load('atom_dump', ref_filename, symbols='Al')
    ne11 = am.NeighborList(system=d0, cutoff=4.07)
    ne12 = am.NeighborList(system=d0, cutoff=3.52)
    ne23 = am.NeighborList(system=d0, cutoff=2.86)
    d0.atoms.slip_vec = np.zeros([1,3])
    d0.atoms.slip_mag = 0.0

    ones = []
    twos = []
    threes = []
    for i in range(len(d0.atoms)):
        atype = d0.atoms[i].atype
        
        if atype == [1]:
            ones.append(i)
        if atype == [2]:
            twos.append(i)
        if atype == [3]:
            threes.append(i)
    

    return d0,ne11,ne12,ne23,ones,twos,threes

def parse_dump(descriptor,timestep,d0,ne11,ne12,ne23,ones,twos,threes):
    st = time.time()
    dump_filename = 'out.{}.{}'.format(descriptor,timestep)
    dd = am.load('atom_dump', dump-filename, symbols='Al')
    dd.atoms.slip_vec = np.zeros([1,3])
    dd.atoms.slip_mag = 0.0

    for one in ones:
        get_slip_1s(one)

    for two in twos:
        get_slip_2s(two)
        
    for thr in threes:
        get_slip_3s(thr)

    atom_dump = d5.dump('atom_dump')
    slip_filename = 'out.{}.slip.{}'.format(descriptor,timestep)
    f = open(filename, "w")
    f.write(atom_dump)
    done = time.time
    run_time = round(done-st,3)
    print("Wrote: {} in {} s".format(filename,run_time))

    


def driver_func(ref_filename, end_filename, timestep):
    
    split = end_filename.split(".")
    descriptor = split[1]
    end = int(split[-1])
    start = int(ref_filename.split(".")[-1])
    times = range(timestep,end+timestep,timestep)

    t_0 = time.time()
    d0,ne11,ne12,ne23,ones,twos,threes = parse_ref(ref_filename)
    t_ref = time.time()
    print("Ref Load time: ",rounad(t_ref-t_0,3))


    for t in times:
        parse_dump(descriptor,timestep,d0,ne11,ne12,ne23,ones,twos,threes)







if len(sys.argv) != 2:
    print('Incorrect arguments: [ref filename, end filename, timestep]')
else:
    ref_filename = sys.argv[1]
    end_filename = sys.argv[2]
    timestep = sys.argv[3]
    driver_func(ref_filename,end_filename,timestep)