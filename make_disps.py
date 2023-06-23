import atomman as am
import numpy as np
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


ref_file = 'out.thermal.0'
r = am.load('atom_dump', ref_file, symbols='Al')
ids = np.arange(0,len(r.atoms))
pos_r = r.atoms.prop(key='pos',index=ids)


r.atoms.disp = np.zeros([1,3])


atom_dump = r.dump('atom_dump')
disp_filename = 'out.thermal.disp.0'
f = open(disp_filename, "w")
f.write(atom_dump)
f.close()
print('Wrote: {}'.format(disp_filename))


for t in times:
    def_file = 'out.thermal.{}'.format(t)

    d = am.load('atom_dump', def_file, symbols='Al')

    pos_d = d.atoms.prop(key='pos',index=ids)

    del_pos = am.dvect(pos_d,pos_r,r.box,r.pbc)

    print(t)
    mags = np.linalg.norm(del_pos,axis=1)
    mmin = np.min(mags)
    mmax = np.max(mags)
    print('Mag: {:.3f} to {:.3f}'.format(mmin, mmax))

    zs = del_pos[:,2]
    zmin = np.min(zs)
    zmax = np.max(zs)
    print('Z_coord: {:.3f} to {:.3f}'.format(zmin, zmax))

    d.atoms.disp = np.zeros([1,3])
    d.atoms.prop('disp',ids,del_pos)

    atom_dump = d.dump('atom_dump')
    disp_filename = 'out.thermal.disp.{}'.format(t)
    f = open(disp_filename, "w")
    f.write(atom_dump)
    f.close()
    print('Wrote: {}'.format(disp_filename))
