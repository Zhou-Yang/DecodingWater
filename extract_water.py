"""
Yang Zhou

To extract water molecules 0.5 nm around the protein from a traj
"""

import mdtraj as md
import numpy as np

TRAJ = 'md_noPBC_center.xtc'
TP = 'npt.gro'
S_FRAME = 0
OUT_NAME = "w"
AROUND_P = 0.5 # 0.5 nm around protein

def main():
    """main funtion to extract water"""
    print "extract water"
    #load file
    trj = md.load(TRAJ, top=TP)
    ref = md.load(TP)
    s_frame = S_FRAME

    print "traj file:%s, topology file:%s, outfile:%s" %(TRAJ, TP, OUT_NAME)

    #align protein
    proindice = ref.top.select('protein')
    trj = trj.superpose(ref, 0,
                        atom_indices=proindice,
                        ref_atom_indices=proindice,
                        parallel=True)

    #extract water around the center
    waterindice = trj.top.select('water')
    twater = trj.atom_slice(waterindice)

    all_water_o = trj.top.select('water and name O')


    while s_frame < trj.n_frames:
        print "extracting %s" %(s_frame)
        wat_indice = []
        water_o = md.compute_neighbors(trj[s_frame],
                                       AROUND_P,
                                       proindice,
                                       all_water_o)
        wat_indice = []
        for i in water_o[0]:
            wat_indice.append(i)
            wat_indice.append(i+1)
            wat_indice.append(i+2)
        wateringrid = trj[s_frame].atom_slice(wat_indice)
        wateringrid.save_pdb('%s_%s.pdb' %(OUT_NAME, s_frame))
        s_frame += 1

if __name__ == "__main__":
    main()
