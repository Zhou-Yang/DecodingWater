"""
Yang Zhou

to extract waters at a given coodinate (X,Y,Z) around 1.2 A
"""

import mdtraj as md
import numpy as np

TRAJ = 'md_noPBC.xtc'
TP = 'npt.gro'
RADIUS = 0.14 #Unit is nm
XYZ=np.array([2.4538, 4.1146, 3.6754], dtype='float32')


def distance(XYZ1=np.array([0, 0, 0], dtype='float32'),
             XYZ2=np.array([1, 1, 1], dtype='float32')):
    """the distance between XYZ1 and XYZ2"""
    a=XYZ2-XYZ1
    b=a**2
    c=b.sum()
    return np.sqrt(c)

def indice_in_Radius(t,
                     water_indice,
                     r=RADIUS):
    """return the indice of water in a radius"""
        
    for i in water_indice:
        wat_indice=[]
        if distance(t.xyz[0,i,:],XYZ) < RADIUS:
            wat_indice.append(i)
            wat_indice.append(i+1)
            wat_indice.append(i+2)
            print "HOH"
    return wat_indice    

def main():
    print "Loading topology %s and traj %s" %(TP,TRAJ)
    t = md.load(TRAJ, top = TP)
    water_oxygen_indice = t.top.select('water and name O')

    print "There are %d waters and %d frames" %(water_oxygen_indice.size,t.n_frames)
    
    f=0
    count=0
    while(f < t.n_frames):
        print "left %s frame" %(t.n_frames-f)
        water_oxygen_indice = t[f].top.select('water and name O')
        wat_indice = []
        for i in water_oxygen_indice:
            if distance(t.xyz[f,i,:],XYZ) < RADIUS:
                wat_indice.append(i)
                wat_indice.append(i+1)
                wat_indice.append(i+2)
                count += 1
                print "HOH: %s,%s" %(distance(t.xyz[f,i,:],XYZ),count)
        if wat_indice:
            wateringrid=t[f].atom_slice(wat_indice)
            wateringrid.save_pdb('w%3.2f_%3.2f_%3.2f_%s.pdb' %(XYZ[0],XYZ[1],XYZ[2],f))
        f+=1
        print "\n"
    print "XYZ: %s, occupationcy: %s" %(XYZ,count/t.n_frames)

if __name__ == "__main__":
    main()
