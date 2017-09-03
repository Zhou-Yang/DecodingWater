"""
Yang Zhou

Identify high occupationcy water site
"""

import mdtraj as md
import numpy as np

IN_PDB = 'all_water.pdb'
GRID_SIZE = 0.14
CUT_OFF = 10
SPHERE_RIDIUS = 0.14
CUT_OFF2 = 100

def distance(XYZ1=np.array([0, 0, 0], dtype='float32'), 
             XYZ2=np.array([1, 1, 1], dtype='float32')):
    """the distance between XYZ1 and XYZ2"""
    a=XYZ:2-XYZ1
    b=a**2
    c=b.sum()
    return np.sqrt(c)


def XYZ_count(hist,binedges,cutoff=10):
    """
    count the number of waters in a given XYZ position its above cutoff
    """
    Xgridsize=len(binedges[0])-1
    Ygridsize=len(binedges[1])-1
    Zgridsize=len(binedges[2])-1
    X=[]
    Y=[]
    Z=[]
    count=[]
    for i in range(Xgridsize):
        for j in range(Ygridsize):
            for k in range(Zgridsize):
                if hist[i][j][k] > cutoff:
                    X.append(binedges[0][i])
                    Y.append(binedges[1][j])
                    Z.append(binedges[2][k])
                    count.append(hist[i][j][k])
    return X,Y,Z,count


def main():

    print "Loading pdb file ..."
    
    t = md.load(IN_PDB)
    points=t.xyz[0]
    waterO = t.topology.select('name O')
    size_v = GRID_SIZE**3
    size = GRID_SIZE

    print "pdb file:%s, grid size:%s nm^3" %(IN_PDB,size_v)
     
    
    #size of the big box
    xmax,xmin=t.xyz[0,waterO,0].max(),t.xyz[0,waterO,0].min()
    ymax,ymin=t.xyz[0,waterO,1].max(),t.xyz[0,waterO,1].min()
    zmax,zmin=t.xyz[0,waterO,2].max(),t.xyz[0,waterO,2].min()
    sizeofthebox=(xmax-xmin)*(ymax-ymin)*(zmax-zmin)

    #bin number    
    xbins=int((xmax-xmin)/size)
    ybins=int((ymax-ymin)/size)
    zbins=int((zmax-zmin)/size)
    
    #distribute XYZ coordinates in pdb
    hist1,binedges1=np.histogramdd(points, bins= (xbins,ybins,zbins) ,normed=False)
    
    ####    
    coordlistX,coordlistY,coordlistZ,numlist=XYZ_count(hist1,binedges1,CUT_OFF)    
    
    ####
    atom_id=[i for i,x in enumerate(numlist) if x >= CUT_OFF]

    atomX=[]
    atomY=[]
    atomZ=[]
    ncount=[]
    
    for i in atom_id:
        atomX.append(coordlistX[i])
        atomY.append(coordlistY[i])
        atomZ.append(coordlistZ[i])
        ncount.append(numlist[i])
    
    s = SPHERE_RIDIUS

    hsx=[]
    hsy=[]
    hsz=[]
    hsn=[]
    
    for i in range(len(ncount)):
        print "identifying occupationcy above %s" %(CUT_OFF2)
        XYZ=np.array([atomX[i],atomY[i],atomZ[i]], dtype='float32')
        x=[]
        y=[]
        z=[]
        n=0
        wat_indice=[]
        for j in waterO:
            XYZj=t.xyz[0,j,:]
            if (distance(XYZj,XYZ)<=s):
                x.append(t.xyz[0,j,0])
                y.append(t.xyz[0,j,1])
                z.append(t.xyz[0,j,2])
                n+=1
                wat_indice.append(j)
                wat_indice.append(j+1)
                wat_indice.append(j+2)
        if n >= CUT_OFF2:
            #wateringrid=t.atom_slice(wat_indice)
            #wateringrid.save_pdb('water_site_%s.pdb' %(i))
            hsx.append(sum(x)/len(x))
            hsy.append(sum(y)/len(y))
            hsz.append(sum(z)/len(z))
            hsn.append(n)
    
    for i in range(len(hsx)):
        print ("%s,%s,%s,%s") %(hsx[i],hsy[i],hsz[i],hsn[i])




if __name__ == "__main__":
    main()



