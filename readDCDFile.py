import numpy
import sys
sys.path.append('/home/zy/anaconda2/lib/python2.7/site-packages/')
from MDAnalysis import Universe, Writer
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis

DCD='md_noPBC.xtc'
PSF='npt.gro'

distanceMat=open('distance.txt','w')
rho=Universe(PSF,DCD)
##print rho
##print list(rho.residues)
p=rho.select_atoms('protein and not backbone and not(name H*)')
w=rho.select_atoms('resname SOL and not(name H*)')
##print list(p)
#pc=p.coordinates()
pc=p.positions
print len(pc)
proteResid=p.resids
waterResid=w.resids
proteResnu=p.resnames
waterResna=w.resnames
waterResnu=w.resnums
atomInf=[]
for i in w.atoms:
    atomid= str(i).split()[2]
    atomseg=str(i).split()[-1]
    atomidandseg=[]
    atomidandseg.append(atomid)
    atomidandseg.append(atomseg)
    atomInf.append(atomidandseg)
print len(atomInf)
##print len(waterResnu)



   
##print list(w)

box = rho.trajectory.ts.dimensions[:3]

frameNo=1
for ts in rho.trajectory:
    distanceMat.write('frame'+' '+str(frameNo)+'\n')
    title='resname'+'    '+'atomid'+'    '+'resnumber'+'    X    Y     Z   '+'   '+'segname'+'\n'
    distanceMat.write(title)
    wc=w.positions
##    print str(wc[0])[1:-1]
##    print list(wc)[0]
##    print pc
    d=distance_array(pc,wc,box)
##    print len(d)
##    print len(d[0])
    print ts
    waterNP=[]
    for i in d:
##        print len(i)
        waterid=0
        
        for j in i:
##            print j
            if j>6:
                waterid+=1
            elif j<=6:
                water=str(waterResna[waterid])+'    '+atomInf[waterid][0]+'    '+str(waterResnu[waterid])+'    '+str(wc[waterid])[1:-1]+'   '+atomInf[waterid][1]+'\n'
                if water not in waterNP:
                    waterNP.append(water)
                
                waterid+=1
                break
    for waterN in waterNP:
        distanceMat.write(waterN)
##                distanceMat.write(str(waterid)+'    '+str(j)+'\n')


##            print waterid
    frameNo+=1
                
##                
##    for i in d:
####        print d[1],'@@@@@@@@@@@@@@@@@@@@@@@@@'
##        line=line+str(i)+'  '
##    line=line+'\n'
##    distanceMat.write(line)

distanceMat.close()
    
##    for wac in wc:
##        for pac in pc:
##            dist=numpy.linalg.norm(wac-pac)
##            print dist
##            if dist<3:
##                break
##    
            
##    print box
##    for wac in wc:
##    d=distance_array(wc,pc,box)
##    a.append(d)
####    len(a)
##    print ts
##    for distance in a:
##    distance=str(distance)[1:-2]+'\n'
##    distanceMat.write(distance)
##
##distanceMat.close()
    
    
    
##        for pac in pc:
##            a=numpy.array(wac)
##            b=numpy.array(pac)

            
            
    



##bb=u.select_atoms('protein and backbone')
##print list(bb)
##print bb.coordinates()

