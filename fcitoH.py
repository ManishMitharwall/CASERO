#!/storage/praha1/home/manishkumar/.conda/envs/kumarpython3/bin/python3.8

import sys

def FCIDUMP2Hfile(s,u,v,w,z):
    FCIDUMP = [s,u,v,w,z]  
    Hfile_line = list()
    if w==z and z==0:
        a = FCIDUMP
        Hfile_line += [(a[0],a[2],a[3],a[4],a[1])]
        if u ==v and v== w and w==z and u==0:
            Hfile_line = []
    else:
        a = FCIDUMP
        b = [(a[0],a[1],a[3],a[4],a[2]),(a[0],a[2],a[3],a[4],a[1]),(a[0],a[1],a[4],a[3],a[2]),(a[0],a[2],a[4],a[3],a[1])]
        c = sorted(set(b), key=b.index)
        toremove = []
        for k in c:
            allsign = [(k[0],k[1],k[2],k[3],k[4]),(k[0],k[1],-1*k[2],-1*k[3],k[4]),(k[0],-1*k[1],k[2],k[3],-1*k[4]),(k[0],-1*k[1],-1*k[2],-1*k[3],-1*k[4])]            
            for l in allsign:    
                if l[-1]!=l[-2] and l[-3]!=l[-4] and l not in toremove:
                    Hfile_line += [l]
                toremove += [l]
                toremove += [(l[0],l[2],l[1],l[4],l[3])]
    return Hfile_line

# print(FCIDUMP2Hfile(0.4,2,1,2,1))
Hfile = list()
with open(sys.argv[1]) as f1:
    f1.readline();f1.readline();f1.readline();f1.readline()
    for line in f1:
        split = line.split()
        c=format(float(split[0]), '.8f')
        coeff =27.2114*float(c) ;d = int(split[1]);e=int(split[2]);f=int(split[3]);g=int(split[4])
        Hfile_line= FCIDUMP2Hfile(coeff,d,e,f,g)
        Hfile += Hfile_line
with open('Hfile.txt','w') as f2:
    for j in Hfile:
        f2.write('%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n ' % (j[0], j[1], j[2],j[3],j[4]))
