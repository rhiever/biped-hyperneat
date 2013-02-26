import sys

fname=sys.argv[1] #"archive.dat"
num=int(sys.argv[2]) #20
a=open(fname).read()
a=a.split("/* Novelty")

out="temp2.dat"
out=open(out,"w")
if num==0:
   out.write(a[0])
else:
   towrite=a[num]
   f=towrite.find("*/")
   f=towrite.find("*/",f+1)
   out.write(towrite[f+3:])
out.close()

import os
os.system('grep -v \* temp2.dat > temp.dat') 
os.system("./biped display temp.dat substrategenes 0")
#os.system("gnuplot plotstuff.gnu -persist")
