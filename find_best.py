import os
import sys
fn=sys.argv[1]
#fn="gen33"
os.system("grep Novelty %s > out.txt" % fn) 
a=open("out.txt").read().split("\n")[:-1]
b=[float(k.split()[4]) for k in a]
num=b.index(max(b))
print max(b)
raw_input("?")
os.system("python extract.py %s %d" % (fn,num))

