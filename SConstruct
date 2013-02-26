import glob

env = Environment(CCFLAGS = '-DdTRIMESH_ENABLED -DdDOUBLE -DGRAPHICS -g ')
env.AppendENVPath('CPLUS_INCLUDE_PATH', '/mnt/home/olsonran/biped-hyperneat/ode-0.11.1/include/')

current=['biped.cpp'] #glob.glob('*.cpp')
simplega=glob.glob('simplega/*.cpp')
rtneat=glob.glob('rtneat/*.cpp')

allsrc=current+simplega+rtneat

env.Program('biped', allsrc,LIBS=['m','SM','ICE','pthread','GL','GLU','drawstuff','ode'],LIBPATH=['.','/usr/lib/','/usr/local/lib'])

