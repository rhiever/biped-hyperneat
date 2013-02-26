Biped HyperNEAT project

Copyright 2013 Randal S. Olson, Joel Lehman, Kenneth Stanley.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

Dependencies:
-------------
Requires the OpenGL & GLU headers and libraries, 
as well as Open Dynamics Engine libraries and headers.

The version of Open Dynamics Engine I used was 0.11

The build tool I used was SCons, which should
work with windows, although the SConstruct
file will have to be modified. 


Using:
---------
To invoke the program for evolving:
./biped evolve <output path> <fitness/novelty> <cppn genes> <substrate genes> <novelty function (only for novelty)>
e.g.
fot fitness:
./biped evolve ./o/ fitness bipedstartgenes substrategenes

for novelty:
./biped evolve ./o/ novelty bipedstartgenes substrategenes 7

This will place output files in the ./o/ directory


To invoke the program for displaying:
./biped display <CPPN file name> <substrate genes file name> 0
e.g. ./biped display temp.dat substrategenes 0

To show the best from a population snapshot:
python find_best.py ./o/<file name>

This will display the first genome in the file
