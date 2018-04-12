#
#
import os 
import sys
import shutil

shutil.copy("../src/test","../bin/test")
os.system("mpiexec -n 3 ../bin/test")
