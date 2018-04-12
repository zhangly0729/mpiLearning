#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from struct import unpack 
import sys
import json

nx=680
nz=400

blocksize=nx*nz
fmt = str(blocksize)+"f"

f=open("result.dat","rb")
V=unpack(fmt,f.read(blocksize*4))
f.close()

V=np.reshape(V,[nx,nz])
vmax=np.max(V)
vmin=np.min(V)

plt.figure(figsize=(8.1,6))

plt.pcolormesh(V.transpose(),cmap="jet")
plt.xlabel("x(i)")
plt.ylabel("y(i)")
plt.gca().set_aspect("equal")
plt.colorbar()
plt.savefig("xxx.png")
plt.show()




