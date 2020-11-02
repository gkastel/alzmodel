import numpy as np

mm = np.random.rand(100, 10000);

bb = (mm>0.95).astype(int)

print(np.sum(bb, axis=1)/10.)

np.savetxt('inspikes.txt', bb, '%d', ' ')
