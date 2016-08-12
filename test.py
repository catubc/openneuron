import numpy as np


filename = '/media/cat/12TB/in_vivo/tim/yuki/IA1/tif_files/IA1am_Mar10_30Hz/IA1am_Mar10_30Hz_aligned.npy'

data = np.load(filename,  mmap_mode='r+')

print data.shape
