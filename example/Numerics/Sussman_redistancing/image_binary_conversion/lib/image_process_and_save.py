import os
import numpy as np

def padding(img, x_pad, y_pad, z_pad):
    z, y, x = img.shape[0], img.shape[1], img.shape[2]
    zp, yp, xp = z + 2*z_pad, y + 2*y_pad, x + 2*x_pad
    x_padded = np.zeros((zp, yp, xp ))
    x_padded[z_pad:zp-z_pad, y_pad:yp-y_pad, x_pad:xp-x_pad] = img
    return x_padded
    
def save_zslices_toBin(numpy_array, filename):
    import os
    curdir = os.getcwd()
    directory = '%s/%s'% (curdir, filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    # loop through z dim and save every slice as separate binary file bc. savetxt only works for 1d or 2d arrays
    z_dim = numpy_array.shape[0]
    for z in range(z_dim):
        numpy_array[z].astype('int8').tofile('%s/zslice_%d.bin' % (directory, z))

def save_array_toBin(numpy_array, path, filename):
    if not os.path.exists(path):
        os.makedirs(path)
    numpy_array.astype('int8').tofile('%s/%s.bin' % (path, filename))

