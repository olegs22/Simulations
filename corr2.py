import numpy as np
import ctypes
import os

_file = 'corr_2p.so'
_path = os.path.join(*(os.path.split(__file__)[:-1] + (_file,)))
_mod = ctypes.cdll.LoadLibrary(_path)

class DoubleArrayType:
    def from_param(self, param):
        typename = type(param).__name__
        if hasattr(self, 'from_' + typename):
            return getattr(self, 'from_' + typename)(param)
        elif isinstance(param, ctypes.Array):
            return param
        else:
            raise TypeError("Can't convert %s" % typename)

    def from_ndarray(self, param):
        return param.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

DoubleArray = DoubleArrayType()
_corr_2p = _mod.corr_2p
_corr_2p.argtypes = (ctypes.c_int, DoubleArray, DoubleArray, DoubleArray, ctypes.c_int, ctypes.c_int, ctypes.c_int)
_corr_2p.restype = ctypes.POINTER(ctypes.c_double)

def corr2(data, box_size, bins, model):
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    result = np.zeros(bins)

    c_function = _corr_2p(len(x), x, y, z, box_size, bins, model)

    for i in range(bins):
        result[i] = c_function[i]
    _mod.free(c_function)
    return result

_corr_2p_norm =_mod.corr_2p_norm
_corr_2p_norm.argtypes = (ctypes.c_int, ctypes.c_int, DoubleArray, DoubleArray,DoubleArray,\
                        DoubleArray, DoubleArray, DoubleArray, ctypes.c_int, ctypes.c_int,\
                        ctypes.c_int)
_corr_2p_norm.restype = ctypes.POINTER(ctypes.c_double)

def corr_2p_n(positions, randoms, box_size, bins, model):
    p_x = positions[:,0]
    p_y = positions[:,1]
    p_z = positions[:,2]

    r_x = randoms[:,0]
    r_y = randoms[:,1]
    r_z = randoms[:,2]

    result = np.zeros(bins)

    function = _corr_2p_norm(len(p_x), len(r_x), p_x, p_y, p_z, r_x, r_y, r_z,\
                            box_size, bins, model)
    for i in range(bins):
        result[i] = function[i]
    _mod.free(function)    
    return result

def binning(box_size, bins):
    max_distance = np.sqrt(3.0) * box_size
    b= np.linspace(0.0, max_distance, bins)
    return b
