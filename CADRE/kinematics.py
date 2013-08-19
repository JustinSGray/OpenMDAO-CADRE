import numpy as np

def modulo(a, p):
    return  a - int(a / p) * p

def fixangles(n, azimuth0, elevation0):
    azimuth, elevation = np.zeros(azimuth0.shape), np.zeros(elevation0.shape)
    for i in xrange(n):
        azimuth[i] = modulo(azimuth0[i], 2*np.pi)
        elevation[i] = modulo(elevation0[i], 2*np.pi)
        if (elevation[i] > np.pi):
            elevation[i] = 2*np.pi - elevation[i]
            azimuth[i] = np.pi + azimuth[i]
            azimuth[i] = modulo(azimuth[i], 2*np.pi)
    return azimuth, elevation

def computepositionrotd(n, vects, mat):
    result = np.empty(vects.shape)
    for i in xrange(n):
        result[:,i] = np.dot(mat[:,:,i], vects[:,i])
    return result

def computepositionrotdjacobian(n, v1, O_21):
    J1, J2 = np.zeros((n, 3, 3, 3)), np.zeros((n, 3, 3))
    eye = np.eye(3)
    for i in xrange(n):
        for k in xrange(3):
            for u in xrange(3):
                for v in xrange(3):
                    J1[i,k,u,v] = eye[k,u]*v1[v,i]
                J2[i,k,u] = O_21[k,u,i]
    return J1, J2