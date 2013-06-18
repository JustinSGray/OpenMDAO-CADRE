from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_GSposEarth(Component):
    
    lon = Float(0, iotype="in")
    lat = Float(0, iotype="in")
    alt = Float(0, iotype="in")
    
    def __init__(self, n):
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        self.add('r_e2g_E', Array(iotype='out', shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computejacobiangs(self.lon, self.lat, self.alt)
        self.dr_dlon, self.dr_dlat, self.dr_dalt = result

    def execute(self):
        self.r_e2g_E = self.lib.computegs(self.n, self.lon, self.lat, self.alt)

    def applyDer(self, arg, result):
        if 'lon' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] = self.dr_dlon[k] * arg['lon']
        if 'lat' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dlat[k] * arg['lat']
        if 'alt' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dalt[k] * arg['alt']
        return result

    def applyDerT(self, arg, result):
        if 'r_e2g_E' in arg:
            for k in xrange(3):
                result['lon'] = self.dr_dlon[k] * numpy.sum(arg['r_e2g_E'][k,:])
                result['lat'] = self.dr_dlat[k] * numpy.sum(arg['r_e2g_E'][k,:])
                result['alt'] = self.dr_dalt[k] * numpy.sum(arg['r_e2g_E'][k,:])
        return result
