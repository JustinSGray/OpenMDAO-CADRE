from openmdao.lib.datatypes.api import Float, Array
from openmdao.main.api import Component
import numpy as np

class Orbit_Initial(Component):
	altPerigee = Float(0, iotype="in", copy=None)
	altApogee = Float(0, iotype="in", copy=None)
	RAAN = Float(0, iotype="in", copy=None)
	Inc = Float(0, iotype="in", copy=None)
	argPerigee = Float(0, iotype="in", copy=None)
	trueAnomaly = Float(0, iotype="in", copy=None)
	
	def __init__(self):
		super(Orbit_Initial, self).__init__()
		self.add('r_e2b_I0', Array(np.ones(3)), size=3, dtype=np.float, iotype='out')
		self.add('v_e2b_I0', Array(np.ones(3)), , size=3, dtype=np.float,iotype='out')
		
	def compute(self):
		Re=6378.137
		mu=398600.44
		
		def S(v):
			S = np.zeros((3,3),complex)
			S[0,:] = [0, -v[2], v[1]]
			S[1,:] = [v[2], 0, -v[0]]
			S[2,:] = [-v[1], v[0], 0]
			return S
			
		def getRotation(axis,angle):
			R = np.eye(3,dtype=complex) + S(axis)*np.sin(angle) + (1 - np.cos(angle)) * (np.outer(axis,axis) - np.eye(3,dtype=complex))
			return R
			
	d2r = np.pi/180.0
	r_perigee = Re + altPerigee
	r_apogee = Re + altApogee
	e = (r_apogee-r_perigee)/(r_apogee+r_perigee)
	a = (r_perigee+r_apogee)/2
	p = a*(1-e**2)
	h = np.sqrt(p*mu)

	rmag0 = p/(1+e*np.cos(d2r*trueAnomaly))
	r0_P = np.array([rmag0*np.cos(d2r*trueAnomaly), rmag0*np.sin(d2r*trueAnomaly), 0], complex)
	v0_P = np.array([-np.sqrt(mu/p)*np.sin(d2r*trueAnomaly), np.sqrt(mu/p)*(e+np.cos(d2r*trueAnomaly)), 0], complex)

	O_IP = np.eye(3, dtype=complex)
	O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),RAAN*d2r))
	O_IP = np.dot(O_IP, getRotation(np.array([1,0,0]),Inc*d2r))
	O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),argPerigee*d2r))

	r0_ECI = np.dot(O_IP,r0_P)
	v0_ECI = np.dot(O_IP,v0_P)

	return r0_ECI, v0_ECI
		
	def linearize(self):
		h = 1e-16
		ih = complex(0,h)
		v = np.zeros(6,complex)
		v[:] = [self.altPerigee[0]. self.altApogee[0], slef.RAAN[0], self.Inc[0], self.argPerigee[0], self.trueAnomaly[0])
		self.J = np.zeros((6,6))
		for i in range(6):
			v[i] += ih
			r0_ECI, v0_ECI = self.compute(v[0], v[1], v[2], v[3], v[4], v[5])
			v[i] -= ih
			self.J[:3,i] = r0_ECI.imag/h
			self.J[3:,i] = v0_ECI.imag/h
	
	def execute(self):
		r0_ECI, v0_ECI = self.compute(self.altPerigee[0], self.altApogee[0], self.RAAN[0], self.Inc[0], self.argPerigee[0], self.trueAnomaly[0])
		self.r_e2b_I0 = r0_ECI.real
		self.v_e2b_I0 = v0_ECI.real
		
	def applyDer(salf, arg, result):
		if not result['r_e2b_I0']:
			result['r_e2b_I0'] = np.zeros(3)
		if not result['v_e2b_I0']:
			result['v_e2b_I0'] = np.zeros(3)
		J = self.J
		result['r_e2b_I0'] = J[:3,0]*arg['altPerigee'][0] + J[:3,1]*arg['altApogee'][0] + J[:3,2]*arg['RAAN'][0] + J[:3,3]*arg['Inc'][0] + J[:3,4]*arg['argPerigee'][0] + J[:3,5]*arg['trueAnomaly'][0]
		result['v_e2b_I0'] = J[3:,0]*arg['altPerigee'][0] + J[3:,1]*arg['altApogee'][0] + J[3:,2]*arg['RAAN'][0] + J[3:,3]*arg['Inc'][0] + J[3:,4]*arg['argPerigee'][0] + J[3:,5]*arg['trueAnomaly'][0]
	
		return result
		
	def applyDerT(self, arg, result):
		J = self.J
		result['altPerigee'][0] += sum(J[:3,0]*arg['r_e2b_I0'][:]) + sum(J[3:,0]*arg['v_e2b_I0'][:])
		result['altApogee'][0] += sum(J[:3,1]*arg['r_e2b_I0'][:]) + sum(J[3:,1]*arg['v_e2b_I0'][:])
		result['RAAN'][0] += sum(J[:3,2]*arg['r_e2b_I0'][:]) + sum(J[3:,2]*arg['v_e2b_I0'][:])
		result['Inc'][0] += sum(J[:3,3]*arg['r_e2b_I0'][:]) + sum(J[3:,3]*arg['v_e2b_I0'][:])
		result['argPerigee'][0] += sum(J[:3,4]*arg['r_e2b_I0'][:]) + sum(J[3:,4]*arg['v_e2b_I0'][:])
		result['trueAnomaly'][0]+= sum(J[:3,5]*arg['r_e2b_I0'][:]) + sum(J[3:,5]*arg['v_e2b_I0'][:])
