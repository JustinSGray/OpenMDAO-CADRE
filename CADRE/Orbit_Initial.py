from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np

class Orbit_Initial(Companent):
	altPerigee = Float(0, iotype="in", copy=None)
	altApogee = Float(0, iotype="in", copy=None)
	RAAN = Float(0, iotype="in", copy=None)
	Inc = Float(0, iotype="in", copy=None)
	argPerigee = Float(0, iotype="in", copy=None)
	trueAnomaly = Float(0, iotype="in", copy=None)
	
    def __init__(self):