from phenom.waveform.waveform import Waveform

class PhenomDInternals(object):
    """docstring for PhenomDInternals"""
    def __init__(self, arg):
        self.arg = arg

class PhenomD(PhenomDInternals, Waveform):
    """docstring for PhenomD"""
    def __init__(self, **p):
        self.p = p
    def amp(self):
        pass
