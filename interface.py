#!/usr/bin/env python2

class PDBAnalyzer(object):
    ''' This class should not be modified directly, except for adding
        more methods to the interface. Instead, this class should
        be subclassed and the methods overriden with an actual implementation'''

    def shortest_contacts(self, model, chainColors, stickColors):
        raise NotImplementedError()
    
