"""Particle Data Class - A class that can store TLorentzVectors, PID, MomPID"""

# Python Class with extra attributes for data storage
# Subclass of ROOT.TLorentzVector
# Created by Nicholas Zube, 8.22.2014
# Tested in Python 2.6.4

from ROOT import TLorentzVector
from ROOT import TVector3

class particleData(TLorentzVector):

    def __init__(self):
        super(particleData, self).__init__()
        self.sPID = 0
        self.sMomPID = 0
        self.status = 0
        self.mcParentage = 0
        self.vertex = TVector3()

    def MomPID(self):
        return self.sMomPID
                
    def PID(self):
        return self.sPID
    
    def Status(self):
        return self.status

    def McParentage(self):
        return self.mcParentage
                
    def SetMomPID(self, iMomPID):
        self.sMomPID = iMomPID

    def SetPID(self, iPID):
        self.sPID = iPID
        
    def SetStatus(self, iStatus):
        self.status = iStatus

    def SetMcParentage(self, iMcParentage):
        self.mcParentage = iMcParentage

    def SetVertex(self, xx, yy,zz):
        self.vertex.SetXYZ(xx,yy,zz)
