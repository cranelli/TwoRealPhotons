"""Parent Cuts - Makes Parentage Selection Cuts for LNuAA Analysis"""

# Python Module For Parentage Cuts
# Created by Nicholas Zube, 8.22.2014
# Tested in Python 2.6.4

from ROOT import TLorentzVector
#from ROOT import TH1F, TH2F

import histogramBuilder

def selectOnPhotonParent(photons):
    #if photons != None:
    filter_photons = filter(lambda photon: abs(photon.MomPID()) < 25, photons)
    histogramBuilder.fillPtHistograms(filter_photons, 'Photon_Pt_PostParentCut')
    histogramBuilder.fillEtaHistograms(filter_photons, 'Photon_Eta_PostParentCut')
    return filter_photons
    
def selectOnElectronParent(electrons):
    #if electrons != None:
    electrons = filter(lambda electron: abs(electron.MomPID()) == 24, electrons)
    histogramBuilder.fillPtHistograms(electrons, 'Electron_Pt_PostWParentCut')
    histogramBuilder.fillEtaHistograms(electrons, 'Electron_Eta_PostWParentCut')
    return electrons

def selectOnMuonParent(muons):
    #if muons != None:
    muons = filter(lambda muon: abs(muon.MomPID()) == 24, muons)
    histogramBuilder.fillPtHistograms(muons, 'Muon_Pt_PostWParentCut')
    histogramBuilder.fillEtaHistograms(muons, 'Muon_Eta_PostWParentCut')
    return muons
