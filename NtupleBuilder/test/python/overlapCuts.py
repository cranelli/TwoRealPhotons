"""Overlap Cuts - Makes DeltaR Selection Cuts for LNuAA Analysis"""

# Python Module For Selection Cuts
# Created by Christopher Anelli, 8.4.2014
# Edits by Nicholas Zube, 8.21.2014

from ROOT import TLorentzVector
from ROOT import TH1F, TH2F

import histogramBuilder

# Delta R Cuts are Part of the Event Cuts

# Removes Photon if it is too close to an electron
def selectOnPhotonElectronDeltaR(photons, electrons):
    for electron in electrons:
        photons = filter(lambda photon:
                         photon.DeltaR(electron) > minPhotonElectronDeltaR, photons)
    histogramBuilder.fillDeltaRHistograms(photons, electrons, 'PhotonElectron_DeltaR_PostPhotonElectronDeltaRCut')
    return photons

def selectOnPhotonMuonDeltaR(photons, muons):
    for muon in muons:
        photons = filter(lambda photon:
                         photon.DeltaR(muon) > minPhotonMuonDeltaR, photons)
    histogramBuilder.fillDeltaRHistograms(photons, muons, 'PhotonMuon_DeltaR_PostPhotonMuonDeltaRCut')
    return photons

# Removes Photon if it is too close to another 
def selectOnPhotonPhotonDeltaR(photons):
    for photon1 in photons:
        photons = filter(lambda photon2:
                         photon1.DeltaR(photon2) > minPhotonPhotonDeltaR
                         or photon2.Pt() > photon1.Pt() # Remove Lower P Photon 
                         or photon1 == photon2, # Do Not Remove itself 
                         photons)
    # Currently Double Counts
    histogramBuilder.fillDeltaRHistograms(photons, photons, 'PhotonPhoton_DeltaR_PostPhotonPhotonDeltaRCut')
    return photons
        