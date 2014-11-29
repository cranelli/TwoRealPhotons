"""Selection Cuts - Makes Selection Cuts for LNuAA Analysis"""

# Python Module For Selection Cuts
# Created by Christopher Anelli, 8.4.2014
# Tested in Python 2.6.4

from ROOT import TLorentzVector
from ROOT import TH1F, TH2F

import histogramBuilder

# Selection Cut Values
# Photons
minPhotonPt = 15
maxPhotonEta = 2.5
minPhotonEndCapEta = 1.47
maxPhotonEndCapEta = 1.57
# Electrons
minElectronPt = 30
maxElectronEta = 2.5
# Muons
minMuonPt = 25
maxMuonEta = 2.4

def selectOnPhotonKinematics(photons):
    #if photons != None:
    # Post Cut Histograms
    photons = filter(lambda photon: photon.Pt() > minPhotonPt, photons)
    histogramBuilder.fillPtHistograms(photons, 'Photon_Pt_PostPhotonPtCut')
    photons = filter(lambda photon: abs(photon.Eta()) < maxPhotonEta, photons)
    #photons = filter(lambda photon: abs(photon.Eta()) < minPhotonEndCapEta
    #                 or abs(photon.Eta()) > maxPhotonEndCapEta, photons)
    histogramBuilder.fillEtaHistograms(photons, 'Photon_Eta_PostPhotonEtaCut')
    return photons

def selectOnElectronKinematics(electrons):
    #if electrons != None:
    electrons = filter(lambda electron: electron.Pt() > minElectronPt, electrons)
    histogramBuilder.fillPtHistograms(electrons, 'Electron_Pt_PostElectronPtCut')
    electrons = filter(lambda electron: abs(electron.Eta()) < maxElectronEta, electrons)
    histogramBuilder.fillEtaHistograms(electrons, 'Electron_Eta_PostElectronEtaCut')
    return electrons

def selectOnMuonKinematics(muons):
    muons = filter(lambda muon: muon.Pt() > minMuonPt, muons)
    histogramBuilder.fillPtHistograms(muons, 'Muon_Pt_PostMuonPtCut')
    muons = filter(lambda muon: abs(muon.Eta()) < maxMuonEta, muons)
    histogramBuilder.fillEtaHistograms(muons, 'Muon_Eta_PostMuonEtaCut')
    return muons
