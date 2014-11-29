"""Histogram Builder - fills a dictionary of histograms"""

# Python Module For Histogram Filling
# Created by Christopher Anelli, 8.4.2014
# Tested in Python 2.6.4

from ROOT import TLorentzVector
from ROOT import TH1F, TH2F

Histograms = {}
 
# Function to Create and Fill Count Histograms
def fillCountHistograms(particles, key, bins=2, xmin=0, xmax=2):
    if not key in Histograms:
        Histograms[key] = TH1F(key,key, bins, xmin, xmax)
    for particle in particles: Histograms[key].Fill(1)

#Function to Create and Fill Number of Particles Histograms
def fillNumParticleHistograms(particles, key, bins=20, xmin=0, xmax=20):
    if not key in Histograms:
        Histograms[key] = TH1F(key,key, bins, xmin, xmax)
    Histograms[key].Fill(len(particles))

# Function to Create and Fill Pt Histograms
def fillPtHistograms(particles, key, bins=400, xmin=0, xmax=200):
    if not key in Histograms:
        Histograms[key] = TH1F(key,key, bins, xmin, xmax)
    for particle in particles: Histograms[key].Fill(particle.Pt())

# Function To Create and Fill Eta Histograms
def fillEtaHistograms(particles, key, bins=120, xmin=-3, xmax=3):
    if not key in Histograms:
        Histograms[key] = TH1F(key,key, bins, xmin, xmax)
    for particle in particles: Histograms[key].Fill(particle.Eta())

# Watch out for double counting when particles1 = particles2
def fillDeltaRHistograms(particles1, particles2, key, bins=200, xmin=0, xmax=10):
    if not key in Histograms:
        Histograms[key] = TH1F(key, key, bins, xmin, xmax)
    for particle1 in particles1:
        for particle2 in particles2:
            if particle1 != particle2: # Don't Histogram comparison with itself
                Histograms[key].Fill(particle1.DeltaR(particle2))
        
def fillStatusHistograms(particles, key, bins=8, xmin=0, xmax=4):
    if not key in Histograms:
        Histograms[key] = TH1F(key,key, bins, xmin, xmax)
    for particle in particles: Histograms[key].Fill(particle.Status())

def fillMHistograms(m, key, bins=300, xmin=0, xmax=300):
    if not key in Histograms:
        Histograms[key] = TH1F(key, key, bins, xmin, xmax)
    Histograms[key].Fill(m)
    
def fillStandardHistograms(photons, electrons, muons, suffix):
    # Photons
    fillCountHistograms(photons, 'Photon_'+suffix)
    fillPtHistograms(photons, 'Photon_Pt_' + suffix)
    fillEtaHistograms(photons, 'Photon_Eta_' + suffix)
    fillNumParticleHistograms(photons, 'Photon_Num_' + suffix)
    fillStatusHistograms(photons, "Photon_Status_" + suffix)
    # Electrons
    fillCountHistograms(electrons, 'Electron_' + suffix)
    fillPtHistograms(electrons, 'Electron_Pt_' + suffix)
    fillEtaHistograms(electrons, 'Electron_Eta_' + suffix)
    fillNumParticleHistograms(electrons, 'Electron_Num_' + suffix)
    fillStatusHistograms(electrons, "Photon_Status_" + suffix)
    # Muons
    fillCountHistograms(muons, 'Muon_' + suffix)
    fillPtHistograms(muons, 'Muon_Pt_' + suffix)
    fillEtaHistograms(muons, 'Muon_Eta_' + suffix)
    fillNumParticleHistograms(muons, 'Muon_Num_' + suffix)
    fillStatusHistograms(muons, "Photon_Status_" + suffix)
    # Delta R
    fillDeltaRHistograms(photons, electrons, 'DeltaR(Ae)_' + suffix)
    fillDeltaRHistograms(photons, muons, 'DeltaR(AMu)_' + suffix)
    fillDeltaRHistograms(photons, photons, 'DeltaR(AA)_' + suffix)
    

