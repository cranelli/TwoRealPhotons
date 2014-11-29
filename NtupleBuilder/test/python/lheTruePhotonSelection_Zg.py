from ROOT import gSystem
from ROOT import TFile
from ROOT import TTree
from ROOT import TDirectory, gDirectory
from ROOT import TLorentzVector
from ROOT import TH1F, TH2F
from ROOT import TBranchElement
from math import sqrt, cos
#import sys
#My Own Helper Modules
import particleIdentification
import objectCuts
import parentCuts
import overlapCuts
import eventCuts
import histogramBuilder

#
# Define Root Files and Tree Locations and Names
#

workDirLoc = '/home/cranelli/WGamGam/Anomolous_QGC/CMSSW_5_3_12/src/Anomolous_QGC/Analysis/python/'
inRootFileDir = '../'
inRootFileName = 'GeneratorsTree_Zg.root'
treeLoc = "demo/gens_tree"
outRootFileName = 'lheTruePhoton_Zg.root'
#print sys.path

#MC Parentage Masks
qcd_parent_mask = 2 #Binary Representation
non_prompt_mask = 4
boson_parent_mask = 8
lep_parent_mask = 16

#
# Begining of Analysis Code
#

def lheTruePhotonSelection():

    inRootFileLoc = inRootFileDir + inRootFileName
    inRootFile = TFile(inRootFileLoc, "READ")
    # Set Histograms
   
    analysis_tree = inRootFile.Get(treeLoc)
    print "Tree ", analysis_tree
    
    num_entries = analysis_tree.GetEntries()
    print "Number of Entries: ", num_entries
    
    #
    # Loop Over Entries
    #
    
    for entry in xrange(num_entries):
        analysis_tree.GetEntry(entry)
        # print analysis_tree.__dict__
            
        # Select Event, should be only one per entry.
        # event = analysis_tree.Event[0]
        # weight = event.Weight

        # Assign Particles
        # Will be filled with TLorentz Vector, for now.
        lhe_photons = []
        lhe_photons = []
        lhe_electrons = []
        lhe_muons = []
        lhe_nu_es = []
        lhe_nu_ms = []
        lhe_ws = []
        
        mc_photons = []
        mc_electrons = []
        mc_muons = []
        mc_nu_es = []
        mc_nu_ms = []
        mc_ws = []

        #
        # Particle Identification
        #

        particleIdentification.assignParticles(analysis_tree, mc_photons, mc_electrons, mc_muons,
                                               mc_nu_es, mc_nu_ms, mc_ws)
        particleIdentification.assignLHEParticles(analysis_tree, lhe_photons, lhe_electrons, lhe_muons,
                                                  lhe_nu_es, lhe_nu_ms, lhe_ws)

        histogramBuilder.fillNumParticleHistograms(mc_photons, 'MCPhotonMultiplicity_PreCuts')
        histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_PreCuts')
        
        
        # photons_PtLT5 = filter(lambda photon: photon.Pt() < 5, photons)
        # histogramBuilder.fillStandardHistograms(photons_PtLT5, electrons, muons, "PhotonPtLT5")

        #
        # Total (Fiducial) Cuts
        #

        # Events after Cutting on MC Photons
        mc_fiducialphotons = filter(lambda mc_photon: mc_photon.Pt() > 15, mc_photons)
        mc_fiducialelectrons = filter(lambda mc_electron: mc_electron.Pt() > 30, mc_electrons)
        mc_fiducialmuons = filter(lambda mc_muon: mc_muon.Pt > 25, mc_muons)
        #truePhotons = parentCuts.selectOnPhotonParent(mc_fiducialphotons)
        # Only Take Prompt Photons
        truePhotons = filter(lambda photon: (photon.McParentage() & non_prompt_mask) != non_prompt_mask, mc_fiducialphotons)

        if not eventCuts.passPhotonPhotonDeltaR(mc_fiducialphotons): continue
        #Separate Into Electron, Muon Channels
        if len(lhe_electrons) == 2:
            if not eventCuts.passPhotonElectronDeltaR(mc_fiducialphotons, mc_fiducialelectrons): continue
            histogramBuilder.fillNumParticleHistograms(mc_fiducialphotons, 'MCPhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(truePhotons, 'MCTruePhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            if len(truePhotons) > 1:
                histogramBuilder.fillDeltaRHistograms(truePhotons, lhe_photons, "PhotondR_MCLHE_MCTruePhotonnGT1_ElectronChannel_MCPtdRCuts")
                new_photons = truePhotons
                for lhe_photon in lhe_photons:
                    new_photons = filter(lambda new_photon: new_photon.DeltaR(lhe_photon) > 0.4, new_photons)
                histogramBuilder.fillStandardHistograms(new_photons, mc_electrons, mc_fiducialmuons, "NewPhotons_ElectronChannel_MCPtdRCuts")
            
        if len(lhe_muons) == 2:
            if not eventCuts.passPhotonMuonDeltaR(lhe_fiducialphotons, lhe_fiducialmuons): continue
            histogramBuilder.fillNumParticleHistograms(mc_photons, 'MCPhotonMultiplicity_MuonChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_MuonChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(truePhotons, 'MCMomPIDPhotonMultiplicity_MuonChannel_MCPtdRCuts')
            if len(truePhotons) > 1:
                histogramBuilder.fillDeltaRHistograms(truePhotons, lhe_photons, "PhotondR_MCLHE_MCMomPIDnPhoGT1_MuonChannel_MCPtdRCuts")
                new_photons = truePhotons
                for lhe_photon in lhe_photons:
                    new_photons = filter(lambda new_photon: new_photon.DeltaR(lhe_photon) > 0.4, new_photons)
                histogramBuilder.fillStandardHistograms(truePhotons, mc_electrons, mc_fiducialmuons, "MCMomPIDnPhoGT1_MuonChannel_MCPtdRCuts")
           
        
        # Events after Cutting on the LHE Events

        lhe_fiducialphotons = filter(lambda lhe_photon: lhe_photon.Pt() > 15, lhe_photons)
        lhe_fiducialelectrons = filter(lambda lhe_electron: lhe_electron.Pt() > 30, lhe_electrons)
        lhe_fiducialmuons = filter(lambda lhe_muon: lhe_muon.Pt > 25, lhe_muons)
        if not eventCuts.passPhotonPhotonDeltaR(lhe_fiducialphotons): continue

        #Separate Into Electron, Muon Channels
        if len(lhe_fiducialelectrons) == 2:
            if not eventCuts.passPhotonElectronDeltaR(lhe_fiducialphotons, lhe_fiducialelectrons): continue
            histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_ElectronChannel_LHEPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(mc_photons, 'MCPhotonMultiplicity_ElectronChannel_LHEPtdRCuts')

        if len(lhe_fiducialmuons) == 2:
            if not eventCuts.passPhotonMuonDeltaR(lhe_fiducialphotons, lhe_fiducialmuons): continue
            histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_MuonChannel_LHEPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(mc_photons, 'MCPhotonMultiplicity_MuonChannel_LHEPtdRCuts')
        
        #if len(lhe_photons) > 0: continue
        

        #Select on dR
        
        #photons_QCD = filter(lambda photon: (photon.McParentage() & qcd_parent_mask) == qcd_parent_mask, photons)
        #histogramBuilder.fillStandardHistograms(photons_QCD, electrons, muons, "PtCut_mcParentage_QCD_Cut")

        #photons_Boson = filter(lambda photon: (photon.McParentage() & boson_parent_mask) == boson_parent_mask, photons)
        #histogramBuilder.fillStandardHistograms(photons_Boson, electrons, muons, "PtCut_mcParentage_Boson_Cut")

        #photons_Lep = filter(lambda photon: (photon.McParentage() & lep_parent_mask) == lep_parent_mask, photons)
        #histogramBuilder.fillStandardHistograms(photons_Lep, electrons, muons, "PtCut_mcParentage_Lep_Cut")

        
                         
    outRootFile = TFile(outRootFileName, 'RECREATE')
    outRootFile.cd()
    
    #Iterate Over All Histograms
    for key, Histogram in histogramBuilder.Histograms.iteritems():
        Histogram.Write()
    inRootFile.Close()
    outRootFile.Close()

if __name__=="__main__":
    lheTruePhotonSelection()
    


