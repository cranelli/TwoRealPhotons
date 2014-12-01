from ROOT import gSystem
from ROOT import TFile
from ROOT import TTree
from ROOT import TDirectory, gDirectory
from ROOT import TLorentzVector
from ROOT import TH1F, TH2F
from ROOT import TBranchElement
from ROOT import TVector3
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
inRootFileDir = '../NtupleRootFiles/'
inRootFileName = 'GeneratorsTree_Wg_test.root'
treeLoc = "demo/gens_tree"
outRootFileName = 'lheTruePhotonHistograms_Wg.root'
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
        if entry%10000 == 0: print entry
        analysis_tree.GetEntry(entry)
        # print analysis_tree.__dict__
            
        # Select Event, should be only one per entry.

        # Assign Particles
        # Will be filled with TLorentz Vector, for now.
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
        
        lhe_candidate_photons = []
        lhe_candidate_electrons = []
        lhe_candidate_muons = []
        lhe_candidate_nu_es = []
        lhe_candidate_nu_ms = []
        lhe_candidate_ws = []
        
        mc_candidate_photons = []
        mc_candidate_electrons = []
        mc_candidate_muons = []
        mc_candidate_nu_es = []
        mc_candidate_nu_ms = []
        mc_candidate_ws = []
        
        prompt_photons = []
        
        #
        # Particle Identification
        #

        particleIdentification.assignParticles(analysis_tree, mc_photons, mc_electrons, mc_muons,
                                               mc_nu_es, mc_nu_ms, mc_ws)
                                               
        particleIdentification.assignLHEParticles(analysis_tree, lhe_photons, lhe_electrons, lhe_muons,
                                                  lhe_nu_es, lhe_nu_ms, lhe_ws)

        # histogramBuilder.fillNumParticleHistograms(mc_photons, 'MCPhotonMultiplicity_PreCuts')
        # histogramBuilder.fillNumParticleHistograms(lhe_photons, 'LHEPhotonMultiplicity_PreCuts')

        passCandidateCuts(lhe_candidate_photons, lhe_candidate_electrons, lhe_candidate_muons,
                                lhe_photons, lhe_electrons, lhe_muons)

        if not passCandidateCuts(mc_candidate_photons, mc_candidate_electrons, mc_candidate_muons,
                                mc_photons, mc_electrons, mc_muons): continue
        
        # Select True Photons based on Parentage
        promptPhotons = filter(lambda photon: (photon.McParentage() & non_prompt_mask) != non_prompt_mask,
                             mc_candidate_photons)

        #print "Num True Photons: ", len(promptPhotons)
        # photons_PtLT5 = filter(lambda photon: photon.Pt() < 5, photons)
    
        # Electron Channel
        if len(mc_candidate_electrons) == 1:
            #histogramBuilder.fillNumParticleHistograms(mc_candidate_photons,     'MCPhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(lhe_candidate_photons,     'LHEPhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(promptPhotons, 'PromptPhotonMultiplicity_ElectronChannel_MCPtdRCuts')
            
            if len(promptPhotons) > 1:
                print "More than One Prompt Photon - Electron Channel"
                histogramBuilder.fillDeltaRHistograms(promptPhotons, lhe_photons, "dRPromptLHEPhotons_PromptPhotonnGT1_ElectronChannel_MCPtdRCuts")
                new_photons = makeNewPhotons(lhe_photons,promptPhotons)
                histogramBuilder.fillParentHistograms(new_photons, "MomPID_NewPhoton_ElectronChannel")
                # histogramBuilder.fillStandardHistograms(new_photons, mc_electrons, mc_muons, "NewPhotons_ElectronChannel_MCPtdRCuts")
                histogramBuilder.fillVertexSeparationHistograms(new_photons,mc_electrons, "vertexSeparation_NewPhotonvMC_ElectronChannel")
                displacement=  TVector3()
                for new_photon in new_photons:
                    for mc_electron in mc_electrons:
                        displacement = new_photon.Vertex()-mc_electron.Vertex()
                        print "Separation ", displacement.Mag()
        
        # Muon Channel
        if len(mc_candidate_muons) == 1:
            #histogramBuilder.fillNumParticleHistograms(mc_candidate_photons,     'MCPhotonMultiplicity_MuonChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(lhe_candidate_photons, 'LHEPhotonMultiplicity_MuonChannel_MCPtdRCuts')
            histogramBuilder.fillNumParticleHistograms(promptPhotons, 'PromptPhotonMultiplicity_MuonChannel_MCPtdRCuts')

            if len(promptPhotons) > 1:
                print "More than One True Photon --Muon Channel"
                histogramBuilder.fillDeltaRHistograms(promptPhotons, lhe_photons, 'dR_PromptLHEPhotons_PromptPhotonGT1_MuonChannel_MCPtdRCuts')
                new_photons = makeNewPhotons(lhe_photons, promptPhotons)
                histogramBuilder.fillParentHistograms(new_photons, "MomPID_NewPhotons_MuonChannel")
                # histogramBuilder.fillStandardHistograms(new_photons, mc_electrons, mc_muons, 'PromptPhotonGT1_MuonChannel_MCPtdRCuts')

                displacement=  TVector3()
                for new_photon in new_photons:
                    for mc_muon in mc_muons:
                        displacement = new_photon.Vertex()-mc_muon.Vertex()
                        print "Separation ", displacement.Mag()
        

    
        
        # Candidate Cuts on the LHE Events
        # if passCandidateCuts(lhe_candidate_photons, lhe_candidate_electrons, lhe_candidate_muons,
        #                    lhe_photons, lhe_electrons, lhe_muons):
        #    histogramBuilder.fillStandardHistograms(lhe_candidate_photons, lhe_candidate_electrons, lhe_candidate_muons, "LHEPtdRCuts")
                
       
        #photons_Boson = filter(lambda photon: (photon.McParentage() & boson_parent_mask) == boson_parent_mask, photons)
                         
    outRootFile = TFile(outRootFileName, 'RECREATE')
    outRootFile.cd()
    
    #Iterate Over All Histograms
    for key, Histogram in histogramBuilder.Histograms.iteritems():
        Histogram.Write()
    
    inRootFile.Close()
    outRootFile.Close()

# Defines the Candidate Photons, Electrons, and Muons (python passes by reference)
# Returns True if it passes selection criteria, false otherwise.
# Currently Cuts are Made on Photon and Lepton Pt as well as their dR.

def passCandidateCuts(candidate_photons, candidate_electrons, candidate_muons,
                     photons, electrons, muons):

    #print "Filtered Photons", 
    #print filter(lambda photon: photon.Pt() > 15, photons)

    candidate_photons.extend(filter(lambda photon: photon.Pt() > 15, photons))
    candidate_electrons.extend( filter(lambda electron: electron.Pt() > 30, electrons))
    candidate_muons.extend(filter(lambda muon: muon.Pt > 25, muons))

    if not ((len(candidate_electrons)== 1 and len(candidate_muons)==0) or
            (len(candidate_electrons)== 0 and len(candidate_muons)==1)): return False
    
    # Photon Photon dR Cut
    if not eventCuts.passPhotonPhotonDeltaR(candidate_photons): return False
        
    # Electron Channel
    if len(candidate_electrons) == 1:
        if not eventCuts.passPhotonElectronDeltaR(candidate_photons, candidate_electrons): return False
        #print "working1"

    # Muon Channel
    if len(candidate_muons) == 1:
        if not eventCuts.passPhotonMuonDeltaR(candidate_photons, candidate_muons): return False

    return True

                                                            
# Select New Photons, by taking those that are not dR matched to an lhe photon.
def makeNewPhotons(lhe_photons, mc_photons):
    for lhe_photon in lhe_photons:
        new_photons = filter(lambda mc_photon: mc_photon.DeltaR(lhe_photon) > 0.4, mc_photons)
    return new_photons



if __name__=="__main__":
    lheTruePhotonSelection()




