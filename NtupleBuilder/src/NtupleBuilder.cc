#include "TwoRealPhotons/NtupleBuilder/interface/NtupleBuilder.h"
#include "TwoRealPhotons/NtupleBuilder/interface/GenParticleParentage.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// -*- C++ -*-
//
// Package:    NtupleBuilder
// Class:      NtupleBuilder
// 
/**\class NtupleBuilder NtupleBuilder.cc Analysis/NtupleBuilder/src/NtupleBuilder.cc

 Implementation:
     [Notes on implementation]
*/
//
// Original Author: Christopher  Anelli
//         Created:  Thu Nov 20 20:14:14 EST 2014

//
// constructors and destructor
//
NtupleBuilder::NtupleBuilder(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  gens_tree = _fileService->make<TTree>("gens_tree", 
					"LHE and General Generator Info");

  SetBranches();

}


NtupleBuilder::~NtupleBuilder()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleBuilder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   ClearVectors();
   nLHE = 0;
   nMC_ = 0;
   
   // Get the LHE Info
   Handle<LHEEventProduct> lhe_event;
   iEvent.getByLabel("source",lhe_event);
   
   const lhef::HEPEUP hepeup = lhe_event->hepeup();
   nLHE  = hepeup.NUP;
   lhePID = hepeup.IDUP;
   
   vector<lhef::HEPEUP::FiveVector> pup = hepeup.PUP;
   for( unsigned int i =0; i < pup.size(); i++){
     lhePx.push_back(pup[i][0]);
     lhePy.push_back(pup[i][1]);
     lhePz.push_back(pup[i][2]);
     lheE.push_back(pup[i][3]);
     lheM.push_back(pup[i][4]);
   }
   
   // Get Generator Info (Following Same Steps as in the ggNtuplizer)
   
   Handle<std::vector<reco::GenParticle> > genParticlesHandle;
   iEvent.getByLabel("genParticles", genParticlesHandle);
   
   
   if(genParticlesHandle.isValid()) {
     //cout << "Gen Particle Handle is Valid" << endl;
     
     int genIndex = 0;
     
     for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin();
	  ip != genParticlesHandle->end(); ++ip) {
       genIndex++;
       int status = ip->status() - 10*(ip->status()/10);
       bool stableFinalStateParticle = status == 1 && ip->pt() > 5.0;
       
       // keep all the photons with pT > 5.0 and all leptons;                                                                              
       bool photonOrLepton =
	 (status == 1 && ip->pdgId() == 22 && ip->pt() > 5.0 ) ||
	 (status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ))  ||
	 (status < 10 && abs(ip->pdgId()) == 15 );
	 
       // select also Z, W, H, and top                                                                                                     
       bool heavyParticle =
	 (ip->pdgId() == 23 || abs(ip->pdgId()) == 24 || ip->pdgId() == 25 ||
	  abs(ip->pdgId()) == 6 || abs(ip->pdgId()) == 5);
       
       if ( stableFinalStateParticle || heavyParticle || photonOrLepton ) {
	 const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	 //if (!runOnParticleGun_ && !p->mother()) continue;
	 if (!p->mother()) continue;  //removed ParticleGun Case, it is ggNtulizer
	 reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle,
							     ip-genParticlesHandle->begin());
	 genpartparentage::GenParticleParentage particleHistory(partRef);

	 mcPID.push_back(p->pdgId());
	 mcVtx_x.push_back(p->vx());
	 mcVtx_y.push_back(p->vy());
	 mcVtx_z.push_back(p->vz());
	 mcPt.push_back(p->pt());
	 mcMass.push_back(p->mass());
	 mcEta.push_back(p->eta());
	 mcPhi.push_back(p->phi());
	 mcE.push_back(p->energy());
	 mcEt.push_back(p->et());

	 
	 mcParentage.push_back(
			       particleHistory.hasLeptonParent()*16   +
			       particleHistory.hasBosonParent()*8     +
			       particleHistory.hasNonPromptParent()*4 +
			       particleHistory.hasQCDParent()*2       +
			       particleHistory.hasExoticParent());
	 mcStatus.push_back(p->status());
	 
	 int mcGMomPID_ = -999;
	 int mcMomPID_  = -999;
	 float mcMomPt_    = -999.;
	 float mcMomMass_  = -999.;
	 float mcMomEta_   = -999.;
	 float mcMomPhi_   = -999.;
	 if ( particleHistory.hasRealParent() ) {
	   reco::GenParticleRef momRef = particleHistory.parent();
	   if ( momRef.isNonnull() && momRef.isAvailable() ) {
	     mcMomPID_  = momRef->pdgId();
	     mcMomPt_   = momRef->pt();
	     mcMomMass_ = momRef->mass();
	     mcMomEta_  = momRef->eta();
	     mcMomPhi_  = momRef->phi();

	     // get Granny                                                                                                         
	     genpartparentage::GenParticleParentage motherParticle(momRef);
	     if ( motherParticle.hasRealParent() ) {
	       reco::GenParticleRef granny = motherParticle.parent();
	       mcGMomPID_ = granny->pdgId();
	     }
	   }
	 }
	 
	 mcGMomPID.push_back(mcGMomPID_);
	 mcMomPID.push_back(mcMomPID_);
	 mcMomPt.push_back(mcMomPt_);
	 mcMomMass.push_back(mcMomMass_);
	 mcMomEta.push_back(mcMomEta_);
	 mcMomPhi.push_back(mcMomPhi_);

	 mcIndex.push_back(genIndex-1);
	 nMC_++; 
       }
     }
   }
   //   ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
   
   gens_tree->Fill();
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
NtupleBuilder::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleBuilder::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
NtupleBuilder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
NtupleBuilder::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
NtupleBuilder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
NtupleBuilder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleBuilder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void NtupleBuilder::SetBranches(){
  gens_tree->Branch("nLHE", &nLHE);
  gens_tree->Branch("lhePID",
                       "std::vector<int>", &lhePID);
  gens_tree->Branch("lhePx",
                       "std::vector<float>", &lhePx);
  gens_tree->Branch("lhePy",
                       "std::vector<float>", &lhePy);
  gens_tree->Branch("lhePz",
                       "std::vector<float>", &lhePz);
  gens_tree->Branch("lheE",
                       "std::vector<float>", &lheE);
  gens_tree->Branch("lheM",
                       "std::vector<float>", &lheM);
  
  gens_tree->Branch("nMC", &nMC_);
  gens_tree->Branch("mcPID", &mcPID);
  gens_tree->Branch("mcVtx_x", &mcVtx_x);
  gens_tree->Branch("mcVtx_y", &mcVtx_y);
  gens_tree->Branch("mcVtx_z", &mcVtx_z);
  gens_tree->Branch("mcPt", &mcPt);
  gens_tree->Branch("mcMass", &mcMass);
  gens_tree->Branch("mcEta", &mcEta);
  gens_tree->Branch("mcPhi", &mcPhi);
  gens_tree->Branch("mcE", &mcE);
  gens_tree->Branch("mcEt", &mcEt);
  gens_tree->Branch("mcGMomPID", &mcGMomPID);
  gens_tree->Branch("mcMomPID", &mcMomPID);
  gens_tree->Branch("mcMomPt", &mcMomPt);
  gens_tree->Branch("mcMomMass", &mcMomMass);
  gens_tree->Branch("mcMomEta", &mcMomEta);
  gens_tree->Branch("mcMomPhi", &mcMomPhi);
  gens_tree->Branch("mcIndex", &mcIndex);
  gens_tree->Branch("mcDecayType", &mcDecayType); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau 
  gens_tree->Branch("mcParentage", &mcParentage); // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
  gens_tree->Branch("mcStatus", &mcStatus); // status of the particle
}

/*
 * Function to Clear the Vectors
 */
void NtupleBuilder::ClearVectors(){
   lhePID.clear();
   lhePx.clear();
   lhePy.clear();
   lhePz.clear();
   lheE.clear();
   lheM.clear();

   mcPID.clear();
   mcVtx_x.clear();
   mcVtx_y.clear();
   mcVtx_z.clear();
   mcPt.clear();
   mcMass.clear();
   mcEta.clear();
   mcPhi.clear();
   mcE.clear();
   mcEt.clear();
   mcGMomPID.clear();
   mcMomPID.clear();
   mcMomPt.clear();
   mcMomMass.clear();
   mcMomEta.clear();
   mcMomPhi.clear();
   mcIndex.clear();
   mcDecayType.clear();
   mcParentage.clear();
   mcStatus.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleBuilder);
