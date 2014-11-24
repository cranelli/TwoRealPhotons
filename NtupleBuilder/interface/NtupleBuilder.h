#ifndef __HOMUON_HOMUONTREEBUILDER_H__
#define __HOMUON_HOMUONTREEBUILDER_H__


// system include files                                                                                            
#include <memory>

// user include files                                                                                              
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TTree.h"
#include "vector"
#include "math.h"

using namespace::std;

/*
Description: Ntuple Builder to get the LHE and Gen Information from our background samples.
*/

//
// class declaration                                                                
//  

class NtupleBuilder : public edm::EDAnalyzer {
 public:
  explicit NtupleBuilder(const edm::ParameterSet&);
  ~NtupleBuilder();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void SetBranches();
  void ClearVectors();

  edm::Service<TFileService> _fileService;
  TTree *gens_tree;         

  int nLHE;
  vector<int> lhePID;
  vector<float> lhePx;
  vector<float> lhePy;
  vector<float> lhePz;
  vector<float> lheE;
  vector<float> lheM;

  int nMC_;
  vector<int> mcPID;
  vector<float> mcVtx_x;
  vector<float> mcVtx_y;
  vector<float> mcVtx_z;
  vector<float> mcPt;
  vector<float>  mcMass;
  vector<float>  mcEta;
  vector<float>  mcPhi;
  vector<float>  mcE;
  vector<float>  mcEt;
  vector<int>    mcGMomPID;
  vector<int>    mcMomPID;
  vector<float>  mcMomPt;
  vector<float>  mcMomMass;
  vector<float>  mcMomEta;
  vector<float>  mcMomPhi;
  vector<int>    mcIndex;
  vector<int>    mcDecayType;
  vector<int>    mcParentage;
  vector<int>    mcStatus;

  // ----------member data ---------------------------                                                         
};

// 
// constants, enums and typedefs
//                                                                                                                 

//                                                                                    
// static data member definitions                                                  
//



#endif
