import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string('GeneratorsTree.root')
                                   )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data/users/cranelli/WGamGam/AODSIM/ZGToLLG_AODSIM_test.root'
    )
)

process.demo = cms.EDAnalyzer('NtupleBuilder')


process.p = cms.Path(process.demo)
