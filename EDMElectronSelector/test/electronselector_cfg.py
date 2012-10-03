import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/user/b2g12006/isandova/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12-PU_S6_START52_V9-v1_TLBSM_52x_v5/2bcc6fdd1e664d93e9026c3764d0b403/ttbsm_52x_mc_77_0_3iY.root''
    )
)

process.myProducerLabel = cms.EDProducer('EDMElectronSelector'
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
