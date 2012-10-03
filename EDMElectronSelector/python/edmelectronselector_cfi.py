import FWCore.ParameterSet.Config as cms

electronSelector = cms.EDProducer('EDMElectronSelector',
                                  electronSrc = cms.InputTag( "selectedPatElectronsPFlow" ),
                                  electronIsolationType = cms.string("relIsodb"),
                                  maximumElectronIsolation = cms.double(0.2)
                                  )
