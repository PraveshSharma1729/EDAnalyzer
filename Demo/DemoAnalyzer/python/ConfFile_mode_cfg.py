import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = "124X_mcRun3_2023_realistic_v12"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_0.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_1.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_2.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_3.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_4.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_5.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_6.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_7.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_8.root', 'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_Interference/step1_9.root'
           ),
            skipBadFiles = cms.untracked.bool(True),
            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output/ZGJJtoNuNuGJJ_Interference.root")
)

process.demo = cms.EDAnalyzer('DemoAnalyzer',
    jets = cms.InputTag("ak4GenJets"),
    mets = cms.InputTag("genMetTrue"),
    photons = cms.InputTag("genParticles")
)

process.p = cms.Path(process.demo)
