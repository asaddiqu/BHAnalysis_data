import FWCore.ParameterSet.Config as cms 

process = cms.Process('ANA')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag.globaltag = 'PHYS14_25_V1::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
        fileName=cms.string("ntuple_output.root"),
        closeFileFast = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(True),
)

process.out = cms.OutputModule('PoolOutputModule',                                                                                                           fileName = cms.untracked.string('edm_output.root'),                                                                              
        outputCommands = cms.untracked.vstring('keep *')
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/*.root*')
    #fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/mc/Phys14DR/QCD_HT_1000ToInf_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/C6923693-3274-E411-B984-002481E94C7A.root')
    fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FE26BEB8-D575-E411-A13E-00266CF2AE10.root')
)

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules_ph = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules_ph:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.bhana = cms.EDAnalyzer('BHAnalyzerTLBSM',
  beamSpot = cms.InputTag('offlineBeamSpot'),
  electronTag = cms.InputTag("slimmedElectrons"),
  muonTag = cms.untracked.InputTag("slimmedMuons"),
  jetTag = cms.untracked.InputTag("slimmedJets"),
  tauTag = cms.untracked.InputTag("slimmedTaus"),
  metTag = cms.untracked.InputTag("slimmedMETs"),
  photonTag = cms.InputTag("slimmedPhotons"),
  rho_lable    = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
  ebRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEBRecHits"),
  eeRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEERecHits"),
  primaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  triggerTag = cms.untracked.InputTag("TriggerResults","","HLT"),
  
  verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
  conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

  eleVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
  eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
  eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
 
  phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
  phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
 
  MCLabel = cms.untracked.bool(True),                               
  DEBUG = cms.untracked.bool(False)                               
)


process.p = cms.Path((process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) * process.bhana)

