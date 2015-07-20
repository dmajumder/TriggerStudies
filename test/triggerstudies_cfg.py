import FWCore.ParameterSet.Config as cms

process = cms.Process("Trig")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "MCRUN2_74_V9" 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-5000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'root://cmsxrootd.fnal.gov///store/data/Run2015B/MuonEG/RAW/v1/000/251/162/00000/20F6A4B7-AD25-E511-B30F-02163E01358B.root' 
      #'file:/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_7_4_6_patch6/src/HLTrigger/Configuration/test/GENSIMRAWHLT.root' 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_1.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_2.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_3.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_4.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_5.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_6.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_7.root', 
      '/store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_8.root', 
      )
    )

process.hltPFHT800emu0 = cms.EDFilter("MiniAODTrigEmu",
    bits = cms.InputTag("TriggerResults","","MYHLT"),
    objects = cms.InputTag("selectedPatTrigger"),
    origpath = cms.string("HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_BTagCSV0p45_v2"),# original path to use as a base
    newthresh = cms.double(600),# new threshold to use
    triggertype = cms.int32(89),# look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    )

process.hltPFHT800emu1 = cms.EDFilter("MiniAODTrigEmu",
    bits = cms.InputTag("TriggerResults","","MYHLT"),
    objects = cms.InputTag("selectedPatTrigger"),
    origpath = cms.string("HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_v2"),# original path to use as a base
    newthresh = cms.double(650),# new threshold to use
    triggertype = cms.int32(89),# look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    )

process.hltPFHT800emu2 = cms.EDFilter("MiniAODTrigEmu",
    bits = cms.InputTag("TriggerResults","","MYHLT"),
    objects = cms.InputTag("selectedPatTrigger"),
    origpath = cms.string("HLT_AK8DiPFJet200_200_TrimMass30_BTagCSV0p45_v2"),# original path to use as a base
    newthresh = cms.double(250),# new threshold to use
    triggertype = cms.int32(85),# look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    )

process.trig = cms.EDAnalyzer('TriggerStudies',
    bits = cms.InputTag("TriggerResults","","MYHLT"),
    prescales = cms.InputTag("patTrigger"), 
    hltProcName = cms.string("MYHLT"), 
    ak8jetLabel  = cms.InputTag('slimmedJetsAK8'),
    ak4jetLabel  = cms.InputTag('slimmedJets'),
    hltPaths = cms.vstring(
      "HLT_PFHT800_v1", 
      "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v2", 
      "HLT_AK8PFJet360_TrimMass30_v2", 
      "HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_v2", # threshold changed to 650 GeV 
      "HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_BTagCSV0p45_v2", # threshold changed to 600 GeV 
      "HLT_AK8DiPFJet200_200_TrimMass30_BTagCSV0p45_v2", # threshold changed to 250 GeV 
      )
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      "singleTprime_triggerStudies_20July2015.root" 
      )
    )

process.p0 = cms.Path(process.hltPFHT800emu0*process.trig) 
process.p1 = cms.Path(process.hltPFHT800emu1*process.trig) 
process.p2 = cms.Path(process.hltPFHT800emu2*process.trig) 

process.schedule = cms.Schedule(process.p0, process.p1, process.p2) 

