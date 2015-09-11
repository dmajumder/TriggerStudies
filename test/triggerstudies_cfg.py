import FWCore.ParameterSet.Config as cms

process = cms.Process("Trig")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('cerr', 'cout'), 
#    cerr = cms.untracked.PSet(threshold = cms.untracked.string('ERROR') ), 
#    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "MCRUN2_74_V9" 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/data/Run2015B/SingleMuon/MINIAOD/17Jul2015-v1/30000/16B50792-172E-E511-B0C8-0025905C43EC.root'
      #'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/160C08A3-4227-E511-B829-02163E01259F.root'
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_1.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_2.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_3.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_4.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_5.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_6.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_7.root', 
      #'root://cmsxrootd.fnal.gov///store/user/devdatta/TprimeBToTH_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring15DR74-Asympt25ns-HLT_users_drankin_FullMenuForBOOST_V2_MCRUN2_74_V9_MINIAODSIM-v1/150714_123545/0000/MINIAOD_8.root', 
      )
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      "triggers_data_Run2015B.root" 
      #"singleTprime_triggerStudies_22July2015.root" 
      )
    )

process.trig = cms.EDAnalyzer('TriggerStudies',
    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("selectedPatTrigger"),
    prescales = cms.InputTag("patTrigger"), 
    origpath = cms.string(""),
    hltProcName = cms.string("HLT"), 
    ak8jetLabel  = cms.InputTag('slimmedJetsAK8'),
    ak4jetLabel  = cms.InputTag('slimmedJets'),
    hltPaths = cms.vstring(
      "HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_v2", 
      "HLT_AK8PFHT500_TrimR0p1PT0p03Mass50_BTagCSV0p45_v2", 
      "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v2", 
      "HLT_PFHT800_v1", 
      "HLT_AK8DiPFJet200_200_TrimMass30_BTagCSV0p45_v2", 
      "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v2", 
      "HLT_AK8PFJet360_TrimMass30_v2", 
      ), 
    triggertypes = cms.vint32(89, 89, 89, 89, 85, 85, 85),# look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    newthreshs = cms.vdouble(650, 600, 700, 800, 250, 280, 360),
    )

process.p = cms.Path(process.trig)

