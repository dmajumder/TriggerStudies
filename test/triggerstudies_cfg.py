import FWCore.ParameterSet.Config as cms

process = cms.Process("Trig")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('cerr', 'cout'), 
#    cerr = cms.untracked.PSet(threshold = cms.untracked.string('ERROR') ), 
#    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/0067D1EA-EE6F-E511-B561-0050560207C5.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/00729FDB-706F-E511-8E12-0050560207C5.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/00EE27AF-B16F-E511-A5F4-00259073E382.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/00F77227-DB6F-E511-8221-00259073E3D0.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/00FFECBB-A16F-E511-91F0-00259073E390.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/021DD636-B46F-E511-B192-00259073E382.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/02C56D74-8C6F-E511-8989-00505602078D.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/02D961B7-E16F-E511-B069-00259073E456.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/040E688D-A66F-E511-AF47-00259073E3FA.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/042532FE-0D70-E511-BB65-002590574604.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/0477FC28-DB6F-E511-B47A-00505602077B.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/048B088C-A36F-E511-A203-00259073E34A.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/048B9885-0570-E511-9BBD-00259073E4C2.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/0490F6B5-C26F-E511-8663-00259073E2F2.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/049843BA-E16F-E511-9605-00505602077B.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/04C695B3-D16F-E511-BDF3-00259073E3AC.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/04DAB0B2-F36F-E511-9836-00259073E4CA.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/0601378A-B86F-E511-B1FC-00259073E3CA.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/0659EC49-1370-E511-9370-00259073E34E.root',
      '/store/data/Run2015D/JetHT/MINIAOD/05Oct2015-v1/50000/06AFBAEF-D06F-E511-A081-0025907277BE.root',
      )
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
      "triggers_data_Run2015D_05Oct2015.root" 
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
      "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_v",
      "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v",
      "HLT_PFHT800_v",
      "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_v",
      "HLT_AK8PFJet360_TrimMass30_v",
      ), 
    triggertypes = cms.vint32(89, 89, 89, 85, 85),# look at https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    newthreshs = cms.vdouble(600, 650, 800, 250, 360),
    )

process.p = cms.Path(process.trig)
