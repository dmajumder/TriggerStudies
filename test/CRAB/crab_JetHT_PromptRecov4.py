from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'JetHT_PromptRecov4'
config.General.workArea = 'TrigEff_13Jan2016/'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'triggerstudies_cfg.py' 
config.JobType.pyCfgParams = []

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v4/MINIAOD' 
config.Data.inputDBS = 'global'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt' 
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/devdatta/TrigEff_13Jan2016/'
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'

config.section_('User')

