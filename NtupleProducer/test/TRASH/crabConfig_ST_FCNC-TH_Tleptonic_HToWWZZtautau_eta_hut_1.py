from WMCore.Configuration import Configuration
config = Configuration()

argv = ['/eos/cms/store/caf/user/mdjordje/Cirkovic/FCNC6', 'ST_FCNC-TH_Tleptonic_HToWWZZtautau_eta_hut', 'output_1.root']

config.section_("General")
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'pset.py'
config.JobType.scriptExe = 'batch_job_ny_crab.sh'
#config.JobType.scriptArgs = ['Nevt='+argv[1], 's='+argv[2], 'c='+argv[3], 'deltaEtaCut='+argv[4], 'energy='+argv[5]]
config.JobType.scriptArgs = [ argv[0]+'/'+argv[1]+'/'+argv[2], argv[0]+'/'+argv[1]+'_Ntuple/'+argv[2]]
config.JobType.inputFiles = ['batch_job_ny_crab.sh', 'NtupleProducer', 'table_ny.txt']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = [ argv[2] ]
config.JobType.sendPythonFolder = True
#config.JobType.numCores = 8
config.JobType.maxMemoryMB = 3000

config.section_("Data")
config.Data.publication = False
config.Data.splitting = 'EventBased'
#config.Data.totalUnits = int(argv[0])
config.Data.totalUnits = 1
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T2_HU_Budapest'

