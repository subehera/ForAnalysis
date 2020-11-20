from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'gen_sim_coherentJpsi_starlight'
config.General.workArea = 'Nov13'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '1st_STEP_GEN_SIM.py'
config.Data.outputPrimaryDataset = 'gen_sim_coherentJpsi_starlight_MC'
config.Data.userInputFiles = open('input.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 2000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

#config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_heavyions/subehera/gen_sim_Stalight_Jpsi_Coherent/'
config.Data.allowNonValidInputDataset = True
config.Data.publication = True
config.Data.outputDatasetTag = 'gen_sim_starlight_coherentJpsi'
config.Site.storageSite = 'T2_CH_CERN'
