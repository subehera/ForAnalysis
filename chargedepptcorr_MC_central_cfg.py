import FWCore.ParameterSet.Config as cms

process = cms.Process("CmsTwoPartCorrAnalysis")


# __________________ General _________________                                                                                                

# Configure the logger                                                                                                                        
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Configure the number of maximum event the analyser run on in interactive mode                                                               
# -1 == ALL                                                                                                                                   
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(-1)                                                                                                          
    input = cms.untracked.int32(200)
    )


# __________________ I/O files _________________                                                                                              

# Define the input file to run on in interactive mode                                                                                         
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/himc/HINPbPbWinter16DR/Hydjet_Quenched_MinBias_5020GeV_750/AODSIM/NoPU_75X_mcRun2_HeavyIon_v13\
_75X_mcRun2_HeavyIon_v13-v1/80000/001E4607-5BBA-E611-9A99-0CC47A7E6A2C.root'                                                                  


        #'/store/himc/HINPbPbAutumn18DR/MinBias_Hydjet_Drum5F_2018_5p02TeV/AODSIM/NoPUmva98_103X_upgrade2018_realistic_HI_v11-v1/270000/F87DB\
3E2-4AF9-4D40-9B91-707A101B1399.root'                                                                                                         
        '/store/user/caber/PixelTracks3pletsPlus4pletsTightHitCleaning_MCHydjet_v1/MinBias_Hydjet_Drum5F_2018_5p02TeV/crab_PixelTracks3pletsP\
lus4pletsTightHitCleaning_MCHydjet_v1/201121_234129/0000/out_MC_1.root'

    )
)

# Define output file name                                                                                                                     
import os
process.TFileService = cms.Service("TFileService",
     #fileName = cms.string(os.getenv('CMSSW_BASE') + '/src/Analyzers/ChargeDepAndPtCorr/test/chargeptdepcorr.root')                          
     fileName = cms.string('chargeptdepcorr_MC_central.root')
)


# __________________ Detector conditions _________________                                                                                    
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')                                                    
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic_hi', '')
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag                                                                                
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
   cms.PSet(record = cms.string("HeavyIonRcd"),
            #tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x02_mc"),                                                        
            tag = cms.string("CentralityTable_HFtowers200_HydjetCymbal5Ev8_v1020x01_mc"),
      connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
      label = cms.untracked.string("HFtowers")
   ),
])



# __________________ Event selection _________________                                                                                        
# Load centrality producer for centrality calculation                                                                                         
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.newCentralityBin = process.centralityBin.clone()

# Load HI event selection modules                                                                                                             
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
# __________________ Analyse Sequence _________________                                                                                       

# Load you analyzer with initial configuration                                                                                                
process.load("Analyzers.CmsTwoPartCorrAnalysis.chargedepptcorr_cff")
process.defaultAnalysis_05   = process.CPDC05.clone()
#process.defaultAnalysis_510   = process.CPDC510.clone()                                                                                      
process.p = cms.Path(#process.hfCoincFilter3 *             # Requier HF coincidence with 3 GeV                                                
                     process.primaryVertexFilter *        # Clean up on vertices                                                              
                     #process.clusterCompatibilityFilter * # Clean up on pileup                                                               
                     process.centralityBin *              # Compute centrality                                                                
                     process.defaultAnalysis_05)



