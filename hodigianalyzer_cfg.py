import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#HBHE event-level noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'

#FileService for histograms
process.TFileService = cms.Service("TFileService",
    fileName=cms.string(
    'SingleMu_Run2012A_22Jan2013_v1_RECO_HODIGI.root'
    #'GRIN_analysis.root'
    )
)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames=cms.untracked.vstring(
                                    #'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_BeniMod_all.root',
                                    #'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO.root'
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_11_1_0y4.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_12_2_k3w.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_15_1_Dz2.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_16_1_wMo.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_25_2_bD3.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_28_1_eTO.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_33_1_d80.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_49_1_JRh.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_62_1_4Gn.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_64_2_9a2.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_65_1_unW.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_66_2_wfJ.root',
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_68_1_ybY.root'
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_00644AB6.root',
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_04D268BB.root',
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_083E3CA8.root'
                                    )
                            )

#Beam background removal
#process.noscraping = cms.EDFilter("FilterOutScraping",
#                                  applyfilter=cms.untracked.bool(True),
#                                  debugOn=cms.untracked.bool(False),
#                                  numtrack=cms.untracked.uint32(10),
#                                  thresh=cms.untracked.double(0.25)
#                                  )

#Primary vertex requirement
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection=cms.InputTag('offlinePrimaryVertices'),
#                                           minimumNDOF=cms.uint32(4) ,
#                                           maxAbsZ=cms.double(24),
#                                           maxd0=cms.double(2) 
#                                           )

process.MuonSelector = cms.EDFilter("MuonSelector",
                                     src=cms.InputTag("muons"),
                                     #cut=cms.string("isStandAloneMuon && eta < 0.8 && eta > -0.8 && pt >= 26")
                                     #cut=cms.string("isStandAloneMuon && pt >= 5")
                                     cut=cms.string("pt >= 5")
                                     )

#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#JSONfile = '/user/erdogan/workspace/mtt_ho_studies/CMSSW_5_3_14_patch2/src/rereco_2012data/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_FilterRun2012A.txt'
#myLumis = LumiList.LumiList(filename=JSONfile).getCMSSWString().split(',')
#process.source.lumisToProcess.extend(myLumis)

process.demo = cms.EDAnalyzer('HODIGIAnalyzer',
                              selection=cms.string(""),
 #                             beamSpot=cms.InputTag("offlineBeamSpot"),
 #                             primaryVertex=cms.InputTag('offlinePrimaryVertices'),
                              TrackAssociatorParameters=TrackAssociatorParameterBlock.TrackAssociatorParameters
                              )


process.p = cms.Path(#process.noscraping
                     #* process.primaryVertexFilter
                     process.MuonSelector
                     * process.demo)

