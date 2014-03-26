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
    )
)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames=cms.untracked.vstring( 
                                    #'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_withBeniCalib.root'
                                    'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/BeniCalib/GRIN_RAW2DIGI_RECO_BeniMod_all.root'
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_00644AB6.root',
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_04D268BB.root',
                                    #'file:/user/erdogan/samples/RAW_RECO/RunA/ZMU/reprocessed/SingleMu_Run2012A_ZMu22Jan2013_v1_RAW2DIGI_RECO_083E3CA8.root'
                                    #'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_10_2_JNQ.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_11_1_a71.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_1_1_yln.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_12_1_7zy.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_13_1_11L.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_14_1_Ept.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_15_1_dYE.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_16_1_KUZ.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_17_1_dp4.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_18_1_0PR.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_19_1_ywm.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_20_2_nCK.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_21_1_4Jr.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_2_1_YFU.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_22_2_PJu.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_23_1_Txo.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_25_1_6Wj.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_26_1_9Ia.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_27_2_j8q.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_28_1_8mI.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_29_2_D1f.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_30_1_8Rt.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_31_1_XnI.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_33_1_Or7.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_35_2_EgS.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_37_1_EkA.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_38_1_krx.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_39_1_2Mg.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_40_2_pps.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_41_1_A8O.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_4_1_CNh.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_42_1_Xhv.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_43_1_vEG.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_45_1_3bF.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_46_1_wmk.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_48_1_Yh6.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_49_2_GDD.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_50_1_37V.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_5_1_A5S.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_59_1_hVm.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_60_2_g86.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_61_1_oTi.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_6_1_3tm.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_62_1_6KB.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_63_1_flP.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_64_1_Atg.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_65_1_fN0.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_66_1_EJy.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_67_1_L2c.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_68_1_SuN.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_7_1_K54.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_8_1_Db5.root'
                                    #,'file:/net/scratch_cms/institut_3b/erdogan/GRIN_RAW2DIGI_RECO/GRIN_RAW2DIGI_RECO_9_1_Xe3.root'
                                    #'file:/user/erdogan/workspace/mtt_ho_studies/CMSSW_5_3_14_patch2/src/rereco_2012data/step2_RAW2DIGI.root'
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
                              #beamSpot=cms.InputTag("offlineBeamSpot"),
                              #primaryVertex=cms.InputTag('offlinePrimaryVertices'),
                              TrackAssociatorParameters=TrackAssociatorParameterBlock.TrackAssociatorParameters
                              )


process.p = cms.Path(#process.noscraping
                     #* process.primaryVertexFilter
                      process.MuonSelector
                     * process.demo)

