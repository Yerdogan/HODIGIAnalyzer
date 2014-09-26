// -*- C++ -*-
//
// Package:    HODIGIAnalyzer
// Class:      HODIGIAnalyzer
// 
/**\class HODIGIAnalyzer HODIGIAnalyzer.cc HODIGIAnalyzer/HODIGIAnalyzer/src/HODIGIAnalyzer.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Yusuf Erdogan
//         Created:  Wed Feb 19 15:38:21 CET 2014
// $Id$
//
//


// system include files
#include <memory>
#include <math.h>

// user include files
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/GeometryVector/interface/Point3DBase.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/RefVectorIterator.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/GeometryVector/interface/Point3DBase.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

// hcal calorimeter info
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/MuonIdentification/interface/MuonHOAcceptance.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"

//
// class declaration
//

class HODIGIAnalyzer: public edm::EDAnalyzer {
public:
	explicit HODIGIAnalyzer(const edm::ParameterSet&);
	~HODIGIAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	//reco::Vertex getvertex(const edm::Event&);
	//void getCalibHODigis(const edm::EventSetup&, HcalDetId, edm::SortedCollection<HODataFrame>::const_iterator, CaloSamples);
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);

	// ----------member data ---------------------------
	std::map<std::string, TH1*> histos_;
	std::map<std::string, TGraph*> graphs_;

	edm::Service<TFileService> theFileService;
	TFile* theOutputFile;

	StringCutObjectSelector<reco::Muon> selector_;

	TrackAssociatorParameters assocParams;
	TrackDetectorAssociator assoc;
	MuonHOAcceptance* theMuonHOAcceptance;

};

//
// constants, enums and typedefs
//

typedef edm::SortedCollection < HODataFrame > HODataFrameCollection;
double twopi = 2.*3.14159265358979323846;

//edm::InputTag beamspotLabel_;
//edm::InputTag primvertexLabel_;

//
// static data member definitions
//

static int const etaBounds = 3;
static int const phiSectorsR0 = 7;
static int const phiSectorsR12 = 2;

static double const etaMin[etaBounds] = { -0.3017, 0.3425, 0.8796 };
static double const etaMax[etaBounds] = { 0.3017, 0.8542, 1.2544 };

//static double const phiMinR0[phiSectors] = { -0.16172, 0.3618786, 0.8854773, 1.409076116, 1.932674892, 2.456273667, 2.979872443, 3.503471219, 4.027069994, 4.55066877, 5.074267545, 5.597866321 };
static double const phiMinR0[phiSectorsR0] = { -0.16172, 1.409076116, 1.932674892, 2.979872443, 3.503471219, 4.027069994, 5.597866321 };
//static double const phiMaxR0[phiSectors] = { 0.317395374, 0.84099415, 1.364592925, 1.888191701, 2.411790477, 2.935389252, 3.458988028, 3.982586803, 4.506185579, 5.029784355, 5.55338313, 6.076981906 };
static double const phiMaxR0[phiSectorsR0] = { 0.317395374, 1.888191701, 2.411790477, 3.458988028, 3.982586803, 4.506185579, 6.076981906 };

static double const phiMinR12[phiSectorsR12] = { 0.88093347, 1.404532245 };
static double const phiMaxR12[phiSectorsR12] = { 1.391186172, 1.914784947 };

//
// constructors and destructor
//
HODIGIAnalyzer::HODIGIAnalyzer(const edm::ParameterSet& iConfig) :
	selector_(iConfig.getParameter<std::string> ("selection"))

{
	//now do what ever initialization is needed

	assoc.useDefaultPropagator();

	assocParams = iConfig.getParameter<edm::ParameterSet> ("TrackAssociatorParameters");
	//beamspotLabel_ = iConfig.getParameter<edm::InputTag> ("beamSpot");
//	primvertexLabel_ = iConfig.getParameter<edm::InputTag> ("primaryVertex");

	//histos_["nevent_alive_after_cut"] = theFileService->make<TH1F> ( "nevent_alive_after_cut", "nevent_alive_after_cut", 7, -0.5, 6.5);

	histos_["n_muons"] = theFileService->make<TH1F> ("n_muons", "n_muons", 5, -0.5, 4.5);
	histos_["n_cosmic_muons"] = theFileService->make<TH1F> ("n_cosmic_muons", "n_cosmic_muons", 5, -0.5, 4.5);
	histos_["n_hodigis"] = theFileService->make<TH1F> ("n_hodigis", "n_hodigis", 50000, 0., 50000.);
	histos_["n_crossedHOIds"] = theFileService->make<TH1F> ("n_crossedHOIds", "n_crossedHOIds", 11, -0.5, 10.5);
	
	histos_["tenTS_fC_onechannel"] = theFileService->make<TH1F> ("tenTS_fC_onechannel", "tenTS_fC_onechannel", 100, 0., 100.);
	histos_["tempindex_onechannel"] = theFileService->make<TH1I> ("tempindex_onechannel", "tempindex_onechannel", 10, 0, 10);
	histos_["peak_fC_onechannel"] = theFileService->make<TH1F> ("peak_fC_onechannel", "peak_fC_onechannel", 50, 0., 50.);
	histos_["fourTS_fC_onechannel"] = theFileService->make<TH1F> ("fourTS_fC_onechannel", "fourTS_fC_onechannel", 100, 0., 100.);

	histos_["fourTSfC_muonchannel_exgrid"] = theFileService->make<TH1F> ("fourTSfC_muonchannel_exgrid", "fourTSfC_muonchannel_exgrid", 1000, 0., 1000.);

	histos_["tenTS_fC_muonchannel"] = theFileService->make<TH1F> ("tenTS_fC_muonchannel", "tenTS_fC_muonchannel", 500, 0., 500.);
	histos_["tempindex_muonchannel"] = theFileService->make<TH1I> ("tempindex_muonchannel", "tempindex_muonchannel", 10, 0, 10);
	histos_["peak_fC_muonchannel"] = theFileService->make<TH1F> ("peak_fC_muonchannel", "peak_fC_muonchannel", 50, 0., 50.);

	histos_["fourTS_fC_muonchannel_etaminus"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaminus", "fourTS_fC_muonchannel_etaminus", 1000, 0., 1000.);
	histos_["fourTS_fC_muonchannel_etaplus"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaplus", "fourTS_fC_muonchannel_etaplus", 1000, 0., 1000.);
	histos_["fourTS_fC_onechannel_muonchannel"] = theFileService->make<TH1F> ("fourTS_fC_onechannel_muonchannel", "fourTS_fC_onechannel_muonchannel", 1000, 0., 1000.);
	
	histos_["fourTS_fC_muonchannel_60to1000"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_60to1000", "fourTS_fC_muonchannel_60to1000", 1000, 0., 1000.);

	histos_["fourTS_fC_muonchannel_etaminus_grid"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaminus_grid", "fourTS_fC_muonchannel_etaminus_grid", 1000, 0., 1000.);
	histos_["fourTS_fC_muonchannel_etaplus_grid"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaplus_grid", "fourTS_fC_muonchannel_etaplus_grid", 1000, 0., 1000.);
	histos_["fourTS_fC_onechannel_muonchannel_grid"] = theFileService->make<TH1F> ("fourTS_fC_onechannel_muonchannel_grid", "fourTS_fC_onechannel_muonchannel_grid", 1000, 0., 1000.);

	histos_["index_vs_nomfc_onechannel"] = theFileService->make<TH1F> ("index_vs_nomfc_onechannel", "index_vs_nomfc_onechannel", 10, -0.5, 9.5);
	histos_["index_vs_nomfc_muonchannel"] = theFileService->make<TH1F> ("index_vs_nomfc_muonchannel", "index_vs_nomfc_muonchannel", 10, -0.5, 9.5);

	histos_["index_vs_fourfc"] = theFileService->make<TH2F> ("index_vs_fourfc", "index_vs_fourfc", 10, 0., 10., 500, 0., 500.);

	graphs_["muon_track_yz"] = theFileService->make<TGraph> ();
	graphs_["muon_track_xy"] = theFileService->make<TGraph> ();

	histos_["angle"] = theFileService->make<TH1F> ("angle", "angle", 100, -1.570795, 1.570795);

	histos_["theID"] = theFileService->make<TH1F> ("theID", "theID", 10000, 1174470000, 1174480000);

	histos_["cosmu_HOeta"] = theFileService->make<TH1F> ("cosmu_innereta", "cosmu_innereta", 120, -1.5, 1.5);
	histos_["cosmu_HOphi"] = theFileService->make<TH1F> ("cosmu_outereta", "cosmu_outereta", 120, -1.5, 1.5);
}

HODIGIAnalyzer::~HODIGIAnalyzer() {

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

//void HODIGIAnalyzer::getCalibHODigis(const edm::EventSetup& iSetup, HcalDetId aId, edm::SortedCollection<HODataFrame>::const_iterator aIt, CaloSamples aTool){
//	// get calibration info
//	edm::ESHandle < HcalDbService > HCalconditions;
//	iSetup.get<HcalDbRecord> ().get(HCalconditions);
//
//	if (!HCalconditions.isValid()) return;
//
//	const HcalQIEShape *shape = HCalconditions->getHcalShape();
//
//	const HcalCalibrations& calibrations = HCalconditions->getHcalCalibrations(aId);
//	const HcalQIECoder *channelCoder = HCalconditions->getHcalCoder(aId);
//	HcalCoderDb coder(*channelCoder, *shape);
//	coder.adc2fC(*aIt, aTool);
//
//	for (int ii = 0; ii < aTool.size(); ++ii) {
//		// default ped is 4.5
//		int capid = (*aIt)[ii].capid();
//		fDigiSum1 += (aTool[ii] - calibrations.pedestal(capid));
//		fDigiSum_ped1 += aTool[ii];
//	}
//}

// ------------ method called for each event  ------------
void HODIGIAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;
	uint32_t tempID1 =        1174470824;
	//uint32_t tempID2 =        1174470952;
	//uint32_t tempID3 =        1174471080;
	//uint32_t tempID4 =        1174471208;
	//uint32_t tempID5 =        1174471336;
	//uint32_t tempID = 1174470814; //for run 1 analysis ieta -1, iphi 30
	//uint32_t temp_onechannel_muonID = 1174479121; //for GRIN analysis ieta 2, iphi 17
			// 1174479116; //for GRIN analysis ieta 2, iphi 12

	Handle <HODataFrameCollection> hodigis;
	iEvent.getByLabel("hcalDigis", hodigis);

	Handle < std::vector < reco::Muon >> muons;
	iEvent.getByLabel("MuonSelector", muons);
	
	Handle < std::vector < reco::Track >> cosmic_muons;
	iEvent.getByLabel("cosmicMuons", cosmic_muons);

	if (!MuonHOAcceptance::Inited()) MuonHOAcceptance::initIds(iSetup);

	histos_["n_muons"]->Fill(muons->size());
	histos_["n_cosmic_muons"]->Fill(cosmic_muons->size());
	histos_["n_hodigis"]->Fill(hodigis->size());

	HODataFrameCollection::const_iterator it1;

	for (it1 = hodigis->begin(); it1 != hodigis->end(); ++it1) {// loop over all HODataFrames
		HcalDetId digiID(it1->id());

		if (digiID.rawId() == tempID1){ //go on if digiID = desired ID
			double tenTSfc_onechannel = 0;
			double tempfc_onechannel = 0.;
			int tempindex_onechannel = 0;
			double fourTSfC_onechannel = 0.;
			
			//plot the nominal fC
			for (int ii = 0; ii < (*it1).size(); ++ii) {
				histos_["index_vs_nomfc_onechannel"]->SetBinContent(ii+1,(*it1)[ii].nominal_fC());

				tenTSfc_onechannel += (*it1)[ii].nominal_fC();
				if (tempfc_onechannel < (*it1)[ii].nominal_fC()){
					tempindex_onechannel = ii;
					tempfc_onechannel = (*it1)[ii].nominal_fC();
				}
			}
			
			if(tempindex_onechannel > 0 && tempindex_onechannel <= 7){
				fourTSfC_onechannel = (*it1)[tempindex_onechannel-1].nominal_fC()+(*it1)[tempindex_onechannel].nominal_fC()+
						(*it1)[tempindex_onechannel+1].nominal_fC()+(*it1)[tempindex_onechannel+2].nominal_fC();
			} else fourTSfC_onechannel = (*it1)[3].nominal_fC()+(*it1)[4].nominal_fC()+(*it1)[5].nominal_fC()+(*it1)[6].nominal_fC();

			histos_["tempindex_onechannel"]->Fill(tempindex_onechannel);
			histos_["peak_fC_onechannel"]->Fill(tempfc_onechannel);
			histos_["tenTS_fC_onechannel"]->Fill(tenTSfc_onechannel);
			histos_["fourTS_fC_onechannel"]->Fill(fourTSfC_onechannel);
		}
	}
	//std::cout << "xxxxxxxxxxxxxxxxxxx" << std::endl;
	
	//if (muons->size() > 0) { // do this if at least one muon there
	//	for (unsigned int i = 0; i < muons->size(); ++i) { //loop over all muons
	if (cosmic_muons->size() >= 1) { // do this if exactly one cosmic muon there
		for (unsigned int i = 0; i < cosmic_muons->size(); ++i) { //loop over all muons

	//		reco::Muon const * muon_iter = &(muons->at(i));

	//		//analyze only tight muons
	//		if (!muon::isTightMuon(*muon_iter, getvertex(iEvent))) return;
			
	//		if (muon_iter->track().isNull()) continue;
	//		reco::Track const * track = muon_iter->track().get();

	//		TrackDetMatchInfo * muMatch = new TrackDetMatchInfo(assoc.associate(iEvent, iSetup, *track, assocParams));
			TrackDetMatchInfo * muMatch = new TrackDetMatchInfo(assoc.associate(iEvent, iSetup, cosmic_muons->at(i), assocParams));

			/*
			double hophi = muMatch->trkGlobPosAtHO.Phi();
			double hoeta = muMatch->trkGlobPosAtHO.Eta();

			// muons in geom.accept of HO
			//if (!theMuonHOAcceptance->inGeomAccept(hoeta2, hophi2, 0.04, 0.017)) inAccept = 0;
			if (!theMuonHOAcceptance->inGeomAccept(hoeta, hophi, 0., 0.)) continue;
			// muons in non dead cells of HO
			//if (!theMuonHOAcceptance->inNotDeadGeom(hoeta2, hophi2, 0.04, 0.017)) inNotDead = 0;
			if (!theMuonHOAcceptance->inNotDeadGeom(hoeta, hophi, 0., 0.)) continue;
			// muons in SiPM tiles of HO
			//if (!theMuonHOAcceptance->inSiPMGeom(hoeta2, hophi2, 0.04, 0.017)) inSiPM = 0;
			 */
			histos_["n_crossedHOIds"]->Fill(muMatch->crossedHOIds.size());

			//int cosmu_HOieta = -999;
			//int cosmu_HOiphi = -999;

			// if muon didn't cross any HO tile go on with the next muon
			if (muMatch->crossedHOIds.size() == 0) continue;

			double hophi = muMatch->trkGlobPosAtHO.Phi();
			double hoeta = muMatch->trkGlobPosAtHO.Eta();

			for (int ieta = 0; ieta<etaBounds; ++ieta) {
				if ( (hoeta > etaMin[ieta]) && (hoeta < etaMax[ieta]) ) {
					int const phiSectors = ((ieta == 1) ? phiSectorsR0 : phiSectorsR12);
					for (int iphi = 0; iphi<phiSectors; ++iphi) {
						double const * mins = ((ieta == 1) ? phiMinR0 : phiMinR12);
						double const * maxes = ((ieta == 1) ? phiMaxR0 : phiMaxR12);
						while (hophi < mins[0]) hophi += twopi;
						while (hophi > mins[0]+twopi) hophi -= twopi;
						if ((hophi > mins[iphi]) && (hophi < maxes[iphi])) {
							int tile_ieta = -99.;

							double fourTSfC_muonchannel = 0.;
							double temp_fourTSfC_muonchannel = 0.;
							int tempindex_muonchannel;
							int index_muonchannel = -99.;

							double fourTSfC_muonchannel_grid = 0.;
							double fourTSfC_muonchannel_grid1 = 0.;
							double temp_fourTSfC_muonchannel_grid = 0.;

							double temp_fourTSfC_muonchannel_exgrid = 0.;

							HODataFrameCollection::const_iterator it2;

							for (std::vector<DetId>::const_iterator aid = muMatch->crossedHOIds.begin(); aid != muMatch->crossedHOIds.end(); ++aid) {// loop over all HOIds crossed by the muon
								HcalDetId mId(aid->rawId());

								for (it2 = hodigis->begin(); it2 != hodigis->end(); ++it2) {
									HcalDetId cell(it2->id());

									if (cell.rawId() == mId.rawId()){ // do something when the crossedId is the same as the Id of the current digi

										histos_["theID"]->Fill(mId.rawId());

										/*
										double angle_for_corr = 0.;
										double theta = 30.;
										const math::XYZVector& cosmic_mom = cosmic_muons->at(i).momentum();
										//const math::XYZVector& cosmic_mom = muon_iter->momentum();
										math::XYZVector normal_vec;

										if(mId.iphi() >= 1 && mId.iphi() <= 6){
											normal_vec.SetXYZ(sin(10*theta),cos(10*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 7 && mId.iphi() <= 12){
											normal_vec.SetXYZ(sin(11*theta),cos(11*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 13 && mId.iphi() <= 18){
											normal_vec.SetXYZ(0.,1.,0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 19 && mId.iphi() <= 24){
											normal_vec.SetXYZ(sin(theta),cos(theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 25 && mId.iphi() <= 30){
											normal_vec.SetXYZ(sin(2*theta),cos(2*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 31 && mId.iphi() <= 36){
											normal_vec.SetXYZ(sin(3*theta),cos(3*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 37 && mId.iphi() <= 42){
											normal_vec.SetXYZ(sin(4*theta),cos(4*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 43 && mId.iphi() <= 48){
											normal_vec.SetXYZ(sin(5*theta),cos(5*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 49 && mId.iphi() <= 54){
											normal_vec.SetXYZ(sin(6*theta),cos(6*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 55 && mId.iphi() <= 60){
											normal_vec.SetXYZ(sin(7*theta),cos(7*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 61 && mId.iphi() <= 66){
											normal_vec.SetXYZ(sin(8*theta),cos(8*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										} else if(mId.iphi() >= 67 && mId.iphi() <= 72){
											normal_vec.SetXYZ(sin(9*theta),cos(9*theta),0.);
											angle_for_corr = angle(normal_vec,cosmic_mom);
											if(angle_for_corr > 1.570795) angle_for_corr -= 1.570795;
										}

										histos_["angle"]->Fill(angle_for_corr);
										 */
										double tempfc_muonchannel = 0.;
										tempindex_muonchannel = 0;

										double tempfc_muonchannel_grid = 0.;
										int tempindex_muonchannel_grid = 0;

										double tempfc_muonchannel_exgrid = 0.;
										int tempindex_muonchannel_exgrid = 0;

										int gridSize = 1;
										int n_gridtiles = 2;

										//plot the nominal fC
										for (int ii = 0; ii < (*it2).size(); ++ii) {
											histos_["index_vs_nomfc_muonchannel"]->SetBinContent(ii+1,(*it2)[ii].nominal_fC());
											if (tempfc_muonchannel < (*it2)[ii].nominal_fC()){
												tempindex_muonchannel = ii;
												tempfc_muonchannel = (*it2)[ii].nominal_fC();
											}
										}

										// integrate over 4 timeslices
										if(tempindex_muonchannel > 0 && tempindex_muonchannel <= 7){
											temp_fourTSfC_muonchannel = (*it2)[tempindex_muonchannel-1].nominal_fC()+(*it2)[tempindex_muonchannel].nominal_fC()+
												(*it2)[tempindex_muonchannel+1].nominal_fC()+(*it2)[tempindex_muonchannel+2].nominal_fC();
										} else temp_fourTSfC_muonchannel = (*it2)[3].nominal_fC()+(*it2)[4].nominal_fC()+(*it2)[5].nominal_fC()+(*it2)[6].nominal_fC();

										if (temp_fourTSfC_muonchannel > fourTSfC_muonchannel) {
											tile_ieta = mId.ieta();
											fourTSfC_muonchannel = temp_fourTSfC_muonchannel;
											index_muonchannel = tempindex_muonchannel;
										}

										//NxN (N=gridSize) neighborhood of the tile where the muon is expected to go through
										HODataFrameCollection::const_iterator it3;

										for (it3 = hodigis->begin(); it3 != hodigis->end(); ++it3) {
											HcalDetId neighborId(it3->id());
											int dEta = abs( (mId.ieta()<0?mId.ieta()+1:mId.ieta() )
															-(neighborId.ieta()<0?neighborId.ieta()+1:neighborId.ieta() ) ) ;
											int dPhi = abs( mId.iphi()-neighborId.iphi() );
											if (abs(72-dPhi) < dPhi) dPhi = 72-dPhi;

											if(dEta <= gridSize && dPhi <= gridSize) {//within the grid
												n_gridtiles++;
												for (int iii = 0; iii < (*it3).size(); ++iii) {
													if (tempfc_muonchannel_grid < (*it3)[iii].nominal_fC()){
														tempindex_muonchannel_grid = iii;
														tempfc_muonchannel_grid = (*it3)[iii].nominal_fC();
													}
												}

												// integrate over 4 timeslices
												if(tempindex_muonchannel_grid > 0 && tempindex_muonchannel_grid <= 7){
													temp_fourTSfC_muonchannel_grid = (*it3)[tempindex_muonchannel_grid-1].nominal_fC()+(*it3)[tempindex_muonchannel_grid].nominal_fC()+
														(*it3)[tempindex_muonchannel_grid+1].nominal_fC()+(*it3)[tempindex_muonchannel_grid+2].nominal_fC();
												} else temp_fourTSfC_muonchannel_grid = (*it3)[3].nominal_fC()+(*it3)[4].nominal_fC()+(*it3)[5].nominal_fC()+(*it3)[6].nominal_fC();

												if (temp_fourTSfC_muonchannel_grid > fourTSfC_muonchannel_grid) fourTSfC_muonchannel_grid = temp_fourTSfC_muonchannel_grid;

											} else {//outside the grid
												for (int iiii = 0; iiii < (*it3).size(); ++iiii) {
													if (tempfc_muonchannel_exgrid < (*it3)[iiii].nominal_fC()){
														tempindex_muonchannel_exgrid = iiii;
														tempfc_muonchannel_exgrid = (*it3)[iiii].nominal_fC();
													}
												}

												// integrate over 4 timeslices
												if(tempindex_muonchannel_exgrid > 0 && tempindex_muonchannel_exgrid <= 7){
													temp_fourTSfC_muonchannel_exgrid = (*it3)[tempindex_muonchannel_exgrid-1].nominal_fC()+(*it3)[tempindex_muonchannel_exgrid].nominal_fC()+
														(*it3)[tempindex_muonchannel_exgrid+1].nominal_fC()+(*it3)[tempindex_muonchannel_exgrid+2].nominal_fC();
												} else temp_fourTSfC_muonchannel_exgrid = (*it3)[3].nominal_fC()+(*it3)[4].nominal_fC()+(*it3)[5].nominal_fC()+(*it3)[6].nominal_fC();

												histos_["fourTSfC_muonchannel_exgrid"]->Fill(temp_fourTSfC_muonchannel_exgrid);
											}
										}

										if (fourTSfC_muonchannel_grid > fourTSfC_muonchannel_grid1) {
											tile_ieta = mId.ieta();
											fourTSfC_muonchannel_grid1 = fourTSfC_muonchannel_grid;
										}
									}
								}
							}

							if (tile_ieta < -4){
								histos_["fourTS_fC_muonchannel_etaminus"]->Fill(fourTSfC_muonchannel);
							} else {
								histos_["index_vs_fourfc"]->Fill(index_muonchannel,fourTSfC_muonchannel);
								histos_["fourTS_fC_muonchannel_etaplus"]->Fill(fourTSfC_muonchannel);
							}

							if (tile_ieta == 10 || tile_ieta == 11){
								histos_["fourTS_fC_onechannel_muonchannel"]->Fill(fourTSfC_muonchannel);
							}


							// plot negative wheels seperately
							if (tile_ieta < -4){
								histos_["fourTS_fC_muonchannel_etaminus_grid"]->Fill(fourTSfC_muonchannel_grid1);
							} else {
								histos_["fourTS_fC_muonchannel_etaplus_grid"]->Fill(fourTSfC_muonchannel_grid1);

								if(fourTSfC_muonchannel_grid1 >= 60.){
									histos_["fourTS_fC_muonchannel_60to1000"]->Fill(fourTSfC_muonchannel-36.);
								}
							}

							if (tile_ieta == 10 || tile_ieta == 11){
								histos_["fourTS_fC_onechannel_muonchannel_grid"]->Fill(fourTSfC_muonchannel_grid1);
							}

							//histos_["tempindex_muonchannel"]->Fill(tempindex_muonchannel);
							//histos_["peak_fC_muonchannel"]->Fill(tempfc_muonchannel);

							delete muMatch;
						}
					}
				}
			}
		}
	}

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}
/*
reco::Vertex HODIGIAnalyzer::getvertex(const edm::Event& iEvent){
	//-------------------------------Vertex information------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	using namespace edm;

	reco::Vertex::Point posVtx;
	reco::Vertex::Error errVtx;

	edm::Handle<reco::VertexCollection> recVtxs;
	iEvent.getByLabel(primvertexLabel_, recVtxs);

	unsigned int theIndexOfThePrimaryVertex = 999.;

	for (unsigned int ind = 0; ind < recVtxs->size(); ++ind) {
		if ((*recVtxs)[ind].isValid() && !((*recVtxs)[ind].isFake())) {
			theIndexOfThePrimaryVertex = ind;
			break;
		}
	}

	if (theIndexOfThePrimaryVertex < 100) {
		posVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).position();
		errVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).error();
	} else {
		LogInfo("RecoMuonValidator")
			<< "reco::PrimaryVertex not found, use BeamSpot position instead\n";
		edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
		iEvent.getByLabel(beamspotLabel_, recoBeamSpotHandle);
		reco::BeamSpot bs = *recoBeamSpotHandle;
		posVtx = bs.position();
		errVtx(0, 0) = bs.BeamWidthX();
		errVtx(1, 1) = bs.BeamWidthY();
		errVtx(2, 2) = bs.sigmaZ();
	}

	const reco::Vertex thePrimaryVertex(posVtx, errVtx);
	const GlobalPoint vertex(thePrimaryVertex.x(), thePrimaryVertex.y(), thePrimaryVertex.z());

	return thePrimaryVertex;
}
*/

// ------------ method called once each job just before starting event loop  ------------
void HODIGIAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void HODIGIAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void HODIGIAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void HODIGIAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void HODIGIAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HODIGIAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HODIGIAnalyzer::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE( HODIGIAnalyzer);
