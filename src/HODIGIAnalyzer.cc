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
	const CaloGeometry *caloGeometry;

};

//
// constants, enums and typedefs
//

typedef edm::SortedCollection < HODataFrame > HODataFrameCollection;

//edm::InputTag beamspotLabel_;
//edm::InputTag primvertexLabel_;

//
// static data member definitions
//

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

	histos_["tenTS_fC_muonchannel"] = theFileService->make<TH1F> ("tenTS_fC_muonchannel", "tenTS_fC_muonchannel", 500, 0., 500.);
	histos_["tempindex_muonchannel"] = theFileService->make<TH1I> ("tempindex_muonchannel", "tempindex_muonchannel", 10, 0, 10);
	histos_["peak_fC_muonchannel"] = theFileService->make<TH1F> ("peak_fC_muonchannel", "peak_fC_muonchannel", 50, 0., 50.);
	histos_["fourTS_fC_muonchannel_etaminus"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaminus", "fourTS_fC_muonchannel_etaminus", 100, 0., 500.);
	histos_["fourTS_fC_muonchannel_etaplus"] = theFileService->make<TH1F> ("fourTS_fC_muonchannel_etaplus", "fourTS_fC_muonchannel_etaplus", 100, 0., 500.);
	
	histos_["index_vs_nomfc_onechannel"] = theFileService->make<TH1F> ("index_vs_nomfc_onechannel", "index_vs_nomfc_onechannel", 10, 0., 10.);
	histos_["index_vs_nomfc_muonchannel"] = theFileService->make<TH1F> ("index_vs_nomfc_muonchannel", "index_vs_nomfc_muonchannel", 10, 0., 10.);

	graphs_["muon_track_yz"] = theFileService->make<TGraph> ();
	graphs_["muon_track_xy"] = theFileService->make<TGraph> ();

	histos_["angle"] = theFileService->make<TH1F> ("angle", "angle", 100, -1.570795, 1.570795);
	histos_["fourTS_fC_onechannel_muonchannel"] = theFileService->make<TH1F> ("fourTS_fC_onechannel_muonchannel", "fourTS_fC_onechannel_muonchannel", 100, 0., 500.);

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

	ESHandle < CaloGeometry > caloGeom;
	iSetup.get<CaloGeometryRecord> ().get(caloGeom);

	caloGeometry = caloGeom.product();

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
				histos_["index_vs_nomfc_onechannel"]->SetBinContent(ii,(*it1)[ii].nominal_fC());
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
	if (cosmic_muons->size() > 0) { // do this if at least one muon there
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

			int cosmu_HOieta = -999;
			int cosmu_HOiphi = -999;

			//alternative matching
			if (cosmic_muons->at(i).innerOk() && cosmic_muons->at(i).outerOk()){

				math::XYZPoint cosmu_innerpoint = cosmic_muons->at(i).innerPosition();
				math::XYZPoint cosmu_outerpoint = cosmic_muons->at(i).outerPosition();

				// inner = entrance, outer = exit
				double cosmu_innerx = cosmu_innerpoint.x();
				double cosmu_innery = cosmu_innerpoint.y();
				double cosmu_innerz = cosmu_innerpoint.z();

				double cosmu_outerx = cosmu_outerpoint.x();
				double cosmu_outery = cosmu_outerpoint.y();
				double cosmu_outerz = cosmu_outerpoint.z();

				graphs_["muon_track_yz"]->SetPoint(0,cosmu_innerz,cosmu_innery);
				graphs_["muon_track_yz"]->SetPoint(1,cosmu_outerz,cosmu_outery);

				graphs_["muon_track_yz"]->Fit("pol1","q");

				graphs_["muon_track_xy"]->SetPoint(0,cosmu_innerx,cosmu_innery);
				graphs_["muon_track_xy"]->SetPoint(1,cosmu_outerx,cosmu_outery);

				graphs_["muon_track_xy"]->Fit("pol1","q");

				TF1 *fit_fuyz = new TF1();
				fit_fuyz = graphs_["muon_track_yz"]->GetFunction("pol1");

				TF1 *fit_fuxy = new TF1();
				fit_fuxy = graphs_["muon_track_xy"]->GetFunction("pol1");

				double HO_y = -99.;

				if(cosmic_muons->at(i).momentum().y() < 0.){
					HO_y = 407.;
				} else HO_y = -407.;

				double HO_z = fit_fuyz->GetX(HO_y,-800.,800.,1.E-10,100,false);
				double HO_x = fit_fuxy->GetX(HO_y,-800.,800.,1.E-10,100,false);

				GlobalPoint cosmu_HOpoint_yz(0.,HO_y,HO_z);
				GlobalPoint cosmu_HOpoint_xy(HO_x,HO_y,0.);

				//eta plane
				cosmu_HOieta = int(cosmu_HOpoint_yz.eta()/0.087) + ((cosmu_HOpoint_yz.eta()>0) ? 1 : -1);

				//phi plane
				cosmu_HOiphi = int(cosmu_HOpoint_xy.phi()/0.087);
			}


			// if muon didn't cross any HO tile go on with the next muon
			if (muMatch->crossedHOIds.size() == 0) continue;
			
			// ----------------------------------------------------------------
/*
			//consider only muons without hits in CSC
			int total = 0;
			bool foundcsc = 0;
			int nAll = muon_iter->numberOfChambers();

			for (int iC = 0; iC < nAll && foundcsc == 0; ++iC) {
				if (muon_iter->matches()[iC].detector() == MuonSubdetId::CSC) {
					foundcsc = 1;
					break;
				}
				if (muon_iter->matches()[iC].detector() == MuonSubdetId::DT)
					total++;
			}

			if (foundcsc || total == 0) continue;

			// ------------------------------------------------------------------
*/
			HODataFrameCollection::const_iterator it2;

			//for (std::vector<DetId>::const_iterator aid = muMatch->crossedHOIds.begin(); aid != muMatch->crossedHOIds.end(); ++aid) {// loop over all HOIds crossed by the muon
				//HcalDetId mId(aid->rawId());

				for (it2 = hodigis->begin(); it2 != hodigis->end(); ++it2) {
					HcalDetId cell(it2->id());

					//if (cell.rawId() == mId.rawId()){
					//if (cell.ieta() == mId.ieta() && cell.iphi() == mId.iphi()){
					if (cell.ieta() == cosmu_HOieta && cell.iphi() == cosmu_HOiphi){
						// do something when the crossedId is the same as the Id of the current digi
						/*
						histos_["theID"]->Fill(mId.rawId());

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
						double tenTSfc_muonchannel = 0;
						double tempfc_muonchannel = 0.;
						int tempindex_muonchannel = 0;
						double fourTSfC_muonchannel = 0.;

						//plot the nominal fC
						for (int ii = 0; ii < (*it2).size(); ++ii) {
							histos_["index_vs_nomfc_muonchannel"]->SetBinContent(ii,(*it2)[ii].nominal_fC());
							tenTSfc_muonchannel += (*it2)[ii].nominal_fC();
							if (tempfc_muonchannel < (*it2)[ii].nominal_fC()){
								tempindex_muonchannel = ii;
								tempfc_muonchannel = (*it2)[ii].nominal_fC();
							}
						}

						// integrate over 4 timeslices
						if(tempindex_muonchannel > 0 && tempindex_muonchannel <= 7){
							fourTSfC_muonchannel = (*it2)[tempindex_muonchannel-1].nominal_fC()+(*it2)[tempindex_muonchannel].nominal_fC()+
								(*it2)[tempindex_muonchannel+1].nominal_fC()+(*it2)[tempindex_muonchannel+2].nominal_fC();
						} else fourTSfC_muonchannel = (*it2)[3].nominal_fC()+(*it2)[4].nominal_fC()+(*it2)[5].nominal_fC()+(*it2)[6].nominal_fC();

						// plot negative wheels seperately
						if (cosmu_HOieta < -4){
							histos_["fourTS_fC_muonchannel_etaminus"]->Fill(fourTSfC_muonchannel);
						} else histos_["fourTS_fC_muonchannel_etaplus"]->Fill(fourTSfC_muonchannel);

						if (cosmu_HOieta == 10 || cosmu_HOieta == 11){
							histos_["fourTS_fC_onechannel_muonchannel"]->Fill(fourTSfC_muonchannel);
						}
						/*else {
							if(fourTSfC_muonchannel > 60.){
								histos_["fourTS_fC_muonchannel_etaplus"]->Fill(fourTSfC_muonchannel*cos(angle_for_corr));
							}

							if (cosmu_HOieta == 4 || cosmu_HOieta == 10 || cosmu_HOieta == 11){
								histos_["fourTS_fC_onechannel_muonchannel"]->Fill(fourTSfC_muonchannel);
							}
						}*/

						histos_["tempindex_muonchannel"]->Fill(tempindex_muonchannel);
						histos_["peak_fC_muonchannel"]->Fill(tempfc_muonchannel);
						histos_["tenTS_fC_muonchannel"]->Fill(tenTSfc_muonchannel);
					}
				}
			//}
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
