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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/MuonIdentification/interface/MuonHOAcceptance.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

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
	reco::Vertex getvertex(const edm::Event&);
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

edm::InputTag beamspotLabel_;
edm::InputTag primvertexLabel_;

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
	//primvertexLabel_ = iConfig.getParameter<edm::InputTag> ("primaryVertex");

	//histos_["nevent_alive_after_cut"] = theFileService->make<TH1F> ( "nevent_alive_after_cut", "nevent_alive_after_cut", 7, -0.5, 6.5);

	histos_["n_muons"] = theFileService->make<TH1F> ("n_muons", "n_muons", 5, -0.5, 4.5);
	histos_["n_cosmic_muons"] = theFileService->make<TH1F> ("n_cosmic_muons", "n_cosmic_muons", 5, -0.5, 4.5);
	histos_["n_hodigis"] = theFileService->make<TH1F> ("n_hodigis", "n_hodigis", 50000, 0., 50000.);
	histos_["n_crossedHOIds"] = theFileService->make<TH1F> ("n_crossedHOIds", "n_crossedHOIds", 11, -0.5, 10.5);
	
	histos_["hodigi_fC_onechannel"] = theFileService->make<TH1F> ("hodigi_fC_onechannel", "hodigi_fC_onechannel", 100, 0., 100.);
	histos_["tempindex_onechannel"] = theFileService->make<TH1F> ("tempindex_onechannel", "tempindex_onechannel", 11, -0.5, 9.5);
	histos_["peak_fC_onechannel"] = theFileService->make<TH1F> ("peak_fC_onechannel", "peak_fC_onechannel", 50, 0., 50.);
	histos_["alternative_fC_onechannel"] = theFileService->make<TH1F> ("alternative_fC_onechannel", "alternative_fC_onechannel", 100, 0., 100.);

	histos_["hodigi_fC_muonchannel"] = theFileService->make<TH1F> ("hodigi_fC_muonchannel", "hodigi_fC_muonchannel", 500, 0., 500.);
	histos_["tempindex_muonchannel"] = theFileService->make<TH1F> ("tempindex_muonchannel", "tempindex_muonchannel", 11, -0.5, 9.5);
	histos_["peak_fC_muonchannel"] = theFileService->make<TH1F> ("peak_fC_muonchannel", "peak_fC_muonchannel", 50, 0., 50.);
	histos_["alternative_fC_muonchannel"] = theFileService->make<TH1F> ("alternative_fC_muonchannel", "alternative_fC_muonchannel", 500, 0., 500.);
	
	histos_["index_vs_nomfc_onechannel"] = theFileService->make<TH2F> ("index_vs_nomfc_onechannel", "index_vs_nomfc_onechannel", 11, -0.5, 10.5, 100, 0., 100.);
	histos_["index_vs_nomfc_muonchannel"] = theFileService->make<TH2F> ("index_vs_nomfc_muonchannel", "index_vs_nomfc_muonchannel", 11, -0.5, 10.5,100, 00., 100.);

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
			double nomfc_onechannel = 0;
			double tempfc_onechannel = 0.;
			int tempindex_onechannel = 0;
			double alternativefC_onechannel = 0.;
			
			//plot the nominal fC
			for (int ii = 0; ii < (*it1).size(); ++ii) {
				histos_["index_vs_nomfc_onechannel"]->Fill(ii,(*it1)[ii].nominal_fC());
				nomfc_onechannel += (*it1)[ii].nominal_fC();
				if (tempfc_onechannel < (*it1)[ii].nominal_fC()){
					tempindex_onechannel = ii;
					tempfc_onechannel = (*it1)[ii].nominal_fC();
				}
			}
			
			if(tempindex_onechannel > 0 && tempindex_onechannel <= 7){
				alternativefC_onechannel = (*it1)[tempindex_onechannel-1].nominal_fC()+(*it1)[tempindex_onechannel].nominal_fC()+
						(*it1)[tempindex_onechannel+1].nominal_fC()+(*it1)[tempindex_onechannel+2].nominal_fC();
			} else alternativefC_onechannel = (*it1)[3].nominal_fC()+(*it1)[4].nominal_fC()+(*it1)[5].nominal_fC()+(*it1)[6].nominal_fC();
				
			histos_["tempindex_onechannel"]->Fill(tempindex_onechannel);
			histos_["peak_fC_onechannel"]->Fill(tempfc_onechannel);
			histos_["hodigi_fC_onechannel"]->Fill(nomfc_onechannel);
			histos_["alternative_fC_onechannel"]->Fill(alternativefC_onechannel);
		}
	}
	//std::cout << "xxxxxxxxxxxxxxxxxxx" << std::endl;
	
	//if (muons->size() > 0) { // do this if at least one muon there
	//	for (unsigned int i = 0; i < muons->size(); ++i) { //loop over all muons
	if (cosmic_muons->size() > 0) { // do this if at least one muon there
		for (unsigned int i = 0; i < cosmic_muons->size(); ++i) { //loop over all muons		
			//reco::Muon const * muon_iter = &(muons->at(i));
			
			//analyze only tight muons
			//if (!muon::isTightMuon(*muon_iter, getvertex(iEvent))) return;

			//if (muon_iter->track().isNull()) continue;
			//reco::Track const * track = muon_iter->track().get();

			//TrackDetMatchInfo * muMatch = new TrackDetMatchInfo(assoc.associate(iEvent, iSetup, *track, assocParams));
			TrackDetMatchInfo * muMatch = new TrackDetMatchInfo(assoc.associate(iEvent, iSetup, cosmic_muons->at(i), assocParams));

			histos_["n_crossedHOIds"]->Fill(muMatch->crossedHOIds.size());

			// if muon didn't cross any HO tile go on with the next muon
			if (muMatch->crossedHOIds.size() == 0) continue;
			
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

			*/

			HODataFrameCollection::const_iterator it2;

			for (it2 = hodigis->begin(); it2 != hodigis->end(); ++it2) { //another loop over all HODataFrames
				HcalDetId cell(it2->id());

				if (cell.subdet() == 3) {// 3 = HO
					for (std::vector<DetId>::const_iterator aid = muMatch->crossedHOIds.begin(); aid != muMatch->crossedHOIds.end(); ++aid) {
						// loop over all HOIds crossed by the muon
						HcalDetId mId(aid->rawId());

						if (mId == cell) { // do something when the crossedId is the same as the Id of the current digi
							double nomfc_muonchannel = 0;
							double tempfc_muonchannel = 0.;
							int tempindex_muonchannel = 0;
							double alternativefC_muonchannel = 0.;
							
							//plot the nominal fC
							for (int ii = 0; ii < (*it2).size(); ++ii) {
								histos_["index_vs_nomfc_muonchannel"]->Fill(ii,(*it2)[ii].nominal_fC());
								nomfc_muonchannel += (*it2)[ii].nominal_fC();
								if (tempfc_muonchannel < (*it2)[ii].nominal_fC()){
									tempindex_muonchannel = ii;
									tempfc_muonchannel = (*it2)[ii].nominal_fC();
								}
							}
							
							if(tempindex_muonchannel > 0 && tempindex_muonchannel <= 7){
								alternativefC_muonchannel = (*it2)[tempindex_muonchannel-1].nominal_fC()+(*it2)[tempindex_muonchannel].nominal_fC()+
										(*it2)[tempindex_muonchannel+1].nominal_fC()+(*it2)[tempindex_muonchannel+2].nominal_fC();
							} else alternativefC_muonchannel = (*it2)[3].nominal_fC()+(*it2)[4].nominal_fC()+(*it2)[5].nominal_fC()+(*it2)[6].nominal_fC();
								
							histos_["tempindex_muonchannel"]->Fill(tempindex_muonchannel);
							histos_["peak_fC_muonchannel"]->Fill(tempfc_muonchannel);
							histos_["hodigi_fC_muonchannel"]->Fill(nomfc_muonchannel);
							histos_["alternative_fC_muonchannel"]->Fill(alternativefC_muonchannel);
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
