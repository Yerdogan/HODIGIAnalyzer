void eff_pur_analyzer(){
	TFile *_file0 = TFile::Open("SingleMu_Run2012A_22Jan2013_v1_RECO_HODIGI.root");

	TH1F *hist1 = (TH1F*)_file0->Get("demo/fourTS_fC_muonchannel_etaplus_grid");
	TH1F *hist2 = (TH1F*)_file0->Get("demo/fourTSfC_muonchannel_exgrid");
	
	TGraph *eff = new TGraph();
	TGraph *pur = new TGraph();

	double hist1_wholeintegral = hist1->Integral(0,1000);
	double hist2_wholeintegral = hist2->Integral(0,1000);
	
	for (int i = 0; i <= 1000; i++){
		double hist1_currintegral = hist1->Integral(i,1000);
		double hist2_currintegral = hist2->Integral(0,i);
		eff->SetPoint(i,i,hist1_currintegral/hist1_wholeintegral);
		pur->SetPoint(i,i,hist2_currintegral/hist2_wholeintegral);
	}
	
	eff->Draw("AC*");
	pur->Draw("same");
}