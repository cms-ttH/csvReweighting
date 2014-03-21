//void overlayCSV () {
{
#include <cmath>
#include <algorithm>

  gStyle->SetOptStat(0);
  TString ofName = "csv_rwt_hf_v0.root";
  TFile * outputFile = new TFile(ofName, "RECREATE");

//   TString WP = "000"; //100x probe jet working point - loose (224), medium (679), tight (898)
  bool EfficiencySF = false; //Set to true if you are measuring the efficiency for different CSV working points
  TString pt_lower = "0";
  TString pt_upper = "10000";
  TString eta_lower = "0.0";
  TString eta_upper = "2.4";
  vector<TString> ptBins; 
  ptBins.push_back("30");
  ptBins.push_back("40");
  ptBins.push_back("60");
  ptBins.push_back("100");
  ptBins.push_back("160");
  ptBins.push_back("10000");
  unsigned int nPt = ptBins.size();
  vector<TString> etaBins;
  etaBins.push_back("0.0");
//   etaBins.push_back("0.8");
//   etaBins.push_back("1.6");
  etaBins.push_back("2.4");
  unsigned int nEta = etaBins.size();

  ///////
  for (unsigned int iPt = 0; iPt < nPt-1; iPt++){
    pt_lower = ptBins[iPt];
    pt_upper = ptBins[iPt+1];
    cout << "pt range is " << pt_lower << " to " << pt_upper << endl;
    for (unsigned int iEta = 0; iEta < nEta-1; iEta++){
      eta_lower = etaBins[iEta];
      eta_upper = etaBins[iEta+1];
      cout << "eta range is " << eta_lower << " to " << eta_upper << endl;

      ///input root files ---- change file path ---
  TFile * sourceFile1a = new TFile ("/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/TwoMuon.root");
  TFile * sourceFile1b = new TFile ("/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/TwoEle.root");
  TFile * sourceFile1c = new TFile ("/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/MuonEle.root");
  TFile * sourceFile2 = new TFile  ("/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/allMC.root");
  TTree * Tree1a = (TTree*) sourceFile1a->Get("summaryTree");
  TTree * Tree1b = (TTree*) sourceFile1b->Get("summaryTree");
  TTree * Tree1c = (TTree*) sourceFile1c->Get("summaryTree");
  TTree * Tree2 = (TTree*) sourceFile2->Get("summaryTree");

  /////////////////////////
  ///// lepton selection
  // for different datasets
  TString lepselection1a = "(TwoMuon && isDoubleMuTriggerPass && PassZmask_met==1 && met_pt>50)"; //Selection for TwoMuon data events
  TString lepselection1b = "(TwoElectron && isDoubleElectronTriggerPass && PassZmask_met==1 && met_pt>50)"; //Selection for TwoEle data events
  TString lepselection1c = "(MuonElectron && isMuEGTriggerPass)"; //Selection for MuonEle data events
  // for MC events
  TString lepselection2 = "( (TwoMuon && isDoubleMuTriggerPass && PassZmask_met==1 && met_pt>50) || (TwoElectron && isDoubleElectronTriggerPass && PassZmask_met==1 && met_pt>50) || (MuonElectron && isMuEGTriggerPass) ) ";

  ///// jet selections
  TString jetselection2a = "(numJets == 2 && jets_by_pt_2_btagCombinedSecVertex > 0.679 "; //Probe jet requirement for first_jet
  jetselection2a += " && jets_by_pt_1_pt >= "+pt_lower+" && jets_by_pt_1_pt < "+pt_upper;
  jetselection2a += " && abs(jets_by_pt_1_eta) >= "+eta_lower+" && abs(jets_by_pt_1_eta) < "+eta_upper+")";
  TString jetselection2b = "(numJets == 2 && jets_by_pt_1_btagCombinedSecVertex > 0.679 "; //Probe jet requirement for second_jet
  jetselection2b += " && jets_by_pt_2_pt >= "+pt_lower+" && jets_by_pt_2_pt < "+pt_upper;
  jetselection2b += " && abs(jets_by_pt_2_eta) >= "+eta_lower+" && abs(jets_by_pt_2_eta) < "+eta_upper+")";

  ///// extra seletions
  TString trigselection = "((dR_leplep > 0.2) && (mass_leplep > 12) && (isCleanEvent == 1) && (oppositeLepCharge == 1))"; //General dilepton selection
  TString firstjetb = "( abs(jets_by_pt_1_flavor)==5 )"; 
  TString secondjetb = "( abs(jets_by_pt_2_flavor)==5 )";

  TString lumi = "19450.0*"; //Nominal luminosity

  /////////// variable  
  TString varName = "jet_CSV_unc"; //Variable which the scale factor is a function of
  TString varName1 = "jets_by_pt_1_btagCombinedSecVertex"; //Variable which the scale factor is a function of
  TString varName2 = "jets_by_pt_2_btagCombinedSecVertex"; //Variable which the scale factor is a function of

  /////////// binning
  double Xmin = 0; //Lower bound for x-axis
  double Xmax = 1.0; //Upper bound for x-axis
  int nBins = 18; //Number of bins
  double xBins[19] = {-10.0, 0.0, 0.122, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.7228, 0.7666, 0.8104, 0.8542, 0.898, 0.9184, 0.9388, 0.9592, 0.9796, 1.01};

  //////////////
  TString canvName = Form("canvas_Pt%i_Eta%i",iPt,iEta);  
  TCanvas * can1a = new TCanvas (canvName, canvName);

  //DoubleMu
  TString proName3 = Form("csv_Data_Pt%i_Eta%i",iPt,iEta);
  TH1F *histTemp3 = new TH1F("histTemp3", "histTemp3", nBins, xBins); //varName1, jetselection2a
  TH1F *histTemp4 = new TH1F("histTemp4", "histTemp4", nBins, xBins); //varName2, jetselection2b
  histTemp3->Sumw2();
  histTemp4->Sumw2();
  TH1 * sourceHisto3 = (TH1*) Tree1a->Draw(varName1+" >> histTemp3","("+lepselection1a+ " && "+trigselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto4 = (TH1*) Tree1a->Draw(varName2+" >> histTemp4","("+lepselection1a+ " && "+trigselection+" && "+jetselection2b+")","");

  //DoubleElectron
  TH1F *histTemp7 = new TH1F("histTemp7", "histTemp7", nBins, xBins); //varName1, jetselection2a
  TH1F *histTemp8 = new TH1F("histTemp8", "histTemp8", nBins, xBins); //varName2, jetselection2b
  histTemp7->Sumw2();
  histTemp8->Sumw2();
  TH1 * sourceHisto7 = (TH1*) Tree1b->Draw(varName1+" >> histTemp7","("+lepselection1b+ " && "+trigselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto8 = (TH1*) Tree1b->Draw(varName2+" >> histTemp8","("+lepselection1b+ " && "+trigselection+" && "+jetselection2b+")","");

  //MuEG
  TH1F *histTemp11 = new TH1F("histTemp11", "histTemp11", nBins, xBins); //varName1, jetselection2a
  TH1F *histTemp12 = new TH1F("histTemp12", "histTemp12", nBins, xBins); //varName2, jetselection2b
  histTemp11->Sumw2();
  histTemp12->Sumw2();
  TH1 * sourceHisto11 = (TH1*) Tree1c->Draw(varName1+" >> histTemp11","("+lepselection1c+ " && "+trigselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto12 = (TH1*) Tree1c->Draw(varName2+" >> histTemp12","("+lepselection1c+ " && "+trigselection+" && "+jetselection2b+")","");

  
  histTemp3->Add(histTemp4); //Add first and second jet distributions in DoubleMu
  histTemp3->Add(histTemp7); //Add first jet DoubleEle 
  histTemp3->Add(histTemp8); //Add second jet DoubleEle
  histTemp3->Add(histTemp11); //Add first jet MuEG
  histTemp3->Add(histTemp12); //Add second jet MuEG
  
  float Data_tagged = histTemp3->Integral(); //Numerator
//   float Data_untagged = histTemp1->Integral(); //Denominator
  if (EfficiencySF) {
    std::cout << "Data tagged: " << Data_tagged << endl;
//     std::cout << "Data untagged: " << Data_untagged << endl;
//     std::cout << "Data tag eff: " << Data_tagged/Data_untagged << endl;
  }
  else std::cout << "Data: " << Data_tagged << endl;
//   if (EfficiencySF) histTemp3->Divide(histTemp3,histTemp1);
  histTemp3->SetMarkerStyle(8);
  histTemp3->SetLineColor(kBlack);
  histTemp3->SetTitle("Two leptons, Z-masked, ==2 jet");
  histTemp3->SetName(proName3);

  //MC
  TH1F *histTemp3a = new TH1F("histTemp3a", "histTemp3a", nBins, xBins);
  TH1F *histTemp4a = new TH1F("histTemp4a", "histTemp4a", nBins, xBins);
  TH1F *histTemp5a = new TH1F("histTemp5a", "histTemp5a", nBins, xBins);
  TH1F *histTemp6a = new TH1F("histTemp6a", "histTemp6a", nBins, xBins);
  histTemp3a->Sumw2();
  histTemp4a->Sumw2();
  histTemp5a->Sumw2();
  histTemp6a->Sumw2();
  TH1 * sourceHisto3a = (TH1*) Tree2->Draw(varName1+" >> histTemp3a",lumi+"weight_PU*topPtWgt*(Xsec/nGen)*("+lepselection2+ " && "+firstjetb + " && "+trigselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto4a = (TH1*) Tree2->Draw(varName2+" >> histTemp4a",lumi+"weight_PU*topPtWgt*(Xsec/nGen)*("+lepselection2+ " && "+secondjetb +" && "+trigselection+" && "+jetselection2b+")","");
  TH1 * sourceHisto5a = (TH1*) Tree2->Draw(varName1+" >> histTemp5a",lumi+"weight_PU*topPtWgt*(Xsec/nGen)*("+lepselection2+ " && !"+firstjetb +" && "+trigselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto6a = (TH1*) Tree2->Draw(varName2+" >> histTemp6a",lumi+"weight_PU*topPtWgt*(Xsec/nGen)*("+lepselection2+ " && !"+secondjetb +" && "+trigselection+" && "+jetselection2b+")","");

  histTemp3a->Add(histTemp4a);
  histTemp5a->Add(histTemp6a);
  double MC_tagged_bjets = histTemp3a->Integral();
  double MC_tagged_nonbjets = histTemp5a->Integral();
  double MC_tagged = MC_tagged_bjets + MC_tagged_nonbjets;
//   float MC_untagged = histTemp1a->Integral();
  if (EfficiencySF) {
    std::cout << "MC tagged: " << MC_tagged << endl;
  }
  else {
    std::cout << "MC: bjets: nonbjets: " << MC_tagged << ", " << MC_tagged_bjets << ", " << MC_tagged_nonbjets << endl;
    std::cout << "" << endl;
    std::cout << "Data/MC : " << Data_tagged/MC_tagged << endl;
    std::cout << "Rescaling MC to match data yields ..." << endl;
  }


  histTemp3a->Scale(Data_tagged/MC_tagged);
  histTemp3a->SetMarkerStyle(8);
  histTemp3a->SetMarkerColor(kRed);
  histTemp3a->SetLineColor(kRed);
  histTemp3a->SetName(Form("csv_MC_bjets_Pt%i_Eta%i",iPt,iEta));

  histTemp5a->Scale(Data_tagged/MC_tagged);
  histTemp5a->SetMarkerStyle(8);
  histTemp5a->SetMarkerColor(kGreen);
  histTemp5a->SetLineColor(kGreen);
  histTemp5a->SetName(Form("csv_MC_nonbjets_Pt%i_Eta%i",iPt,iEta));

  TH1F *hClone = (TH1F*)histTemp3->Clone(Form("csv_ratio_Pt%i_Eta%i",iPt,iEta));
  TH1F *hNonbjets = (TH1F*)histTemp5a->Clone();
  hNonbjets->Scale(-1);
  hClone->Add(hNonbjets);
  TH1F *hsub = (TH1F*)hClone->Clone(Form("csv_sub_Pt%i_Eta%i",iPt,iEta));
  hClone->Divide(hClone,histTemp3a);
  hClone->SetMarkerStyle(8);
  hClone->SetMarkerColor(kBlue);
  hClone->SetLineColor(kBlue);

  TH1F *hratio = (TH1F*)hClone->Clone();
  histTemp3->SetDirectory(outputFile);
  histTemp3a->SetDirectory(outputFile);
  histTemp5a->SetDirectory(outputFile);
  hratio->SetDirectory(outputFile);
  hsub->SetDirectory(outputFile);

  can1a->cd();

  //Overlay shapes
  //    hClone->Rebin(5); 
  hClone->Draw("histE1");

  TString plotName = "_Pt";
  plotName += iPt;
  plotName += "_Eta";
  plotName += iEta;
  //  can1a->Print("csvRWT_"+varName+"_"+plotName +".png");
  can1a->Print(varName+"_"+plotName +".png");



  
    } //end iEta bin 
  } // end iPt bin

  /////
  outputFile->Write();    
  outputFile->Close();  
}
