#include <cmath>
#include <algorithm>

void lightflavorCSVSF(){
  gStyle->SetOptStat(0);

  //// the input you need to change
  int verNum = 0; // ??? different iteration version
  TString MCfile = "/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/allMC.root"; //??? change MC file for different iteration --- 
  TString dataFile1a = "/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/TwoMuon.root";
  TString dataFile1b = "/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/TwoEle.root";
  TString dataFile1c = "/afs/crc.nd.edu/user/w/wluo1/CMSSW_5_3_8_patch1/src/csvReweighting/dataMCRootFiles/MuonEle.root";

  //---------------------
  TString ofName = Form("csv_rwt_lf_v%i.root",verNum);
  TFile * outputFile = new TFile(ofName, "RECREATE");

  // pt eta bins
  TString pt_lower = "0";
  TString pt_upper = "10000";
  TString eta_lower = "0.0";
  TString eta_upper = "2.4";
  vector<TString> ptBins; 
  ptBins.push_back("30");
  ptBins.push_back("40");
  ptBins.push_back("60");
// //   ptBins.push_back("100");
// //   ptBins.push_back("160");
  ptBins.push_back("10000");
  unsigned int nPt = ptBins.size();
  vector<TString> etaBins;
  etaBins.push_back("0.0");
  etaBins.push_back("0.8");
  etaBins.push_back("1.6");
  etaBins.push_back("2.4");
  unsigned int nEta = etaBins.size();

  /////// loop through differet pt and eta bins
  for (unsigned int iPt = 0; iPt < nPt-1; iPt++){
    pt_lower = ptBins[iPt];
    pt_upper = ptBins[iPt+1];
    cout << "pt range is " << pt_lower << " to " << pt_upper << endl;
    for (unsigned int iEta = 0; iEta < nEta-1; iEta++){
      eta_lower = etaBins[iEta];
      eta_upper = etaBins[iEta+1];
      cout << "eta range is " << eta_lower << " to " << eta_upper << endl;

      ///input root files ---- 
  TFile * sourceFile1a = new TFile (dataFile1a);
  TFile * sourceFile1b = new TFile (dataFile1b);
  TFile * sourceFile1c = new TFile (dataFile1c);
  TFile * sourceFile2 = new TFile  (MCfile); 

  /////////////////////////
  ///// lepton selection
  //--- for different datasets
  TString lepselection1a = "(TwoMuon && isDoubleMuTriggerPass && PassZmask_met==0 && met_pt<30 && abs(mass_leplep-91)<10)"; //Selection for TwoMuon data events
  TString lepselection1b = "(TwoElectron && isDoubleElectronTriggerPass && PassZmask_met==0 && met_pt<30 && abs(mass_leplep-91)<10)"; //Selection for TwoEle data events
//   TString lepselection1c = "(MuonEle && isMuEGTriggerPass && TwoMuon && TwoEle)"; //Selection for MuonEle data events
  //--- for MC events
  TString lepselection2 = "( ((TwoMuon && isDoubleMuTriggerPass ) || (TwoElectron && isDoubleElectronTriggerPass )) && PassZmask_met==0 && met_pt<30 && abs(mass_leplep-91)<10) ";

  ///// jet selections
  TString jetselection2a = "(numJets == 2 && jets_by_pt_2_btagCombinedSecVertex < 0.244 "; //Probe jet being the first_jet
  jetselection2a += " && jets_by_pt_1_pt >= "+pt_lower+" && jets_by_pt_1_pt < "+pt_upper;
  jetselection2a += " && abs(jets_by_pt_1_eta) >= "+eta_lower+" && abs(jets_by_pt_1_eta) < "+eta_upper+")";
  TString jetselection2b = "(numJets == 2 && jets_by_pt_1_btagCombinedSecVertex < 0.244 "; //Probe jet being the second_jet
  jetselection2b += " && jets_by_pt_2_pt >= "+pt_lower+" && jets_by_pt_2_pt < "+pt_upper;
  jetselection2b += " && abs(jets_by_pt_2_eta) >= "+eta_lower+" && abs(jets_by_pt_2_eta) < "+eta_upper+")";

  ///// extra seletions
  TString exselection = "((dR_leplep > 0.2) && (mass_leplep > 12) && (isCleanEvent == 1))"; //General dilepton selection // opposite sign applied in treemaking 
  TString firstjetbc = "( abs(jets_by_pt_1_flavour)==5 || abs(jets_by_pt_1_flavour)==4 )"; 
  TString secondjetbc = "( abs(jets_by_pt_2_flavour)==5 || abs(jets_by_pt_2_flavour)==4 )";

  TString lumi = "19450.0*(Xsec/nGen)*"; //Nominal luminosity
  TString wgtStr = "weight_PU*topPtWgt*osTriggerSF*lepIDAndIsoSF*"; // various weights
  if (verNum !=0 ) wgtStr += "csvWgthf*"; // applying hfSFs 

  /////////// variable   ---  can use other btagging algorithm discriminator as well     
  TString varName = "jet_CSV"; //Variable which the scale factor is a function of
  TString varName1 = "jets_by_pt_1_btagCombinedSecVertex"; //Variable which the scale factor is a function of
  TString varName2 = "jets_by_pt_2_btagCombinedSecVertex"; //Variable which the scale factor is a function of

  /////////// CSV binning
//   double Xmin = 0; //Lower bound for x-axis
//   double Xmax = 1.0; //Upper bound for x-axis
  int nBins = 22; //Number of bins
  double xBins[23] = {-10.0, -1, 0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.915, 0.932, 0.949, 0.966, 0.983, 1.01};

  //////////////
  ///// calculate SFs
  //////////////
  
  TString canvName = Form("canvas_Pt%i_Eta%i",iPt,iEta);  
  TCanvas * can1a = new TCanvas (canvName, canvName);

  //DoubleMu
  sourceFile1a->cd();
  TH1F *histTemp3 = new TH1F("histTemp3", "histTemp3", nBins, xBins); //varName1, jetselection2a
  TH1F *histTemp4 = new TH1F("histTemp4", "histTemp4", nBins, xBins); //varName2, jetselection2b
  histTemp3->Sumw2();
  histTemp4->Sumw2();
  TH1 * sourceHisto3 = (TH1*) summaryTree->Draw(varName1+" >> histTemp3","("+lepselection1a+ " && "+exselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto4 = (TH1*) summaryTree->Draw(varName2+" >> histTemp4","("+lepselection1a+ " && "+exselection+" && "+jetselection2b+")","");

  //DoubleElectron
  sourceFile1b->cd();
  TH1F *histTemp7 = new TH1F("histTemp7", "histTemp7", nBins, xBins); //varName1, jetselection2a
  TH1F *histTemp8 = new TH1F("histTemp8", "histTemp8", nBins, xBins); //varName2, jetselection2b
  histTemp7->Sumw2();
  histTemp8->Sumw2();
  TH1 * sourceHisto7 = (TH1*) summaryTree->Draw(varName1+" >> histTemp7","("+lepselection1b+ " && "+exselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto8 = (TH1*) summaryTree->Draw(varName2+" >> histTemp8","("+lepselection1b+ " && "+exselection+" && "+jetselection2b+")","");

  //MuEG
 //  sourceFile1c->cd();
//   TH1F *histTemp11 = new TH1F("histTemp11", "histTemp11", nBins, xBins); //varName1, jetselection2a
//   TH1F *histTemp12 = new TH1F("histTemp12", "histTemp12", nBins, xBins); //varName2, jetselection2b
//   histTemp11->Sumw2();
//   histTemp12->Sumw2();
//   TH1 * sourceHisto11 = (TH1*) summaryTree->Draw(varName1+" >> histTemp11","("+lepselection1c+ " && "+exselection+" && "+jetselection2a+")","");
//   TH1 * sourceHisto12 = (TH1*) summaryTree->Draw(varName2+" >> histTemp12","("+lepselection1c+ " && "+exselection+" && "+jetselection2b+")","");

  
  histTemp3->Add(histTemp4); //Add first and second jet distributions in DoubleMu
  histTemp3->Add(histTemp7); //Add first jet DoubleEle 
  histTemp3->Add(histTemp8); //Add second jet DoubleEle
//   histTemp3->Add(histTemp11); //Add first jet MuEG
//   histTemp3->Add(histTemp12); //Add second jet MuEG
  
  float Data_tagged = histTemp3->Integral();
  std::cout << "Data: " << Data_tagged << endl;

  histTemp3->SetMarkerStyle(8);
  histTemp3->SetLineColor(kBlack);
  histTemp3->SetTitle("Two leptons, Z-masked, ==2 jet");
  TString proName3 = Form("csv_Data_Pt%i_Eta%i",iPt,iEta);
  histTemp3->SetName(proName3);

  //MC
  sourceFile2->cd();
  TH1F *histTemp3a = new TH1F("histTemp3a", "histTemp3a", nBins, xBins);
  TH1F *histTemp4a = new TH1F("histTemp4a", "histTemp4a", nBins, xBins);
  TH1F *histTemp5a = new TH1F("histTemp5a", "histTemp5a", nBins, xBins);
  TH1F *histTemp6a = new TH1F("histTemp6a", "histTemp6a", nBins, xBins);
  histTemp3a->Sumw2();
  histTemp4a->Sumw2();
  histTemp5a->Sumw2();
  histTemp6a->Sumw2();
  TH1 * sourceHisto3a = (TH1*) summaryTree->Draw(varName1+" >> histTemp3a",lumi+wgtStr+"("+lepselection2+ " && "+firstjetbc+" && "+exselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto4a = (TH1*) summaryTree->Draw(varName2+" >> histTemp4a",lumi+wgtStr+"("+lepselection2+ " && "+secondjetbc+" && "+exselection+" && "+jetselection2b+")","");
  TH1 * sourceHisto5a = (TH1*) summaryTree->Draw(varName1+" >> histTemp5a",lumi+wgtStr+"("+lepselection2+ " && !"+firstjetbc+" && "+exselection+" && "+jetselection2a+")","");
  TH1 * sourceHisto6a = (TH1*) summaryTree->Draw(varName2+" >> histTemp6a",lumi+wgtStr+"("+lepselection2+ " && !"+secondjetbc+" && "+exselection+" && "+jetselection2b+")","");


  histTemp3a->Add(histTemp4a); //MC bc
  histTemp5a->Add(histTemp6a); //MC lf
  double MC_tagged_bjets = histTemp3a->Integral();
  double MC_tagged_nonbjets = histTemp5a->Integral();
  double MC_tagged = MC_tagged_bjets + MC_tagged_nonbjets;

  std::cout << "MC: bjets: nonbjets: " << MC_tagged << ", " << MC_tagged_bjets << ", " << MC_tagged_nonbjets << endl;
  std::cout << "" << endl;
  std::cout << "Data/MC : " << Data_tagged/MC_tagged << endl;
  std::cout << "Rescaling MC to match data yields ..." << endl;


  //-----------
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
  TH1F *hbjets = (TH1F*)histTemp3a->Clone();
  hbjets->Scale(-1);
  hClone->Add(hbjets);
  TH1F *hsub = (TH1F*)hClone->Clone(Form("csv_sub_Pt%i_Eta%i",iPt,iEta));
  hClone->Divide(hClone,histTemp5a);
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
  hClone->Draw("histE1");
  TString plotName = Form("lfSF_Pt%i_Eta%i_v%i",iPt,iEta,verNum);
  can1a->Print(varName+plotName +".png"); 



    } //end iEta bin 
  } // end iPt bin

  /////
  outputFile->Write();    
  outputFile->Close();  
  
}
  
