#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TCanvas.h"

#include "Math/MinimizerOptions.h"
#include <TSpline.h>
#include "Math/Interpolator.h"
#define MAXPOINTS 200

//______________________________________________________________________________
void fit_csvSF_HF( bool makePlots=true, bool verbose=false ){


  TH1::SetDefaultSumw2();

  TString dirprefix = "Images/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


  TFile *file = TFile::Open("HistoFiles/csv_rwt_hf_IT_v2.root");
  TFile *file_JESUp = TFile::Open("HistoFiles/csv_rwt_hf_IT_v2_JESUp.root");
  TFile *file_JESDown = TFile::Open("HistoFiles/csv_rwt_hf_IT_v2_JESDown.root");

  TFile *file_mcbtag = TFile::Open("HistoFiles/mc_btag_efficiency_8TeV_53x.root");

  TFile *file_Stats = TFile::Open("HistoFiles/fit_csvSF_HF_statErr.root");
  
  TH1D* hc_ttbar = (TH1D*)file_mcbtag->Get("ttbar_8TeV_hc");
  TH1D* hc_ttbar_jj = (TH1D*)file_mcbtag->Get("ttbar_jj_8TeV_hc");
  TH1D* hc_ttbar_lj = (TH1D*)file_mcbtag->Get("ttbar_lj_8TeV_hc");
  TH1D* hc_ttbar_ll = (TH1D*)file_mcbtag->Get("ttbar_ll_8TeV_hc");

  hc_ttbar->Add(hc_ttbar_jj);
  hc_ttbar->Add(hc_ttbar_lj);
  hc_ttbar->Add(hc_ttbar_ll);

  std::string histofilename = "HistoFiles/csv_rwt_hf_final_IT_v2.root";
  TFile histofile(histofilename.c_str(),"recreate");
  histofile.cd();


  double useUp = 1.2;
  double useDown = 0.8;

  int ncsvbins = 18;
  double csvbins_new[] = { -0.04, 0.0, 0.122, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.723, 0.767, 0.810, 0.854, 0.898, 0.918, 0.939, 0.959, 0.980, 1.01 };


  int ncsvbins_btv = 25;
  double csvbins[] = { 0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.244, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.679, 0.72, 0.76, 0.80, 0.84, 0.898, 0.92, 0.96, 1.001 };

  double btv_csv_sf[] = { 0.46, 0.57, 0.53, 1.28, 1.16, 1.06, 1.77, 1.20, 1.34, 0.90, 1.48, 1.11, 1.42, 1.25, 1.27, 1.12, 1.16, 1.04, 1.27, 0.99, 0.98, 0.91, 0.99, 1.03, 0.97 };
  double btv_csv_sfDown[] = { 0.66, 0.62, 0.43, 0.28, 0.38, 0.25, 0.39, 0.25, 0.27, 0.24, 0.27, 0.17, 0.35, 0.23, 0.28, 0.19, 0.20, 0.24, 0.16, 0.11, 0.15, 0.08, 0.08, 0.04, 0.02 };
  double btv_csv_sfUp[] = { 0.70, 0.39, 0.40, 0.25, 0.38, 0.24, 0.45, 0.19, 0.23, 0.35, 0.33, 0.17, 0.20, 0.18, 0.22, 0.23, 0.17, 0.17, 0.20, 0.11, 0.14, 0.06, 0.06, 0.04, 0.02 };

  double btv_csv_sfDown_x2[ncsvbins_btv];
  double btv_csv_sfUp_x2[ncsvbins_btv];
  for( int iBin=0; iBin<ncsvbins_btv; iBin++ ){
    btv_csv_sfDown_x2[iBin] = 2 * btv_csv_sfDown[iBin];
    btv_csv_sfUp_x2[iBin]   = 2 * btv_csv_sfUp[iBin];
  }

  double csvbins_center[ncsvbins_btv];
  double csvbins_width[ncsvbins_btv];
  for( int iBin=0; iBin<ncsvbins_btv; iBin++ ){
    csvbins_center[iBin] = 0.5 * ( csvbins[iBin] + csvbins[iBin+1] );
    csvbins_width[iBin] = 0.5 * ( csvbins[iBin+1] - csvbins[iBin] );
  }

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(ncsvbins_btv,csvbins_center,btv_csv_sf,csvbins_width,csvbins_width,btv_csv_sfUp,btv_csv_sfDown); 
  gr->SetMarkerSize(1.2);
  gr->SetMarkerStyle(22);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kRed-9);
  gr->SetFillColor(kRed-9);

  TGraphAsymmErrors *gr_x2 = new TGraphAsymmErrors(ncsvbins_btv,csvbins_center,btv_csv_sf,csvbins_width,csvbins_width,btv_csv_sfUp_x2,btv_csv_sfDown_x2); 
  gr_x2->SetMarkerSize(1.2);
  gr_x2->SetMarkerStyle(22);
  gr_x2->SetMarkerColor(kBlack);
  gr_x2->SetLineColor(kRed-9);
  gr_x2->SetFillColor(kRed-9);

  std::vector<TString> hist_name;
  std::vector<TString> data_hist_name;
  std::vector<TString> mc_b_hist_name;
  std::vector<TString> mc_nonb_hist_name;

  int maxPt  = 5;
  int maxEta = 1;

  for( int iPt=0; iPt<maxPt; iPt++ ){
    for( int iEta=0; iEta<maxEta; iEta++ ){
      hist_name.push_back( Form("csv_ratio_Pt%d_Eta%d", iPt, iEta) );

      data_hist_name.push_back( Form("csv_Data_Pt%d_Eta%d", iPt, iEta) );
      mc_b_hist_name.push_back( Form("csv_MC_bjets_Pt%d_Eta%d", iPt, iEta) );
      mc_nonb_hist_name.push_back( Form("csv_MC_nonbjets_Pt%d_Eta%d", iPt, iEta) );
    }
  }

  int NumHists_normal = int( hist_name.size() );
  int numHists = NumHists_normal+2;

  TH1D* h_csv_ratio[numHists];
  TH1D* h_csv_ratio_LF[numHists];
  TH1D* h_csv_ratio_LFUp[numHists];
  TH1D* h_csv_ratio_LFDown[numHists];
  TH1D* h_csv_ratio_JESUp[numHists];
  TH1D* h_csv_ratio_JESDown[numHists];
  TH1D* h_csv_ratio_Stats1Up[numHists];
  TH1D* h_csv_ratio_Stats1Down[numHists];
  TH1D* h_csv_ratio_Stats2Up[numHists];
  TH1D* h_csv_ratio_Stats2Down[numHists];

  TH1D* h_csv_data[NumHists_normal];
  TH1D* h_csv_mc_b[NumHists_normal];
  TH1D* h_csv_mc_nonb[NumHists_normal];

  TH1D* h_csv_mc_b_JESUp[NumHists_normal];
  TH1D* h_csv_mc_nonb_JESUp[NumHists_normal];
  TH1D* h_csv_mc_b_JESDown[NumHists_normal];
  TH1D* h_csv_mc_nonb_JESDown[NumHists_normal];


  TH1D* h_csv_ratio_StatsUp[numHists][ncsvbins];
  TH1D* h_csv_ratio_StatsDown[numHists][ncsvbins];



  TH1D* h_csv_data_all = NULL;
  TH1D* h_csv_mc_b_all = NULL;
  TH1D* h_csv_mc_nonb_all = NULL;


  TH1D* h_csv_mc_b_all_JESUp = NULL;
  TH1D* h_csv_mc_nonb_all_JESUp = NULL;

  TH1D* h_csv_mc_b_all_JESDown = NULL;
  TH1D* h_csv_mc_nonb_all_JESDown = NULL;

  for( int iHist=0; iHist<NumHists_normal; iHist++ ){
    TH1D* h_csv_data_temp0 = (TH1D*)file->Get( data_hist_name[iHist] )->Clone( Form("h_%s_temp0",data_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_b_temp0 = (TH1D*)file->Get( mc_b_hist_name[iHist] )->Clone( Form("h_%s_temp0",mc_b_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_nonb_temp0 = (TH1D*)file->Get( mc_nonb_hist_name[iHist] )->Clone( Form("h_%s_temp0",mc_nonb_hist_name[iHist].Data()) );

    // JES up/down
    TH1D* h_csv_mc_b_JESUp_temp0 = (TH1D*)file_JESUp->Get( mc_b_hist_name[iHist] )->Clone( Form("h_%s_JESUp_temp0",mc_b_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_nonb_JESUp_temp0 = (TH1D*)file_JESUp->Get( mc_nonb_hist_name[iHist] )->Clone( Form("h_%s_JESUp_temp0",mc_nonb_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_b_JESDown_temp0 = (TH1D*)file_JESDown->Get( mc_b_hist_name[iHist] )->Clone( Form("h_%s_JESDown_temp0",mc_b_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_nonb_JESDown_temp0 = (TH1D*)file_JESDown->Get( mc_nonb_hist_name[iHist] )->Clone( Form("h_%s_JESDown_temp0",mc_nonb_hist_name[iHist].Data()) );

    h_csv_data[iHist] = new TH1D( Form("h_%s",data_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_b[iHist] = new TH1D( Form("h_%s",mc_b_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_nonb[iHist] = new TH1D( Form("h_%s",mc_nonb_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );

    h_csv_mc_b_JESUp[iHist] = new TH1D( Form("h_%s_JESUp",mc_b_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_nonb_JESUp[iHist] = new TH1D( Form("h_%s_JESUp",mc_nonb_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_b_JESDown[iHist] = new TH1D( Form("h_%s_JESDown",mc_b_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_nonb_JESDown[iHist] = new TH1D( Form("h_%s_JESDown",mc_nonb_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );

    for( int iBin=0; iBin<ncsvbins; iBin++ ){
      h_csv_data[iHist]->SetBinContent(iBin+1, h_csv_data_temp0->GetBinContent(iBin+1));
      h_csv_data[iHist]->SetBinError(iBin+1, h_csv_data_temp0->GetBinError(iBin+1));

      h_csv_mc_b[iHist]->SetBinContent(iBin+1, h_csv_mc_b_temp0->GetBinContent(iBin+1));
      h_csv_mc_b[iHist]->SetBinError(iBin+1, h_csv_mc_b_temp0->GetBinError(iBin+1));

      h_csv_mc_nonb[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_temp0->GetBinContent(iBin+1));
      h_csv_mc_nonb[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_temp0->GetBinError(iBin+1));

      h_csv_mc_b_JESUp[iHist]->SetBinContent(iBin+1, h_csv_mc_b_JESUp_temp0->GetBinContent(iBin+1));
      h_csv_mc_b_JESUp[iHist]->SetBinError(iBin+1, h_csv_mc_b_JESUp_temp0->GetBinError(iBin+1));
      h_csv_mc_nonb_JESUp[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_JESUp_temp0->GetBinContent(iBin+1));
      h_csv_mc_nonb_JESUp[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_JESUp_temp0->GetBinError(iBin+1));

      h_csv_mc_b_JESDown[iHist]->SetBinContent(iBin+1, h_csv_mc_b_JESDown_temp0->GetBinContent(iBin+1));
      h_csv_mc_b_JESDown[iHist]->SetBinError(iBin+1, h_csv_mc_b_JESDown_temp0->GetBinError(iBin+1));
      h_csv_mc_nonb_JESDown[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_JESDown_temp0->GetBinContent(iBin+1));
      h_csv_mc_nonb_JESDown[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_JESDown_temp0->GetBinError(iBin+1));
    }


    if( iHist==0 ){
      h_csv_data_all = (TH1D*)h_csv_data[iHist]->Clone("h_csv_data_all");
      h_csv_mc_b_all = (TH1D*)h_csv_mc_b[iHist]->Clone("h_csv_mc_b_all");
      h_csv_mc_nonb_all = (TH1D*)h_csv_mc_nonb[iHist]->Clone("h_csv_mc_nonb_all");

      h_csv_mc_b_all_JESUp = (TH1D*)h_csv_mc_b_JESUp[iHist]->Clone("h_csv_mc_b_all_JESUp");
      h_csv_mc_nonb_all_JESUp = (TH1D*)h_csv_mc_nonb_JESUp[iHist]->Clone("h_csv_mc_nonb_all_JESUp");
      h_csv_mc_b_all_JESDown = (TH1D*)h_csv_mc_b_JESDown[iHist]->Clone("h_csv_mc_b_all_JESDown");
      h_csv_mc_nonb_all_JESDown = (TH1D*)h_csv_mc_nonb_JESDown[iHist]->Clone("h_csv_mc_nonb_all_JESDown");
    }
    else{
      h_csv_data_all->Add(h_csv_data[iHist]);
      h_csv_mc_b_all->Add(h_csv_mc_b[iHist]);
      h_csv_mc_nonb_all->Add(h_csv_mc_nonb[iHist]);

      h_csv_mc_b_all_JESUp->Add(h_csv_mc_b_JESUp[iHist]);
      h_csv_mc_nonb_all_JESUp->Add(h_csv_mc_nonb_JESUp[iHist]);
      h_csv_mc_b_all_JESDown->Add(h_csv_mc_b_JESDown[iHist]);
      h_csv_mc_nonb_all_JESDown->Add(h_csv_mc_nonb_JESDown[iHist]);
    }


    h_csv_ratio[iHist]        = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_%d",iHist));
    h_csv_ratio_LFUp[iHist]   = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_LFUp_%d",iHist));
    h_csv_ratio_LFDown[iHist] = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_LFDown_%d",iHist));

    h_csv_ratio_JESUp[iHist]   = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_JESUp_%d",iHist));
    h_csv_ratio_JESDown[iHist] = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_JESDown_%d",iHist));

    TH1D* h_csv_mc_b_temp0_LFUp = (TH1D*)h_csv_mc_b[iHist]->Clone(Form("h_csv_mc_b_temp0_LFUp_%d",iHist));
    TH1D* h_csv_mc_b_temp0_LFDown = (TH1D*)h_csv_mc_b[iHist]->Clone(Form("h_csv_mc_b_temp0_LFDown_%d",iHist));

    TH1D* h_csv_mc_nonb_temp0_LFUp = (TH1D*)h_csv_mc_nonb[iHist]->Clone(Form("h_csv_mc_nonb_temp0_LFUp_%d",iHist));
    TH1D* h_csv_mc_nonb_temp0_LFDown = (TH1D*)h_csv_mc_nonb[iHist]->Clone(Form("h_csv_mc_nonb_temp0_LFDown_%d",iHist));

    h_csv_mc_b_temp0_LFUp->Scale( h_csv_data[iHist]->Integral() / ( h_csv_mc_b[iHist]->Integral() + useUp*h_csv_mc_nonb[iHist]->Integral() ) );
    h_csv_mc_b_temp0_LFDown->Scale( h_csv_data[iHist]->Integral() / ( h_csv_mc_b[iHist]->Integral() + useDown*h_csv_mc_nonb[iHist]->Integral() ) );

    h_csv_mc_nonb_temp0_LFUp->Scale( h_csv_data[iHist]->Integral() / ( h_csv_mc_b[iHist]->Integral() + useUp*h_csv_mc_nonb[iHist]->Integral() ) );
    h_csv_mc_nonb_temp0_LFDown->Scale( h_csv_data[iHist]->Integral() / ( h_csv_mc_b[iHist]->Integral() + useDown*h_csv_mc_nonb[iHist]->Integral() ) );

    h_csv_ratio[iHist]->Add(h_csv_mc_nonb[iHist],-1);
    h_csv_ratio_LFUp[iHist]->Add(h_csv_mc_nonb_temp0_LFUp,-useUp);
    h_csv_ratio_LFDown[iHist]->Add(h_csv_mc_nonb_temp0_LFDown,-useDown);
    h_csv_ratio_JESUp[iHist]->Add(h_csv_mc_nonb_JESUp[iHist],-1);
    h_csv_ratio_JESDown[iHist]->Add(h_csv_mc_nonb_JESDown[iHist],-1);

    h_csv_ratio[iHist]->Divide(h_csv_mc_b[iHist]);
    h_csv_ratio_LFUp[iHist]->Divide(h_csv_mc_b_temp0_LFUp);
    h_csv_ratio_LFDown[iHist]->Divide(h_csv_mc_b_temp0_LFDown);
    h_csv_ratio_JESUp[iHist]->Divide(h_csv_mc_b_JESUp[iHist]);
    h_csv_ratio_JESDown[iHist]->Divide(h_csv_mc_b_JESDown[iHist]);


    for( int iBin=0; iBin<ncsvbins; iBin++ ){
      h_csv_ratio_StatsUp[iHist][iBin] = (TH1D*)file_Stats->Get( Form("h_%s_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1) )->Clone( Form("h_%s_CSVbin%d_StatsUp_temp",hist_name[iHist].Data(),iBin+1) );
      h_csv_ratio_StatsDown[iHist][iBin] = (TH1D*)file_Stats->Get( Form("h_%s_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1) )->Clone( Form("h_%s_CSVbin%d_StatsDown_temp",hist_name[iHist].Data(),iBin+1) );
    }

  }


  TH1D* h_csv_ratio_all = (TH1D*)h_csv_data_all->Clone("h_csv_ratio_all_temp");
  TH1D* h_csv_ratio_all_LFUp   = (TH1D*)h_csv_data_all->Clone("h_csv_ratio_all_LFUp_temp");
  TH1D* h_csv_ratio_all_LFDown = (TH1D*)h_csv_data_all->Clone("h_csv_ratio_all_LFDown_temp");

  TH1D* h_csv_ratio_all_JESUp   = (TH1D*)h_csv_data_all->Clone("h_csv_ratio_all_JESUp_temp");
  TH1D* h_csv_ratio_all_JESDown = (TH1D*)h_csv_data_all->Clone("h_csv_ratio_all_JESDown_temp");

  TH1D* h_csv_ratio_cumulative = (TH1D*)h_csv_ratio_all->Clone("h_csv_ratio_cumulative");
  TH1D* h_csv_ratio_cumulative_LFUp = (TH1D*)h_csv_ratio_all->Clone("h_csv_ratio_cumulative_LFUp");
  TH1D* h_csv_ratio_cumulative_LFDown = (TH1D*)h_csv_ratio_all->Clone("h_csv_ratio_cumulative_LFDown");
  TH1D* h_csv_ratio_cumulative_JESUp = (TH1D*)h_csv_ratio_all->Clone("h_csv_ratio_cumulative_JESUp");
  TH1D* h_csv_ratio_cumulative_JESDown = (TH1D*)h_csv_ratio_all->Clone("h_csv_ratio_cumulative_JESDown");

  int nBins = h_csv_ratio_cumulative->GetNbinsX();

  TH1D* h_csv_mc_b_all_LFUp = (TH1D*)h_csv_mc_b_all->Clone("h_csv_mc_b_all_LFUp");
  TH1D* h_csv_mc_b_all_LFDown = (TH1D*)h_csv_mc_b_all->Clone("h_csv_mc_b_all_LFDown");

  TH1D* h_csv_mc_nonb_all_LFUp = (TH1D*)h_csv_mc_nonb_all->Clone("h_csv_mc_nonb_all_LFUp");
  TH1D* h_csv_mc_nonb_all_LFDown = (TH1D*)h_csv_mc_nonb_all->Clone("h_csv_mc_nonb_all_LFDown");


  h_csv_mc_b_all_LFUp->Scale( h_csv_ratio_all->Integral() / (h_csv_mc_b_all->Integral() +  useUp*h_csv_mc_nonb_all->Integral()) );
  h_csv_mc_b_all_LFDown->Scale( h_csv_ratio_all->Integral() / (h_csv_mc_b_all->Integral() +  useDown*h_csv_mc_nonb_all->Integral()) );

  h_csv_mc_nonb_all_LFUp->Scale( useUp * h_csv_ratio_all->Integral() / (h_csv_mc_b_all->Integral() +  useUp*h_csv_mc_nonb_all->Integral()) );
  h_csv_mc_nonb_all_LFDown->Scale( useDown * h_csv_ratio_all->Integral() / (h_csv_mc_b_all->Integral() +  useDown*h_csv_mc_nonb_all->Integral()) );

  h_csv_ratio_all->Add(h_csv_mc_nonb_all,-1);
  h_csv_ratio_all_LFUp->Add(h_csv_mc_nonb_all_LFUp,-1);
  h_csv_ratio_all_LFDown->Add(h_csv_mc_nonb_all_LFDown,-1);
  h_csv_ratio_all_JESUp->Add(h_csv_mc_nonb_all_JESUp,-1);
  h_csv_ratio_all_JESDown->Add(h_csv_mc_nonb_all_JESDown,-1);


  for( int iBin=0; iBin<nBins; iBin++ ){
    h_csv_ratio_cumulative->SetBinContent( iBin+1, h_csv_ratio_all->Integral(iBin+1,nBins) );
    h_csv_ratio_cumulative_LFUp->SetBinContent( iBin+1, h_csv_ratio_all_LFUp->Integral(iBin+1,nBins) );
    h_csv_ratio_cumulative_LFDown->SetBinContent( iBin+1, h_csv_ratio_all_LFDown->Integral(iBin+1,nBins) );
    h_csv_ratio_cumulative_JESUp->SetBinContent( iBin+1, h_csv_ratio_all_JESUp->Integral(iBin+1,nBins) );
    h_csv_ratio_cumulative_JESDown->SetBinContent( iBin+1, h_csv_ratio_all_JESDown->Integral(iBin+1,nBins) );
  }


  h_csv_ratio_all->Divide(h_csv_mc_b_all);
  h_csv_ratio_all_LFUp->Divide(h_csv_mc_b_all_LFUp);
  h_csv_ratio_all_LFDown->Divide(h_csv_mc_b_all_LFDown);
  h_csv_ratio_all_JESUp->Divide(h_csv_mc_b_all_JESUp);
  h_csv_ratio_all_JESDown->Divide(h_csv_mc_b_all_JESDown);



  TH1D* h_mc_nonb_cumulative = (TH1D*)h_csv_mc_nonb_all->Clone("h_mc_nonb_cumulative");
  TH1D* h_mc_b_cumulative = (TH1D*)h_csv_mc_b_all->Clone("h_mc_b_cumulative");
  TH1D* h_mc_b_cumulative_LFUp = (TH1D*)h_csv_mc_b_all_LFUp->Clone("h_mc_b_cumulative_LFUp");
  TH1D* h_mc_b_cumulative_LFDown = (TH1D*)h_csv_mc_b_all_LFDown->Clone("h_mc_b_cumulative_LFDown");
  TH1D* h_mc_b_cumulative_JESUp = (TH1D*)h_csv_mc_b_all_JESUp->Clone("h_mc_b_cumulative_JESUp");
  TH1D* h_mc_b_cumulative_JESDown = (TH1D*)h_csv_mc_b_all_JESDown->Clone("h_mc_b_cumulative_JESDown");

  for( int iBin=0; iBin<nBins; iBin++ ){
    h_mc_nonb_cumulative->SetBinContent( iBin+1, h_csv_mc_nonb_all->Integral(iBin+1,nBins) );

    h_mc_b_cumulative->SetBinContent( iBin+1, h_csv_mc_b_all->Integral(iBin+1,nBins) );
    h_mc_b_cumulative_LFUp->SetBinContent( iBin+1, h_csv_mc_b_all_LFUp->Integral(iBin+1,nBins) );
    h_mc_b_cumulative_LFDown->SetBinContent( iBin+1, h_csv_mc_b_all_LFDown->Integral(iBin+1,nBins) );
    h_mc_b_cumulative_JESUp->SetBinContent( iBin+1, h_csv_mc_b_all_JESUp->Integral(iBin+1,nBins) );
    h_mc_b_cumulative_JESDown->SetBinContent( iBin+1, h_csv_mc_b_all_JESDown->Integral(iBin+1,nBins) );
  }

  h_csv_ratio_cumulative->Divide(h_mc_b_cumulative);
  h_csv_ratio_cumulative_LFUp->Divide(h_mc_b_cumulative_LFUp);
  h_csv_ratio_cumulative_LFDown->Divide(h_mc_b_cumulative_LFDown);
  h_csv_ratio_cumulative_JESUp->Divide(h_mc_b_cumulative_JESUp);
  h_csv_ratio_cumulative_JESDown->Divide(h_mc_b_cumulative_JESDown);


  for( int iBin=0; iBin<nBins; iBin++ ){
    if( verbose ) {
      printf("\t iBin=%d,\t D=%.0f,\t LF=%.0f,\t HF=%.0f,\t iD=%.0f,\t iLF=%.0f,\t iHF=%.0f,\t iLFUp=%.0f,\t iHFUp=%.0f,\t iLFDown=%.0f,\t iHFDown=%.0f,\t (D-LF)/HF=%.2f,\t i(D-LF)/HF=%.3f,\t LFUp=%.3f,\t LFDown=%.3f,\t cum=%.2f,\t cumLFUp=%.2f,\t cumLFDown=%.2f \n", 
	   iBin, h_csv_data_all->GetBinContent(iBin+1), h_csv_mc_nonb_all->GetBinContent(iBin+1), h_csv_mc_b_all->GetBinContent(iBin+1),
	   h_csv_data_all->Integral(iBin+1,nBins), h_csv_mc_nonb_all->Integral(iBin+1,nBins), h_csv_mc_b_all->Integral(iBin+1,nBins),
	   h_csv_mc_nonb_all_LFUp->Integral(iBin+1,nBins), h_csv_mc_b_all_LFUp->Integral(iBin+1,nBins),
	   h_csv_mc_nonb_all_LFDown->Integral(iBin+1,nBins), h_csv_mc_b_all_LFDown->Integral(iBin+1,nBins),
	   (h_csv_data_all->GetBinContent(iBin+1) - h_csv_mc_nonb_all->GetBinContent(iBin+1))/h_csv_mc_b_all->GetBinContent(iBin+1),
	   (h_csv_data_all->Integral(iBin+1,nBins) - h_csv_mc_nonb_all->Integral(iBin+1,nBins))/h_csv_mc_b_all->Integral(iBin+1,nBins),
	   (h_csv_data_all->Integral(iBin+1,nBins) - h_csv_mc_nonb_all_LFUp->Integral(iBin+1,nBins))/h_csv_mc_b_all_LFUp->Integral(iBin+1,nBins),
	   (h_csv_data_all->Integral(iBin+1,nBins) - h_csv_mc_nonb_all_LFDown->Integral(iBin+1,nBins))/h_csv_mc_b_all_LFDown->Integral(iBin+1,nBins),
	   h_csv_ratio_cumulative->GetBinContent(iBin+1), h_csv_ratio_cumulative_LFUp->GetBinContent(iBin+1), h_csv_ratio_cumulative_LFDown->GetBinContent(iBin+1) );
    }
  }

  double csvL = h_csv_ratio_cumulative->GetBinContent(4);
  double csvL_LFUp = h_csv_ratio_cumulative_LFUp->GetBinContent(4) - csvL;
  double csvL_JESUp = h_csv_ratio_cumulative_JESUp->GetBinContent(4) - csvL;
  double csvL_LFDown = h_csv_ratio_cumulative_LFDown->GetBinContent(4) - csvL;
  double csvL_JESDown = h_csv_ratio_cumulative_JESDown->GetBinContent(4) - csvL;

  double csvL_Up = sqrt( csvL_LFUp*csvL_LFUp + csvL_JESUp*csvL_JESUp );
  double csvL_Down = sqrt( csvL_LFDown*csvL_LFDown + csvL_JESDown*csvL_JESDown );

  double csvL_Err = 0.5 * (csvL_Up + csvL_Down);

  if( verbose ) printf("\t CSVL: %.3f +/- %.3f \n", csvL, csvL_Err );

  double csvM = h_csv_ratio_cumulative->GetBinContent(9);
  double csvM_LFUp = h_csv_ratio_cumulative_LFUp->GetBinContent(9) - csvM;
  double csvM_JESUp = h_csv_ratio_cumulative_JESUp->GetBinContent(9) - csvM;
  double csvM_LFDown = h_csv_ratio_cumulative_LFDown->GetBinContent(9) - csvM;
  double csvM_JESDown = h_csv_ratio_cumulative_JESDown->GetBinContent(9) - csvM;

  double csvM_Up = sqrt( csvM_LFUp*csvM_LFUp + csvM_JESUp*csvM_JESUp );
  double csvM_Down = sqrt( csvM_LFDown*csvM_LFDown + csvM_JESDown*csvM_JESDown );

  double csvM_Err = 0.5 * (csvM_Up + csvM_Down);

  if( verbose ) printf("\t CSVL: %.3f +/- %.3f \n", csvM, csvM_Err );

  double csvT = h_csv_ratio_cumulative->GetBinContent(14);
  double csvT_LFUp = h_csv_ratio_cumulative_LFUp->GetBinContent(14) - csvT;
  double csvT_JESUp = h_csv_ratio_cumulative_JESUp->GetBinContent(14) - csvT;
  double csvT_LFDown = h_csv_ratio_cumulative_LFDown->GetBinContent(14) - csvT;
  double csvT_JESDown = h_csv_ratio_cumulative_JESDown->GetBinContent(14) - csvT;

  double csvT_Up = sqrt( csvT_LFUp*csvT_LFUp + csvT_JESUp*csvT_JESUp );
  double csvT_Down = sqrt( csvT_LFDown*csvT_LFDown + csvT_JESDown*csvT_JESDown );

  double csvT_Err = 0.5 * (csvT_Up + csvT_Down);

  if( verbose ) printf("\t CSVL: %.3f +/- %.3f \n", csvT, csvT_Err );


  hist_name.push_back("csv_ratio_all");

  hist_name.push_back("csv_ratio_cumulative_all");


  h_csv_ratio[numHists-2] = (TH1D*)h_csv_ratio_all->Clone( Form("h_%s",hist_name[numHists-2].Data()) );

  h_csv_ratio_LFUp[numHists-2] = (TH1D*)h_csv_ratio_all_LFUp->Clone( Form("h_%s_LFUp",hist_name[numHists-2].Data()) );
  h_csv_ratio_LFDown[numHists-2] = (TH1D*)h_csv_ratio_all_LFDown->Clone( Form("h_%s_LFDown",hist_name[numHists-2].Data()) );

  h_csv_ratio_JESUp[numHists-2] = (TH1D*)h_csv_ratio_all_JESUp->Clone( Form("h_%s_JESUp",hist_name[numHists-2].Data()) );
  h_csv_ratio_JESDown[numHists-2] = (TH1D*)h_csv_ratio_all_JESDown->Clone( Form("h_%s_JESDown",hist_name[numHists-2].Data()) );


  h_csv_ratio[numHists-1] = (TH1D*)h_csv_ratio_cumulative->Clone( Form("h_%s",hist_name[numHists-1].Data()) );
  h_csv_ratio[numHists-1]->SetMarkerStyle(20);


  for( int iHist=0; iHist<numHists-1; iHist++ ){
    h_csv_ratio_LF[iHist] = (TH1D*)h_csv_ratio[iHist]->Clone( Form("h_%s_LF",hist_name[numHists-2].Data()) );

    for( int iBin=0; iBin<nBins; iBin++ ){
      double content = h_csv_ratio_LF[iHist]->GetBinContent(iBin+1);
      double hfUp = fabs( content - h_csv_ratio_LFUp[iHist]->GetBinContent(iBin+1) );
      double hfDown = fabs( content - h_csv_ratio_LFDown[iHist]->GetBinContent(iBin+1) );

      double hfUnc = 0.5 * ( hfUp + hfDown );
      h_csv_ratio_LF[iHist]->SetBinError(iBin+1,hfUnc);
    }
  }


  for( int iHist=0; iHist<numHists-1; iHist++ ){
    h_csv_ratio_Stats1Up[iHist] = (TH1D*)h_csv_ratio[iHist]->Clone( Form("h_%s_Stats1Up",hist_name[numHists-2].Data()) );
    h_csv_ratio_Stats1Down[iHist] = (TH1D*)h_csv_ratio[iHist]->Clone( Form("h_%s_Stats1Down",hist_name[numHists-2].Data()) );
    h_csv_ratio_Stats2Up[iHist] = (TH1D*)h_csv_ratio[iHist]->Clone( Form("h_%s_Stats2Up",hist_name[numHists-2].Data()) );
    h_csv_ratio_Stats2Down[iHist] = (TH1D*)h_csv_ratio[iHist]->Clone( Form("h_%s_Stats2Down",hist_name[numHists-2].Data()) );

    for( int iBin=0; iBin<nBins; iBin++ ){
      double center = h_csv_ratio[iHist]->GetBinCenter(iBin+1);
      double content = h_csv_ratio[iHist]->GetBinContent(iBin+1);

      double delta = h_csv_ratio[iHist]->GetBinError(iBin+1);
      double stat1Up   = content + delta * ( 1 - 2*center );
      double stat1Down = content + delta * ( 2*center - 1 );

      double stat2Up   = content + delta * ( 1 - 6*center*(1-center) );
      double stat2Down = content - delta * ( 1 - 6*center*(1-center) );

      h_csv_ratio_Stats1Up[iHist]->SetBinContent(iBin+1,stat1Up);
      h_csv_ratio_Stats1Down[iHist]->SetBinContent(iBin+1,stat1Down);

      h_csv_ratio_Stats2Up[iHist]->SetBinContent(iBin+1,stat2Up);
      h_csv_ratio_Stats2Down[iHist]->SetBinContent(iBin+1,stat2Down);
    }
  }



  int n = 100000;
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(n);



  TH1D* h_csv_ratio_final[numHists];
  TH1D* h_csv_ratio_final_LFUp[numHists];
  TH1D* h_csv_ratio_final_LFDown[numHists];
  TH1D* h_csv_ratio_final_JESUp[numHists];
  TH1D* h_csv_ratio_final_JESDown[numHists];
  TH1D* h_csv_ratio_final_Stats1Up[numHists];
  TH1D* h_csv_ratio_final_Stats1Down[numHists];
  TH1D* h_csv_ratio_final_Stats2Up[numHists];
  TH1D* h_csv_ratio_final_Stats2Down[numHists];

  TH1D* c_csv_ratio_final[numHists];


  TH1D* h_csv_ratio_final_StatsUp[numHists][ncsvbins];
  TH1D* h_csv_ratio_final_StatsDown[numHists][ncsvbins];

  int NumFinalBins = 1000;


  TCanvas *c1 = new TCanvas("c1");
  TF1* btagSF = new TF1("btagSF","0.607239*x*x*x*x+-1.49984*x*x*x+1.3473*x*x+-0.599888*x+1.09396",0.15,0.95);
  TF1* btagSFUp   = new TF1("btagSFUp","1.02*(0.607239*x*x*x*x+-1.49984*x*x*x+1.3473*x*x+-0.599888*x+1.09396)",0.15,0.95);
  TF1* btagSFDown = new TF1("btagSFDown","0.98*(0.607239*x*x*x*x+-1.49984*x*x*x+1.3473*x*x+-0.599888*x+1.09396)",0.15,0.95);
  h_csv_ratio_cumulative->SetStats(0);
  h_csv_ratio_cumulative->GetYaxis()->SetRangeUser(0.88,1.1);
  h_csv_ratio_cumulative->GetXaxis()->SetRangeUser(0.1,1.00);
  h_csv_ratio_cumulative->SetMarkerStyle(20);
  h_csv_ratio_cumulative->SetTitle(";CSV;Heavy flavor scale factor");

  TF1* btagSFUp_x2 = new TF1("btagSFUp_x2","2 * (btagSFUp-btagSF) + btagSF");
  TF1* btagSFDown_x2 = new TF1("btagSFDown_x2","-2 * (btagSF-btagSFDown) + btagSF");

  btagSFUp->SetLineStyle(2);
  btagSFDown->SetLineStyle(2);

  btagSFUp->SetLineColor(kMagenta-2);
  btagSFDown->SetLineColor(kMagenta-2);



  btagSFUp_x2->SetLineStyle(2);
  btagSFDown_x2->SetLineStyle(2);

  btagSFUp_x2->SetLineColor(kMagenta-2);
  btagSFDown_x2->SetLineColor(kMagenta-2);


  //TH1D* h_csv_ratio_cumulative_LF = (TH1D*)h_csv_ratio_cumulative_LFUp->Clone("h_csv_ratio_cumulative_LF");

  for( int iBin=0; iBin<nBins; iBin++ ){
    double center = h_csv_ratio_cumulative->GetBinCenter(iBin+1);
    if( center>0.95 || center<0.15 ){
      h_csv_ratio_cumulative->SetBinContent( iBin+1, 0 );
      h_csv_ratio_cumulative_LFUp->SetBinContent( iBin+1, 0 );
      h_csv_ratio_cumulative_LFDown->SetBinContent( iBin+1, 0 );
      h_csv_ratio_cumulative_JESUp->SetBinContent( iBin+1, 0 );
      h_csv_ratio_cumulative_JESDown->SetBinContent( iBin+1, 0 );
    }
  }


  h_csv_ratio_cumulative_LFUp->SetLineColor(kGreen+3);
  h_csv_ratio_cumulative_LFDown->SetLineColor(kGreen+3);
  h_csv_ratio_cumulative_JESUp->SetLineColor(kRed);
  h_csv_ratio_cumulative_JESDown->SetLineColor(kBlue);

  h_csv_ratio_cumulative_LFUp->SetLineWidth(2);
  h_csv_ratio_cumulative_LFDown->SetLineWidth(2);
  h_csv_ratio_cumulative_JESUp->SetLineWidth(2);
  h_csv_ratio_cumulative_JESDown->SetLineWidth(2);


  h_csv_ratio_cumulative->Draw("pe1");
  btagSF->Draw("same");
  btagSFUp->Draw("same");
  btagSFDown->Draw("same");

  h_csv_ratio_cumulative_LFUp->Draw("][histsame");
  h_csv_ratio_cumulative_LFDown->Draw("][histsame");
  h_csv_ratio_cumulative_JESUp->Draw("][histsame");
  h_csv_ratio_cumulative_JESDown->Draw("][histsame");
  
  TLegend *legend_compare_btv = new TLegend(0.2,0.75,0.87,0.89);

  legend_compare_btv->SetFillColor(kWhite);
  legend_compare_btv->SetLineColor(kWhite);
  legend_compare_btv->SetShadowColor(kWhite);
  legend_compare_btv->SetTextFont(42);
  legend_compare_btv->SetTextSize(0.05);

  legend_compare_btv->SetNColumns(3);

  legend_compare_btv->AddEntry(h_csv_ratio_cumulative,"Scale Factor","pe1");
  legend_compare_btv->AddEntry(h_csv_ratio_cumulative_JESUp,"JES Up","l");
  legend_compare_btv->AddEntry(h_csv_ratio_cumulative_JESDown,"JES Down","l");
  legend_compare_btv->AddEntry(h_csv_ratio_cumulative_LFUp,"LF Err","l");
  legend_compare_btv->AddEntry(btagSF,"BTAG POG","l");
  legend_compare_btv->AddEntry(btagSFUp,"BTAG POG Err","l");

  legend_compare_btv->Draw();

  TString myimg = dirprefix + "hfSF_cumulative.png";
  if( makePlots ) c1->Print(myimg);

  myimg = dirprefix + "hfSF_cumulative.pdf";
  if( makePlots ) c1->Print(myimg);

  for( int iHist=0; iHist<numHists-2; iHist++ ){
    //continue;
    h_csv_ratio_final[iHist] = new TH1D( Form("h_%s_final",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_JESUp[iHist] = new TH1D( Form("h_%s_final_JESUp",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_JESDown[iHist] = new TH1D( Form("h_%s_final_JESDown",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_LFUp[iHist] = new TH1D( Form("h_%s_final_LFUp",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_LFDown[iHist] = new TH1D( Form("h_%s_final_LFDown",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_Stats1Up[iHist] = new TH1D( Form("h_%s_final_Stats1Up",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_Stats1Down[iHist] = new TH1D( Form("h_%s_final_Stats1Down",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_Stats2Up[iHist] = new TH1D( Form("h_%s_final_Stats2Up",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );
    h_csv_ratio_final_Stats2Down[iHist] = new TH1D( Form("h_%s_final_Stats2Down",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );


    if( iHist<NumHists_normal ){
      for( int iBin=0; iBin<ncsvbins; iBin++ ){
	h_csv_ratio_final_StatsUp[iHist][iBin] = new TH1D( Form("h_%s_final_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1), ";CSV", NumFinalBins, -0.04, 1.01 );
	h_csv_ratio_final_StatsDown[iHist][iBin] = new TH1D( Form("h_%s_final_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1), ";CSV", NumFinalBins, -0.04, 1.01 );
      }
    }

    c_csv_ratio_final[iHist] = new TH1D( Form("c_%s_final",hist_name[iHist].Data()), ";CSV", NumFinalBins, -0.04, 1.01 );

    std::vector<double> x; x.clear();
    std::vector<double> y; y.clear();
    std::vector<double> y_JESUp; y_JESUp.clear();
    std::vector<double> y_JESDown; y_JESDown.clear();
    std::vector<double> y_LFUp; y_LFUp.clear();
    std::vector<double> y_LFDown; y_LFDown.clear();
    std::vector<double> y_Stats1Up; y_Stats1Up.clear();
    std::vector<double> y_Stats1Down; y_Stats1Down.clear();
    std::vector<double> y_Stats2Up; y_Stats2Up.clear();
    std::vector<double> y_Stats2Down; y_Stats2Down.clear();

    std::vector<std::vector<double> > y_temp_StatsUp;
    std::vector<std::vector<double> > y_temp_StatsDown;


    for( int iBin=0; iBin<nBins; iBin++ ){
      double center = h_csv_ratio[iHist]->GetBinCenter(iBin+1);
      if( center<0 ) continue;
      double sf = h_csv_ratio[iHist]->GetBinContent(iBin+1);
      double sf_JESUp = h_csv_ratio_JESUp[iHist]->GetBinContent(iBin+1);
      double sf_JESDown = h_csv_ratio_JESDown[iHist]->GetBinContent(iBin+1);
      double sf_LFUp = h_csv_ratio_LFUp[iHist]->GetBinContent(iBin+1);
      double sf_LFDown = h_csv_ratio_LFDown[iHist]->GetBinContent(iBin+1);
      double sf_Stats1Up = h_csv_ratio_Stats1Up[iHist]->GetBinContent(iBin+1);
      double sf_Stats1Down = h_csv_ratio_Stats1Down[iHist]->GetBinContent(iBin+1);
      double sf_Stats2Up = h_csv_ratio_Stats2Up[iHist]->GetBinContent(iBin+1);
      double sf_Stats2Down = h_csv_ratio_Stats2Down[iHist]->GetBinContent(iBin+1);

      x.push_back(center);
      y.push_back(sf);
      y_JESUp.push_back(sf_JESUp);
      y_JESDown.push_back(sf_JESDown);
      y_LFUp.push_back(sf_LFUp);
      y_LFDown.push_back(sf_LFDown);
      y_Stats1Up.push_back(sf_Stats1Up);
      y_Stats1Down.push_back(sf_Stats1Down);
      y_Stats2Up.push_back(sf_Stats2Up);
      y_Stats2Down.push_back(sf_Stats2Down);

      std::vector<double> temp_y_StatsUp;
      std::vector<double> temp_y_StatsDown;

      if( iHist<NumHists_normal ){
	for( int jBin=0; jBin<nBins; jBin++ ){
	  double sf_StatsUp = h_csv_ratio_StatsUp[iHist][jBin]->GetBinContent(iBin+1);
	  double sf_StatsDown = h_csv_ratio_StatsDown[iHist][jBin]->GetBinContent(iBin+1);
	
	  temp_y_StatsUp.push_back(sf_StatsUp);
	  temp_y_StatsDown.push_back(sf_StatsDown);
	  //printf("\t\t\t StatsUp,\t iBin = %d,\t jBin = %d,\t sf_StatsUp = %.3f,\t sf_StatsDown = %.3f \n", iBin, jBin, sf_StatsUp, sf_StatsDown );
	}
	y_temp_StatsUp.push_back(temp_y_StatsUp);
	y_temp_StatsDown.push_back(temp_y_StatsDown);
      }
    }


    std::vector<std::vector<double> > y_StatsUp( y_temp_StatsUp[0].size(), std::vector<double>(y_temp_StatsUp.size()) );
    std::vector<std::vector<double> > y_StatsDown( y_temp_StatsDown[0].size(), std::vector<double>(y_temp_StatsDown.size()) );


    if( iHist<NumHists_normal ){
      for( unsigned int i = 0; i < y_temp_StatsUp[0].size(); i++){
	for( unsigned int j = 0; j < y_temp_StatsUp.size(); j++) {
	  y_StatsUp[i][j] = y_temp_StatsUp[j][i];
	  y_StatsDown[i][j] = y_temp_StatsDown[j][i];
	  //printf("\t\t StatsUp,\t i = %d,\t j = %d,\t y_StatsUp = %.3f \n", i, j, y_StatsUp[i][j]);
	}
      }
    }


    //ROOT::Math::Interpolation::Type iType = ROOT::Math::Interpolation::kLINEAR;
    //ROOT::Math::Interpolation::Type iType = ROOT::Math::Interpolation::kCSPLINE;
    //ROOT::Math::Interpolation::Type iType = ROOT::Math::Interpolation::kPOLYNOMIAL;
    ROOT::Math::Interpolation::Type iType = ROOT::Math::Interpolation::kAKIMA;

    ROOT::Math::Interpolator * m_i =  new ROOT::Math::Interpolator(x,y,iType);
    ROOT::Math::Interpolator * m_i_JESUp =  new ROOT::Math::Interpolator(x,y_JESUp,iType);
    ROOT::Math::Interpolator * m_i_JESDown =  new ROOT::Math::Interpolator(x,y_JESDown,iType);
    ROOT::Math::Interpolator * m_i_LFUp =  new ROOT::Math::Interpolator(x,y_LFUp,iType);
    ROOT::Math::Interpolator * m_i_LFDown =  new ROOT::Math::Interpolator(x,y_LFDown,iType);
    ROOT::Math::Interpolator * m_i_Stats1Up =  new ROOT::Math::Interpolator(x,y_Stats1Up,iType);
    ROOT::Math::Interpolator * m_i_Stats1Down =  new ROOT::Math::Interpolator(x,y_Stats1Down,iType);
    ROOT::Math::Interpolator * m_i_Stats2Up =  new ROOT::Math::Interpolator(x,y_Stats2Up,iType);
    ROOT::Math::Interpolator * m_i_Stats2Down =  new ROOT::Math::Interpolator(x,y_Stats2Down,iType);

    ROOT::Math::Interpolator * m_i_StatsUp[nBins];
    ROOT::Math::Interpolator * m_i_StatsDown[nBins];


    if( iHist<NumHists_normal ){
      for( int iBin=0; iBin<nBins; iBin++ ){
	//std::cout << "\t x.size() = " << x.size() << ",\t y_StatsUp[" << iBin << "].size() = " << y_StatsUp[iBin].size() << std::endl;
	m_i_StatsUp[iBin] = new ROOT::Math::Interpolator(x,y_StatsUp[iBin],iType);
	//std::cout << "\t x.size() = " << x.size() << ",\t y_StatsDown[" << iBin << "].size() = " << y_StatsUp[iBin].size() << std::endl;
	m_i_StatsDown[iBin] = new ROOT::Math::Interpolator(x,y_StatsDown[iBin],iType);
      }
    }


    h_csv_ratio_final[iHist]->SetLineWidth(2);
    h_csv_ratio_final_JESUp[iHist]->SetLineWidth(2);
    h_csv_ratio_final_JESDown[iHist]->SetLineWidth(2);
    h_csv_ratio_final_LFUp[iHist]->SetLineWidth(2);
    h_csv_ratio_final_LFDown[iHist]->SetLineWidth(2);
    h_csv_ratio_final_Stats1Up[iHist]->SetLineWidth(2);
    h_csv_ratio_final_Stats1Down[iHist]->SetLineWidth(2);
    h_csv_ratio_final_Stats2Up[iHist]->SetLineWidth(2);
    h_csv_ratio_final_Stats2Down[iHist]->SetLineWidth(2);


    h_csv_ratio_final[iHist]->SetLineColor(kBlack);
    h_csv_ratio_final_JESUp[iHist]->SetLineColor(kBlue);
    h_csv_ratio_final_JESDown[iHist]->SetLineColor(kBlue);
    h_csv_ratio_final_LFUp[iHist]->SetLineColor(kGreen+3);
    h_csv_ratio_final_LFDown[iHist]->SetLineColor(kGreen+3);
    h_csv_ratio_final_Stats1Up[iHist]->SetLineColor(kMagenta-2);
    h_csv_ratio_final_Stats1Down[iHist]->SetLineColor(kMagenta-2);
    h_csv_ratio_final_Stats2Up[iHist]->SetLineColor(kRed-2);
    h_csv_ratio_final_Stats2Down[iHist]->SetLineColor(kRed-2);

    if( iHist<NumHists_normal ){
      for( int iBin=0; iBin<nBins; iBin++ ){
	h_csv_ratio_final_StatsUp[iHist][iBin]->SetLineWidth(2);
	h_csv_ratio_final_StatsDown[iHist][iBin]->SetLineWidth(2);

	h_csv_ratio_final_StatsUp[iHist][iBin]->SetLineColor(kRed);
	h_csv_ratio_final_StatsDown[iHist][iBin]->SetLineColor(kBlue);
      }
    }


    double firstPoint = x[0];
    double lastPoint  = x[int(x.size())-1];
    for( int iBin=0; iBin<n; iBin++ ){
      double center = h_csv_ratio_final[iHist]->GetBinCenter(iBin+1);
      if( center<0 ){
	h_csv_ratio_final[iHist]->SetBinContent(iBin+1,h_csv_ratio[iHist]->GetBinContent(1));
	h_csv_ratio_final_JESUp[iHist]->SetBinContent(iBin+1,h_csv_ratio_JESUp[iHist]->GetBinContent(1));
	h_csv_ratio_final_JESDown[iHist]->SetBinContent(iBin+1,h_csv_ratio_JESDown[iHist]->GetBinContent(1));
	h_csv_ratio_final_LFUp[iHist]->SetBinContent(iBin+1,h_csv_ratio_LFUp[iHist]->GetBinContent(1));
	h_csv_ratio_final_LFDown[iHist]->SetBinContent(iBin+1,h_csv_ratio_LFDown[iHist]->GetBinContent(1));
	h_csv_ratio_final_Stats1Up[iHist]->SetBinContent(iBin+1,h_csv_ratio[iHist]->GetBinContent(1) + h_csv_ratio[iHist]->GetBinError(1));
	h_csv_ratio_final_Stats1Down[iHist]->SetBinContent(iBin+1,h_csv_ratio[iHist]->GetBinContent(1) - h_csv_ratio[iHist]->GetBinError(1));
	h_csv_ratio_final_Stats2Up[iHist]->SetBinContent(iBin+1,h_csv_ratio[iHist]->GetBinContent(1) + h_csv_ratio[iHist]->GetBinError(1));
	h_csv_ratio_final_Stats2Down[iHist]->SetBinContent(iBin+1,h_csv_ratio[iHist]->GetBinContent(1) - h_csv_ratio[iHist]->GetBinError(1));

	if( iHist<NumHists_normal ){
	  for( int jBin=0; jBin<nBins; jBin++ ){
	    h_csv_ratio_final_StatsUp[iHist][jBin]->SetBinContent(iBin+1,h_csv_ratio_StatsUp[iHist][jBin]->GetBinContent(1));
	    h_csv_ratio_final_StatsDown[iHist][jBin]->SetBinContent(iBin+1,h_csv_ratio_StatsDown[iHist][jBin]->GetBinContent(1));
	  }
	}
      }
      else if( center<firstPoint ){
      	h_csv_ratio_final[iHist]->SetBinContent(iBin+1,m_i->Eval( firstPoint ));
      	h_csv_ratio_final_JESUp[iHist]->SetBinContent(iBin+1,m_i_JESUp->Eval( firstPoint ));
      	h_csv_ratio_final_JESDown[iHist]->SetBinContent(iBin+1,m_i_JESDown->Eval( firstPoint ));
      	h_csv_ratio_final_LFUp[iHist]->SetBinContent(iBin+1,m_i_LFUp->Eval( firstPoint ));
      	h_csv_ratio_final_LFDown[iHist]->SetBinContent(iBin+1,m_i_LFDown->Eval( firstPoint ));
      	h_csv_ratio_final_Stats1Up[iHist]->SetBinContent(iBin+1,m_i_Stats1Up->Eval( firstPoint ));
      	h_csv_ratio_final_Stats1Down[iHist]->SetBinContent(iBin+1,m_i_Stats1Down->Eval( firstPoint ));
      	h_csv_ratio_final_Stats2Up[iHist]->SetBinContent(iBin+1,m_i_Stats2Up->Eval( firstPoint ));
      	h_csv_ratio_final_Stats2Down[iHist]->SetBinContent(iBin+1,m_i_Stats2Down->Eval( firstPoint ));

	if( iHist<NumHists_normal ){
	  for( int jBin=0; jBin<nBins; jBin++ ){
	    h_csv_ratio_final_StatsUp[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsUp[jBin]->Eval( firstPoint ));
	    h_csv_ratio_final_StatsDown[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsDown[jBin]->Eval( firstPoint ));
	  }
	}
      }
      else if( center>lastPoint ){
      	h_csv_ratio_final[iHist]->SetBinContent(iBin+1,m_i->Eval( lastPoint ));
      	h_csv_ratio_final_JESUp[iHist]->SetBinContent(iBin+1,m_i_JESUp->Eval( lastPoint ));
      	h_csv_ratio_final_JESDown[iHist]->SetBinContent(iBin+1,m_i_JESDown->Eval( lastPoint ));
      	h_csv_ratio_final_LFUp[iHist]->SetBinContent(iBin+1,m_i_LFUp->Eval( lastPoint ));
      	h_csv_ratio_final_LFDown[iHist]->SetBinContent(iBin+1,m_i_LFDown->Eval( lastPoint ));
      	h_csv_ratio_final_Stats1Up[iHist]->SetBinContent(iBin+1,m_i_Stats1Up->Eval( lastPoint ));
      	h_csv_ratio_final_Stats1Down[iHist]->SetBinContent(iBin+1,m_i_Stats1Down->Eval( lastPoint ));
      	h_csv_ratio_final_Stats2Up[iHist]->SetBinContent(iBin+1,m_i_Stats2Up->Eval( lastPoint ));
      	h_csv_ratio_final_Stats2Down[iHist]->SetBinContent(iBin+1,m_i_Stats2Down->Eval( lastPoint ));

	if( iHist<NumHists_normal ){
	  for( int jBin=0; jBin<nBins; jBin++ ){
	    h_csv_ratio_final_StatsUp[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsUp[jBin]->Eval( lastPoint ));
	    h_csv_ratio_final_StatsDown[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsDown[jBin]->Eval( lastPoint ));
	  }
	}
      }
      else {
	h_csv_ratio_final[iHist]->SetBinContent(iBin+1,m_i->Eval( center ));
	h_csv_ratio_final_JESUp[iHist]->SetBinContent(iBin+1,m_i_JESUp->Eval( center ));
	h_csv_ratio_final_JESDown[iHist]->SetBinContent(iBin+1,m_i_JESDown->Eval( center ));
	h_csv_ratio_final_LFUp[iHist]->SetBinContent(iBin+1,m_i_LFUp->Eval( center ));
	h_csv_ratio_final_LFDown[iHist]->SetBinContent(iBin+1,m_i_LFDown->Eval( center ));
	h_csv_ratio_final_Stats1Up[iHist]->SetBinContent(iBin+1,m_i_Stats1Up->Eval( center ));
	h_csv_ratio_final_Stats1Down[iHist]->SetBinContent(iBin+1,m_i_Stats1Down->Eval( center ));
	h_csv_ratio_final_Stats2Up[iHist]->SetBinContent(iBin+1,m_i_Stats2Up->Eval( center ));
	h_csv_ratio_final_Stats2Down[iHist]->SetBinContent(iBin+1,m_i_Stats2Down->Eval( center ));

	if( iHist<NumHists_normal ){
	  for( int jBin=0; jBin<nBins; jBin++ ){
	    h_csv_ratio_final_StatsUp[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsUp[jBin]->Eval( center ));
	    h_csv_ratio_final_StatsDown[iHist][jBin]->SetBinContent(iBin+1,m_i_StatsDown[jBin]->Eval( center ));
	  }
	}
      }
    }


    h_csv_ratio[iHist]->SetMarkerStyle(20);

    h_csv_ratio[iHist]->SetTitle( Form("HF %s", hist_name[iHist].Data() ) );
    if( iHist==numHists-1 ) h_csv_ratio[iHist]->SetTitle( Form("LF %s", hist_name[iHist].Data() ) );
    h_csv_ratio[iHist]->GetYaxis()->SetTitle("Data/MC SF");
    h_csv_ratio[iHist]->GetXaxis()->SetTitle("CSV");

    h_csv_ratio[iHist]->SetStats(0);

    double maxY = 1.3 * h_csv_ratio[iHist]->GetMaximum();
    maxY = std::min( maxY, 5. );

    h_csv_ratio[iHist]->SetMaximum(maxY);
    h_csv_ratio[iHist]->SetMinimum(0.2);

    TLegend *legend = new TLegend(0.2,0.75,0.87,0.89);


    h_csv_ratio_LFUp[iHist]->SetLineColor(kGreen+1);
    h_csv_ratio_LFDown[iHist]->SetLineColor(kGreen+1);

    //h_csv_ratio_JESUp[iHist]->SetLineColor(kRed);
    h_csv_ratio_JESUp[iHist]->SetLineColor(kBlue);
    h_csv_ratio_JESDown[iHist]->SetLineColor(kBlue);



    legend->SetFillColor(kWhite);
    legend->SetLineColor(kWhite);
    legend->SetShadowColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);

    legend->SetNColumns(3);

    legend->AddEntry(h_csv_ratio[iHist],"HF SF","p");
    legend->AddEntry(h_csv_ratio_final_JESUp[iHist],"JES Unc","l");
    legend->AddEntry(h_csv_ratio_final[iHist],"HF Fit","l");
    legend->AddEntry(h_csv_ratio_final_LFUp[iHist],"LF Unc","l");
    legend->AddEntry(h_csv_ratio_final_Stats1Up[iHist],"Stats1 Unc","l");
    legend->AddEntry(h_csv_ratio_final_Stats2Up[iHist],"Stats2 Unc","l");
    //legend->AddEntry(gr,"BTV HF SF","pf");


    TLegend *legend_reduced = new TLegend(0.2,0.82,0.87,0.89);

    legend_reduced->SetFillColor(kWhite);
    legend_reduced->SetLineColor(kWhite);
    legend_reduced->SetShadowColor(kWhite);
    legend_reduced->SetTextFont(42);
    legend_reduced->SetTextSize(0.035);

    legend_reduced->SetNColumns(3);

    legend_reduced->AddEntry(h_csv_ratio[iHist],"HF SF","p");
    legend_reduced->AddEntry(h_csv_ratio_final[iHist],"HF Fit","l");
    legend_reduced->AddEntry(gr,"BTV HF SF","pf");


    TString img;



    h_csv_ratio[iHist]->Draw("pe1");
    //gr->Draw("pe2same");
    //h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend_reduced->Draw();
 
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_only.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_only.pdf";
    if( makePlots ) c1->Print(img);



    h_csv_ratio[iHist]->Draw("pe1");
    h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio_final_JESUp[iHist]->Draw("histsame");
    h_csv_ratio_final_JESDown[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend->Draw();
 
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_JES.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_JES.pdf";
    if( makePlots ) c1->Print(img);


    h_csv_ratio[iHist]->Draw("pe1");
    h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio_final_LFUp[iHist]->Draw("histsame");
    h_csv_ratio_final_LFDown[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend->Draw();
 
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_LF.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_LF.pdf";
    if( makePlots ) c1->Print(img);


    h_csv_ratio[iHist]->Draw("pe1");
    h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats1Up[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats1Down[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend->Draw();
 
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_Stats1.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_Stats1.pdf";
    if( makePlots ) c1->Print(img);



    h_csv_ratio[iHist]->Draw("pe1");
    h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats2Up[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats2Down[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend->Draw();
 
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_Stats2.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_Stats2.pdf";
    if( makePlots ) c1->Print(img);



    h_csv_ratio_StatsUp[0][0]->SetLineColor(kRed);
    h_csv_ratio_StatsDown[0][0]->SetLineColor(kBlue);

    TLegend *legend_stats = new TLegend(0.2,0.82,0.87,0.89);

    legend_stats->SetFillColor(kWhite);
    legend_stats->SetLineColor(kWhite);
    legend_stats->SetShadowColor(kWhite);
    legend_stats->SetTextFont(42);
    legend_stats->SetTextSize(0.035);

    legend_stats->SetNColumns(3);

    legend_stats->AddEntry(h_csv_ratio[iHist],"HF SF","p");
    legend_stats->AddEntry(h_csv_ratio_StatsUp[0][0],"StatsUp","l");
    legend_stats->AddEntry(h_csv_ratio_StatsDown[0][0],"StatsDown","l");

    if( iHist<NumHists_normal ){
      for( int iBin=0; iBin<ncsvbins; iBin++ ){

	h_csv_ratio_StatsUp[iHist][iBin]->SetLineColor(kRed);
	h_csv_ratio_StatsDown[iHist][iBin]->SetLineColor(kBlue);

	h_csv_ratio[iHist]->Draw("pe1");
	h_csv_ratio_final[iHist]->Draw("histsame");
	h_csv_ratio_final_StatsUp[iHist][iBin]->Draw("histsame");
	h_csv_ratio_final_StatsDown[iHist][iBin]->Draw("histsame");
	// h_csv_ratio_StatsUp[iHist][iBin]->Draw("histsame");
	// h_csv_ratio_StatsDown[iHist][iBin]->Draw("histsame");
	h_csv_ratio[iHist]->Draw("pe1same");

	legend_stats->Draw();

	if( makePlots ){
	  c1->Print(Form("%shfSF_%s_fit_CSVbin%dStats.png",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));
	  c1->Print(Form("%shfSF_%s_fit_CSVbin%dStats.pdf",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));
	}
      }
    }


    h_csv_ratio[iHist]->Draw("pe1");
    //gr->Draw("pe2same");
    h_csv_ratio_final[iHist]->Draw("histsame");
    h_csv_ratio_final_JESUp[iHist]->Draw("histsame");
    h_csv_ratio_final_JESDown[iHist]->Draw("histsame");
    h_csv_ratio_final_LFUp[iHist]->Draw("histsame");
    h_csv_ratio_final_LFDown[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats1Up[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats1Down[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats2Up[iHist]->Draw("histsame");
    h_csv_ratio_final_Stats2Down[iHist]->Draw("histsame");
    h_csv_ratio[iHist]->Draw("pe1same");

    legend->Draw();

    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_All.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "hfSF_" + hist_name[iHist] + "_fit_All.pdf";
    if( makePlots ) c1->Print(img);


    h_csv_ratio_final[iHist]->Write(Form("%s_final",hist_name[iHist].Data()));
    h_csv_ratio_final_JESUp[iHist]->Write(Form("%s_final_JESUp",hist_name[iHist].Data()));
    h_csv_ratio_final_JESDown[iHist]->Write(Form("%s_final_JESDown",hist_name[iHist].Data()));
    h_csv_ratio_final_LFUp[iHist]->Write(Form("%s_final_LFUp",hist_name[iHist].Data()));
    h_csv_ratio_final_LFDown[iHist]->Write(Form("%s_final_LFDown",hist_name[iHist].Data()));
    h_csv_ratio_final_Stats1Up[iHist]->Write(Form("%s_final_Stats1Up",hist_name[iHist].Data()));
    h_csv_ratio_final_Stats1Down[iHist]->Write(Form("%s_final_Stats1Down",hist_name[iHist].Data()));
    h_csv_ratio_final_Stats2Up[iHist]->Write(Form("%s_final_Stats2Up",hist_name[iHist].Data()));
    h_csv_ratio_final_Stats2Down[iHist]->Write(Form("%s_final_Stats2Down",hist_name[iHist].Data()));

    if( iHist<NumHists_normal ){
      for( int iBin=0; iBin<nBins; iBin++ ){
	h_csv_ratio_final_StatsUp[iHist][iBin]->Write(Form("%s_final_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1));
	h_csv_ratio_final_StatsDown[iHist][iBin]->Write(Form("%s_final_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1));
      }
    }

    for( int iBin=0; iBin<c_csv_ratio_final[iHist]->GetNbinsX(); iBin++ ) c_csv_ratio_final[iHist]->SetBinContent(iBin+1,1.0);

    TH1D* h_ratio_JESUp = (TH1D*)h_csv_ratio_final_JESUp[iHist]->Clone("h_ratio_JESUp_"+hist_name[iHist]);
    TH1D* h_ratio_JESDown = (TH1D*)h_csv_ratio_final_JESDown[iHist]->Clone("h_ratio_JESDown_"+hist_name[iHist]);
    TH1D* h_ratio_LFUp = (TH1D*)h_csv_ratio_final_LFUp[iHist]->Clone("h_ratio_LFUp_"+hist_name[iHist]);
    TH1D* h_ratio_LFDown = (TH1D*)h_csv_ratio_final_LFDown[iHist]->Clone("h_ratio_LFDown_"+hist_name[iHist]);
    TH1D* h_ratio_Stats1Up = (TH1D*)h_csv_ratio_final_Stats1Up[iHist]->Clone("h_ratio_Stats1Up_"+hist_name[iHist]);
    TH1D* h_ratio_Stats1Down = (TH1D*)h_csv_ratio_final_Stats1Down[iHist]->Clone("h_ratio_Stats1Down_"+hist_name[iHist]);
    TH1D* h_ratio_Stats2Up = (TH1D*)h_csv_ratio_final_Stats2Up[iHist]->Clone("h_ratio_Stats2Up_"+hist_name[iHist]);
    TH1D* h_ratio_Stats2Down = (TH1D*)h_csv_ratio_final_Stats2Down[iHist]->Clone("h_ratio_Stats2Down_"+hist_name[iHist]);

    TH1D* h_total_cUp   = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_total_cUp_"+hist_name[iHist]);
    TH1D* h_total_cDown = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_total_cDown_"+hist_name[iHist]);

    h_ratio_JESUp->Divide(h_csv_ratio_final[iHist]);
    h_ratio_JESDown->Divide(h_csv_ratio_final[iHist]);
    h_ratio_LFUp->Divide(h_csv_ratio_final[iHist]);
    h_ratio_LFDown->Divide(h_csv_ratio_final[iHist]);
    h_ratio_Stats1Up->Divide(h_csv_ratio_final[iHist]);
    h_ratio_Stats1Down->Divide(h_csv_ratio_final[iHist]);
    h_ratio_Stats2Up->Divide(h_csv_ratio_final[iHist]);
    h_ratio_Stats2Down->Divide(h_csv_ratio_final[iHist]);

    h_ratio_JESUp->SetLineWidth(2);
    h_ratio_JESDown->SetLineWidth(2);
    h_ratio_LFUp->SetLineWidth(2);
    h_ratio_LFDown->SetLineWidth(2);
    h_ratio_Stats1Up->SetLineWidth(2);
    h_ratio_Stats1Down->SetLineWidth(2);
    h_ratio_Stats2Up->SetLineWidth(2);
    h_ratio_Stats2Down->SetLineWidth(2);

    h_ratio_JESUp->SetLineColor(kBlue);
    h_ratio_JESDown->SetLineColor(kBlue);
    h_ratio_LFUp->SetLineColor(kGreen+1);
    h_ratio_LFDown->SetLineColor(kGreen+1);
    h_ratio_Stats1Up->SetLineColor(kMagenta-2);
    h_ratio_Stats1Down->SetLineColor(kMagenta-2);
    h_ratio_Stats2Up->SetLineColor(kRed-2);
    h_ratio_Stats2Down->SetLineColor(kRed-2);


    for( int iBin=0; iBin<c_csv_ratio_final[iHist]->GetNbinsX(); iBin++ ){
      double nominal = 1.;
      double diff_JESUp = h_ratio_JESUp->GetBinContent(iBin+1)     - nominal;
      double diff_JESDown = h_ratio_JESDown->GetBinContent(iBin+1) - nominal;
      double diff_LFUp = h_ratio_LFUp->GetBinContent(iBin+1)     - nominal;
      double diff_LFDown = h_ratio_LFDown->GetBinContent(iBin+1) - nominal;
      double diff_Stats1Up = h_ratio_Stats1Up->GetBinContent(iBin+1)     - nominal;
      double diff_Stats1Down = h_ratio_Stats1Down->GetBinContent(iBin+1) - nominal;
      double diff_Stats2Up = h_ratio_Stats2Up->GetBinContent(iBin+1)     - nominal;
      double diff_Stats2Down = h_ratio_Stats2Down->GetBinContent(iBin+1) - nominal;

      double sum_up2 = 0., sum_down2 = 0.;
      if( diff_JESUp>0 ) sum_up2   += diff_JESUp*diff_JESUp;
      else               sum_down2 += diff_JESUp*diff_JESUp;
      if( diff_JESDown>0 ) sum_up2   += diff_JESDown*diff_JESDown;
      else                 sum_down2 += diff_JESDown*diff_JESDown;

      if( diff_LFUp>0 ) sum_up2   += diff_LFUp*diff_LFUp;
      else              sum_down2 += diff_LFUp*diff_LFUp;
      if( diff_LFDown>0 ) sum_up2   += diff_LFDown*diff_LFDown;
      else                sum_down2 += diff_LFDown*diff_LFDown;

      if( diff_Stats1Up>0 ) sum_up2   += diff_Stats1Up*diff_Stats1Up;
      else                  sum_down2 += diff_Stats1Up*diff_Stats1Up;
      if( diff_Stats1Down>0 ) sum_up2   += diff_Stats1Down*diff_Stats1Down;
      else                    sum_down2 += diff_Stats1Down*diff_Stats1Down;

      if( diff_Stats2Up>0 ) sum_up2   += diff_Stats2Up*diff_Stats2Up;
      else                  sum_down2 += diff_Stats2Up*diff_Stats2Up;
      if( diff_Stats2Down>0 ) sum_up2   += diff_Stats2Down*diff_Stats2Down;
      else                    sum_down2 += diff_Stats2Down*diff_Stats2Down;


      double total_cUp_bin = nominal+sqrt(sum_up2);
      double total_cDown_bin = nominal-sqrt(sum_down2);

      h_total_cUp->SetBinContent(iBin+1,total_cUp_bin);
      h_total_cDown->SetBinContent(iBin+1,total_cDown_bin);

    }



    for( int iBin=0; iBin<h_csv_ratio[iHist]->GetNbinsX(); iBin++ ){
      double sf = h_csv_ratio[iHist]->GetBinContent(iBin+1);
      double diff_JESUp = h_csv_ratio_JESUp[iHist]->GetBinContent(iBin+1) - sf;
      double diff_JESDown = h_csv_ratio_JESDown[iHist]->GetBinContent(iBin+1) - sf;
      double diff_LFUp = h_csv_ratio_LFUp[iHist]->GetBinContent(iBin+1) - sf;
      double diff_LFDown = h_csv_ratio_LFDown[iHist]->GetBinContent(iBin+1) - sf;
      double diff_Stats1Up = h_csv_ratio_Stats1Up[iHist]->GetBinContent(iBin+1) - sf;
      double diff_Stats1Down = h_csv_ratio_Stats1Down[iHist]->GetBinContent(iBin+1) - sf;
      double diff_Stats2Up = h_csv_ratio_Stats2Up[iHist]->GetBinContent(iBin+1) - sf;
      double diff_Stats2Down = h_csv_ratio_Stats2Down[iHist]->GetBinContent(iBin+1) - sf;

      double diff_JES    = 0.5 * ( fabs(diff_JESUp) + fabs(diff_JESDown) );
      double diff_LF     = 0.5 * ( fabs(diff_LFUp) + fabs(diff_LFDown) );
      double diff_Stats1 = 0.5 * ( fabs(diff_Stats1Up) + fabs(diff_Stats1Down) );
      double diff_Stats2 = 0.5 * ( fabs(diff_Stats2Up) + fabs(diff_Stats2Down) );

      double sum_up2 = 0., sum_down2 = 0.;
      if( diff_JESUp>0 ) sum_up2   += diff_JESUp*diff_JESUp;
      else               sum_down2 += diff_JESUp*diff_JESUp;
      if( diff_JESDown>0 ) sum_up2   += diff_JESDown*diff_JESDown;
      else                 sum_down2 += diff_JESDown*diff_JESDown;

      if( diff_LFUp>0 ) sum_up2   += diff_LFUp*diff_LFUp;
      else              sum_down2 += diff_LFUp*diff_LFUp;
      if( diff_LFDown>0 ) sum_up2   += diff_LFDown*diff_LFDown;
      else                sum_down2 += diff_LFDown*diff_LFDown;

      if( diff_Stats1Up>0 ) sum_up2   += diff_Stats1Up*diff_Stats1Up;
      else                  sum_down2 += diff_Stats1Up*diff_Stats1Up;
      if( diff_Stats1Down>0 ) sum_up2   += diff_Stats1Down*diff_Stats1Down;
      else                    sum_down2 += diff_Stats1Down*diff_Stats1Down;

      if( diff_Stats2Up>0 ) sum_up2   += diff_Stats2Up*diff_Stats2Up;
      else                  sum_down2 += diff_Stats2Up*diff_Stats2Up;
      if( diff_Stats2Down>0 ) sum_up2   += diff_Stats2Down*diff_Stats2Down;
      else                    sum_down2 += diff_Stats2Down*diff_Stats2Down;

      double sum_err = 0.5 * ( sqrt(sum_up2) + sqrt(sum_down2) );
      if( verbose ){
	printf(" iHist = %d,\t iBin = %d,\t SF = %.3f,\t TotalErr = %.3f,\t JESErr = %.3f,\t LFErr = %.3f,\t Stats1Err = %.3f,\t Stats2Err = %.3f,\t binErr = %.3f \n", 
	       iHist, iBin, sf, sum_err, diff_JES, diff_LF, diff_Stats1, diff_Stats2, h_csv_ratio[iHist]->GetBinError(iBin+1) );
      }
    }


    h_total_cUp->SetLineWidth(3);
    h_total_cDown->SetLineWidth(3);
    h_total_cUp->SetLineColor(kGreen+3);
    h_total_cDown->SetLineColor(kCyan+3);

    TH1D* h_cErr1Up   = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_cErr1Up_"+hist_name[iHist]);
    TH1D* h_cErr1Down = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_cErr1Down_"+hist_name[iHist]);
    TH1D* h_cErr2Up   = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_cErr2Up_"+hist_name[iHist]);
    TH1D* h_cErr2Down = (TH1D*)c_csv_ratio_final[iHist]->Clone("h_cErr2Down_"+hist_name[iHist]);

    for( int iBin=0; iBin<c_csv_ratio_final[iHist]->GetNbinsX(); iBin++ ){
      double center = h_total_cUp->GetBinCenter(iBin+1);
      double content = 1.0;

      double bSF = h_csv_ratio_final[iHist]->GetBinContent(iBin+1);
      double cUp = content + 2*fabs(h_total_cUp->GetBinContent(iBin+1) - content);
      double cDown = content - 2*fabs(h_total_cDown->GetBinContent(iBin+1) - content);
      if( (bSF - cUp)>0 )   cUp = bSF;
      if( (cDown - bSF)>0 ) cDown = bSF;
      cUp -= content;
      cDown -= content;

      double deltaUp = fabs( cUp );
      double deltaDown = fabs( cDown );

      // double delta = 2*TMath::Max( fabs(cUp), fabs(h_total_cDown->GetBinContent(iBin+1) - content) );
      double delta = TMath::Max( deltaUp, deltaDown );

      double cErr1Up   = ( center<0. ) ? TMath::Max( bSF, content + cUp   ) : content + delta * ( 1 - 2*center );
      double cErr1Down = ( center<0. ) ? TMath::Min( bSF, content + cDown ) : content - delta * ( 1 - 2*center );

      double cErr2Up   = ( center<0. ) ? content : content + delta * ( 1 - 6*center*(1-center) );
      double cErr2Down = ( center<0. ) ? content : content - delta * ( 1 - 6*center*(1-center) );

      //if( center>0 && center<0.03 ) printf("\t\t bin = %d,\t cErr1Up = %.1f,\t cErr1Down = %.1f,\t cErr2Up = %.1f,\t  cErr2Down = %.1f \n", iBin+1, cErr1Up, cErr1Down, cErr2Up, cErr2Down );

      h_cErr1Up->SetBinContent(iBin+1,cErr1Up);
      h_cErr1Down->SetBinContent(iBin+1,cErr1Down);

      h_cErr2Up->SetBinContent(iBin+1,cErr2Up);
      h_cErr2Down->SetBinContent(iBin+1,cErr2Down);
    }


    TLegend *legend3 = new TLegend(0.2,0.8,0.84,0.89);

    legend3->SetFillColor(kWhite);
    legend3->SetLineColor(kWhite);
    legend3->SetShadowColor(kWhite);
    legend3->SetTextFont(42);
    legend3->SetTextSize(0.035);

    legend3->SetNColumns(3);

    legend3->AddEntry(c_csv_ratio_final[iHist],"c SF","l");
    legend3->AddEntry(h_ratio_JESUp,"JES","l");
    legend3->AddEntry(h_ratio_LFUp,"LF","l");
    legend3->AddEntry(h_ratio_Stats1Up,"Stats1","l");
    legend3->AddEntry(h_ratio_Stats2Up,"Stats2","l");
    legend3->AddEntry(h_total_cUp,"Total Up","l");
    legend3->AddEntry(h_total_cDown,"Total Down","l");

    h_ratio_JESUp->SetTitle(";CSV;Uncertainty/Nominal");
    h_ratio_JESUp->SetStats(0);
    h_ratio_JESUp->GetYaxis()->SetRangeUser(0.5,1.5);
    h_ratio_JESUp->GetXaxis()->SetRangeUser(-0.041,0.9489);

    h_ratio_JESUp->Draw("hist");
    h_ratio_JESDown->Draw("histsame");
    h_ratio_LFUp->Draw("histsame");
    h_ratio_LFDown->Draw("histsame");
    h_ratio_Stats1Up->Draw("histsame");
    h_ratio_Stats1Down->Draw("histsame");
    h_ratio_Stats2Up->Draw("histsame");
    h_ratio_Stats2Down->Draw("histsame");
    c_csv_ratio_final[iHist]->Draw("histsame");
    h_total_cUp->Draw("histsame");
    h_total_cDown->Draw("histsame");

    legend3->Draw();

    img = dirprefix + "cSF_" + hist_name[iHist] + "_fit_All.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "cSF_" + hist_name[iHist] + "_fit_All.pdf";
    if( makePlots ) c1->Print(img);



    TH1D* h_total_cUp_x2 = (TH1D*)h_total_cUp->Clone("h_total_cUp_x2_"+hist_name[iHist]);
    TH1D* h_total_cDown_x2 = (TH1D*)h_total_cDown->Clone("h_total_cDown_x2_"+hist_name[iHist]);

    for( int iBin=0; iBin<h_total_cUp_x2->GetNbinsX(); iBin++ ){
      double up_x2   = c_csv_ratio_final[iHist]->GetBinContent(iBin+1) + 2 * ( h_total_cUp->GetBinContent(iBin+1) - c_csv_ratio_final[iHist]->GetBinContent(iBin+1) );
      double down_x2 = c_csv_ratio_final[iHist]->GetBinContent(iBin+1) - 2 * ( c_csv_ratio_final[iHist]->GetBinContent(iBin+1) - h_total_cDown->GetBinContent(iBin+1) );

      double bSF = h_csv_ratio_final[iHist]->GetBinContent(iBin+1);
      double center = h_csv_ratio_final[iHist]->GetBinCenter(iBin+1);
      if( (center>=0.244 && center<=0.679) || true ){
	if( (bSF - up_x2)>0.   ) up_x2 = bSF;
	if( (down_x2 - bSF)>0. ) down_x2 = bSF;
      }

      h_total_cUp_x2->SetBinContent(iBin+1,up_x2);
      h_total_cDown_x2->SetBinContent(iBin+1,down_x2);
    }

    h_total_cUp_x2->SetLineWidth(3);
    h_total_cDown_x2->SetLineWidth(3);
    h_total_cUp_x2->SetLineColor(kGreen+3);
    h_total_cDown_x2->SetLineColor(kGreen+3);


    h_cErr1Up->SetLineWidth(2);
    h_cErr1Down->SetLineWidth(2);
    h_cErr2Up->SetLineWidth(2);
    h_cErr2Down->SetLineWidth(2);

    h_cErr1Up->SetLineColor(kRed);
    h_cErr1Down->SetLineColor(kRed);
    h_cErr2Up->SetLineColor(kBlue);
    h_cErr2Down->SetLineColor(kBlue);

    h_total_cUp_x2->SetTitle(";CSV;Scale Factor");
    h_total_cUp_x2->SetStats(0);
    h_total_cUp_x2->GetYaxis()->SetRangeUser(0.,2.0);
    h_total_cUp_x2->GetXaxis()->SetRangeUser(-0.041,1.0);

    h_total_cUp_x2->Draw("hist");
    h_total_cDown_x2->Draw("histsame");

    h_cErr1Up->Draw("histsame");
    h_cErr1Down->Draw("histsame");
    h_cErr2Up->Draw("histsame");
    h_cErr2Down->Draw("histsame");

    c_csv_ratio_final[iHist]->Draw("histsame");


    TLegend *legend4 = new TLegend(0.3,0.8,0.89,0.89);

    legend4->SetFillColor(kWhite);
    legend4->SetLineColor(kWhite);
    legend4->SetShadowColor(kWhite);
    legend4->SetTextFont(42);
    legend4->SetTextSize(0.035);

    legend4->SetNColumns(4);

    legend4->AddEntry(c_csv_ratio_final[iHist],"c SF","l");
    legend4->AddEntry(h_total_cUp_x2,"c Tot Err","l");
    legend4->AddEntry(h_cErr1Up,"cErr1","l");
    legend4->AddEntry(h_cErr2Up,"cErr2","l");

    legend4->Draw();

    img = dirprefix + "cSF_totalErr_" + hist_name[iHist] + "_fit_All.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "cSF_totalErr_" + hist_name[iHist] + "_fit_All.pdf";
    if( makePlots ) c1->Print(img);



    h_total_cUp_x2->SetStats(0);
    h_total_cUp_x2->GetYaxis()->SetRangeUser(0.,2.0);
    h_total_cUp_x2->GetYaxis()->SetTitle("Scale Factor");
    h_total_cUp_x2->Draw("hist");
    h_total_cDown_x2->Draw("histsame");

    
    h_csv_ratio[iHist]->SetLineWidth(2);
    h_csv_ratio[iHist]->SetLineColor(kBlack);
    h_csv_ratio[iHist]->Draw("pe1same");
    h_total_cUp_x2->Draw("histsame");
    h_total_cDown_x2->Draw("histsame");



    TLegend *legend5 = new TLegend(0.15,0.91,0.84,0.99);

    legend5->SetFillColor(kWhite);
    legend5->SetLineColor(kWhite);
    legend5->SetShadowColor(kWhite);
    legend5->SetTextFont(42);
    legend5->SetTextSize(0.035);

    legend5->SetNColumns(3);

    legend5->AddEntry(c_csv_ratio_final[iHist],"Charm SF","l");
    legend5->AddEntry(h_total_cUp,"Charm SF Unc","l");
    legend5->AddEntry(h_csv_ratio[iHist],"Bottom SF","pl");

    legend5->Draw();


    // h_cErr1Up->Draw("histsame");
    // h_cErr1Down->Draw("histsame");
    // h_cErr2Up->Draw("histsame");
    // h_cErr2Down->Draw("histsame");

    c_csv_ratio_final[iHist]->Draw("histsame");


    img = dirprefix + "cSF_totalErr_compare_" + hist_name[iHist] + "_fit_All.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "cSF_totalErr_compare_" + hist_name[iHist] + "_fit_All.pdf";
    if( makePlots ) c1->Print(img);


    TH1D* hc = (TH1D*)hc_ttbar->Clone("hc_"+hist_name[iHist]);
    TH1D* hc_SF = (TH1D*)hc_ttbar->Clone("hc_SF_"+hist_name[iHist]);
    TH1D* hc_SF_cUp = (TH1D*)hc_ttbar->Clone("hc_SF_cUp_"+hist_name[iHist]);
    TH1D* hc_SF_cDown = (TH1D*)hc_ttbar->Clone("hc_SF_cDown_"+hist_name[iHist]);
    TH1D* hc_SF_cErr1Up = (TH1D*)hc_ttbar->Clone("hc_SF_cErr1Up_"+hist_name[iHist]);
    TH1D* hc_SF_cErr1Down = (TH1D*)hc_ttbar->Clone("hc_SF_cErr1Down_"+hist_name[iHist]);
    TH1D* hc_SF_cErr2Up = (TH1D*)hc_ttbar->Clone("hc_SF_cErr2Up_"+hist_name[iHist]);
    TH1D* hc_SF_cErr2Down = (TH1D*)hc_ttbar->Clone("hc_SF_cErr2Down_"+hist_name[iHist]);

    for( int iBin=0; iBin<hc->GetNbinsX(); iBin++ ){
      double center = hc->GetBinCenter(iBin+1);
      int useBin = c_csv_ratio_final[iHist]->FindBin(center);

      double content = hc->GetBinContent(iBin+1);
      hc_SF->SetBinContent(iBin+1,content*c_csv_ratio_final[iHist]->GetBinContent(useBin));
      hc_SF_cUp->SetBinContent(iBin+1,content*h_total_cUp->GetBinContent(useBin));
      hc_SF_cDown->SetBinContent(iBin+1,content*h_total_cDown->GetBinContent(useBin));
      hc_SF_cErr1Up->SetBinContent(iBin+1,content*h_cErr1Up->GetBinContent(useBin));
      hc_SF_cErr1Down->SetBinContent(iBin+1,content*h_cErr1Down->GetBinContent(useBin));
      hc_SF_cErr2Up->SetBinContent(iBin+1,content*h_cErr2Up->GetBinContent(useBin));
      hc_SF_cErr2Down->SetBinContent(iBin+1,content*h_cErr2Down->GetBinContent(useBin));
    }

    for( int iBin=0; iBin<hc->GetNbinsX(); iBin++ ){
      double sum = hc->Integral(iBin+1,hc->GetNbinsX()+1);
      hc_SF->SetBinContent(iBin+1,hc_SF->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cUp->SetBinContent(iBin+1,hc_SF_cUp->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cDown->SetBinContent(iBin+1,hc_SF_cDown->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cErr1Up->SetBinContent(iBin+1,hc_SF_cErr1Up->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cErr1Down->SetBinContent(iBin+1,hc_SF_cErr1Down->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cErr2Up->SetBinContent(iBin+1,hc_SF_cErr2Up->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
      hc_SF_cErr2Down->SetBinContent(iBin+1,hc_SF_cErr2Down->Integral(iBin+1,hc->GetNbinsX()+1)/sum);
    }

    int binCSVL = hc_SF->FindBin(0.244);
    int binCSVM = hc_SF->FindBin(0.679);
    int binCSVT = hc_SF->FindBin(0.898);
    if( verbose ){
    printf("==> CSVM: %s\t SF = %.2f,\t SF_cUp = %.2f,\t SF_cDown = %.2f,\t SF_cErr1Up = %.2f,\t SF_cErr1Down = %.2f,\t SF_cErr2Up = %.2f,\t SF_cErr2Down = %.2f \n",
	   hist_name[iHist].Data(), hc_SF->GetBinContent(binCSVM), hc_SF_cUp->GetBinContent(binCSVM), hc_SF_cDown->GetBinContent(binCSVM),
	   hc_SF_cErr1Up->GetBinContent(binCSVM), hc_SF_cErr1Down->GetBinContent(binCSVM), hc_SF_cErr2Up->GetBinContent(binCSVM), hc_SF_cErr2Down->GetBinContent(binCSVM) );
    }

    //CSVL
    double cErr1Up_CSVL = hc_SF_cErr1Up->GetBinContent(binCSVL) - hc_SF->GetBinContent(binCSVL);
    double cErr1Down_CSVL = hc_SF_cErr1Down->GetBinContent(binCSVL) - hc_SF->GetBinContent(binCSVL);
    double cErr2Up_CSVL = hc_SF_cErr2Up->GetBinContent(binCSVL) - hc_SF->GetBinContent(binCSVL);
    double cErr2Down_CSVL = hc_SF_cErr2Down->GetBinContent(binCSVL) - hc_SF->GetBinContent(binCSVL);

    double cErrUp_CSVL_2 = 0, cErrDown_CSVL_2 = 0;
    if( cErr1Up_CSVL>0 ) cErrUp_CSVL_2 += cErr1Up_CSVL*cErr1Up_CSVL;
    else                 cErrDown_CSVL_2 += cErr1Up_CSVL*cErr1Up_CSVL;
    if( cErr2Up_CSVL>0 ) cErrUp_CSVL_2 += cErr2Up_CSVL*cErr2Up_CSVL;
    else                 cErrDown_CSVL_2 += cErr2Up_CSVL*cErr2Up_CSVL;
    if( cErr1Down_CSVL>0 ) cErrUp_CSVL_2 += cErr1Down_CSVL*cErr1Down_CSVL;
    else                   cErrDown_CSVL_2 += cErr1Down_CSVL*cErr1Down_CSVL;
    if( cErr2Down_CSVL>0 ) cErrUp_CSVL_2 += cErr2Down_CSVL*cErr2Down_CSVL;
    else                   cErrDown_CSVL_2 += cErr2Down_CSVL*cErr2Down_CSVL;

    double cErrUp_CSVL = sqrt( cErrUp_CSVL_2 );
    double cErrDown_CSVL = sqrt( cErrDown_CSVL_2 );

    //CSVM
    double cErr1Up_CSVM = hc_SF_cErr1Up->GetBinContent(binCSVM) - hc_SF->GetBinContent(binCSVM);
    double cErr1Down_CSVM = hc_SF_cErr1Down->GetBinContent(binCSVM) - hc_SF->GetBinContent(binCSVM);
    double cErr2Up_CSVM = hc_SF_cErr2Up->GetBinContent(binCSVM) - hc_SF->GetBinContent(binCSVM);
    double cErr2Down_CSVM = hc_SF_cErr2Down->GetBinContent(binCSVM) - hc_SF->GetBinContent(binCSVM);

    double cErrUp_CSVM_2 = 0, cErrDown_CSVM_2 = 0;
    if( cErr1Up_CSVM>0 ) cErrUp_CSVM_2   += cErr1Up_CSVM*cErr1Up_CSVM;
    else                 cErrDown_CSVM_2 += cErr1Up_CSVM*cErr1Up_CSVM;
    if( cErr2Up_CSVM>0 ) cErrUp_CSVM_2   += cErr2Up_CSVM*cErr2Up_CSVM;
    else                 cErrDown_CSVM_2 += cErr2Up_CSVM*cErr2Up_CSVM;
    if( cErr1Down_CSVM>0 ) cErrUp_CSVM_2 += cErr1Down_CSVM*cErr1Down_CSVM;
    else                   cErrDown_CSVM_2 += cErr1Down_CSVM*cErr1Down_CSVM;
    if( cErr2Down_CSVM>0 ) cErrUp_CSVM_2 += cErr2Down_CSVM*cErr2Down_CSVM;
    else                   cErrDown_CSVM_2 += cErr2Down_CSVM*cErr2Down_CSVM;

    double cErrUp_CSVM = sqrt( cErrUp_CSVM_2 );
    double cErrDown_CSVM = sqrt( cErrDown_CSVM_2 );


    //CSVT
    double cErr1Up_CSVT = hc_SF_cErr1Up->GetBinContent(binCSVT) - hc_SF->GetBinContent(binCSVT);
    double cErr1Down_CSVT = hc_SF_cErr1Down->GetBinContent(binCSVT) - hc_SF->GetBinContent(binCSVT);
    double cErr2Up_CSVT = hc_SF_cErr2Up->GetBinContent(binCSVT) - hc_SF->GetBinContent(binCSVT);
    double cErr2Down_CSVT = hc_SF_cErr2Down->GetBinContent(binCSVT) - hc_SF->GetBinContent(binCSVT);

    double cErrUp_CSVT_2 = 0, cErrDown_CSVT_2 = 0;
    if( cErr1Up_CSVT>0 ) cErrUp_CSVT_2 += cErr1Up_CSVT*cErr1Up_CSVT;
    else                 cErrDown_CSVT_2 += cErr1Up_CSVT*cErr1Up_CSVT;
    if( cErr2Up_CSVT>0 ) cErrUp_CSVT_2 += cErr2Up_CSVT*cErr2Up_CSVT;
    else                 cErrDown_CSVT_2 += cErr2Up_CSVT*cErr2Up_CSVT;
    if( cErr1Down_CSVT>0 ) cErrUp_CSVT_2 += cErr1Down_CSVT*cErr1Down_CSVT;
    else                   cErrDown_CSVT_2 += cErr1Down_CSVT*cErr1Down_CSVT;
    if( cErr2Down_CSVT>0 ) cErrUp_CSVT_2 += cErr2Down_CSVT*cErr2Down_CSVT;
    else                   cErrDown_CSVT_2 += cErr2Down_CSVT*cErr2Down_CSVT;

    double cErrUp_CSVT = sqrt( cErrUp_CSVT_2 );
    double cErrDown_CSVT = sqrt( cErrDown_CSVT_2 );



    if( verbose ){
    printf("==> CSVL: %s\t SF = %.2f +%.2f -%.2f (+%.3f %.2f)\n",hist_name[iHist].Data(), hc_SF->GetBinContent(binCSVL),
	   cErrUp_CSVL, cErrDown_CSVL, 2*(hc_SF_cUp->GetBinContent(binCSVL)-hc_SF->GetBinContent(binCSVL)), hc_SF_cDown->GetBinContent(binCSVL)-hc_SF->GetBinContent(binCSVL) );
    printf("==> CSVM: %s\t SF = %.2f +%.2f -%.2f (+%.3f %.2f)\n",hist_name[iHist].Data(), hc_SF->GetBinContent(binCSVM),
	   cErrUp_CSVM, cErrDown_CSVM, 2*(hc_SF_cUp->GetBinContent(binCSVM)-hc_SF->GetBinContent(binCSVM)), hc_SF_cDown->GetBinContent(binCSVM)-hc_SF->GetBinContent(binCSVM) );
    printf("==> CSVT: %s\t SF = %.2f +%.2f -%.2f (+%.3f %.2f)\n",hist_name[iHist].Data(), hc_SF->GetBinContent(binCSVT),
	   cErrUp_CSVT, cErrDown_CSVT, 2*(hc_SF_cUp->GetBinContent(binCSVT)-hc_SF->GetBinContent(binCSVT)), hc_SF_cDown->GetBinContent(binCSVT)-hc_SF->GetBinContent(binCSVT) );
    }


    TH1D* hc_SF_cUp_x2 = (TH1D*)hc_SF_cUp->Clone("hc_SF_cUp_x2_"+hist_name[iHist]);
    TH1D* hc_SF_cDown_x2 = (TH1D*)hc_SF_cDown->Clone("hc_SF_cDown_x2_"+hist_name[iHist]);

    for( int iBin=0; iBin<hc_SF_cUp_x2->GetNbinsX(); iBin++ ){
      double up_x2   = hc_SF->GetBinContent(iBin+1) + 2 * ( hc_SF_cUp->GetBinContent(iBin+1) - hc_SF->GetBinContent(iBin+1) );
      double down_x2 = hc_SF->GetBinContent(iBin+1) - 2 * ( hc_SF->GetBinContent(iBin+1) - hc_SF_cDown->GetBinContent(iBin+1) );

      hc_SF_cUp_x2->SetBinContent(iBin+1,up_x2);
      hc_SF_cDown_x2->SetBinContent(iBin+1,down_x2);
    }

    hc_SF->SetLineColor(kBlue);
    hc_SF_cUp_x2->SetLineColor(kGreen+1);
    hc_SF_cDown_x2->SetLineColor(kGreen+1);

    hc_SF->SetLineWidth(2);
    hc_SF_cUp_x2->SetLineWidth(2);
    hc_SF_cDown_x2->SetLineWidth(2);


    hc_SF->SetStats(0);
    hc_SF->GetYaxis()->SetRangeUser(0.86,1.14);
    hc_SF->GetXaxis()->SetRangeUser(0.1,1.00);

    hc_SF->SetTitle(";CSV;Charm Scale Factor");

    hc_SF->Draw("hist");
    hc_SF_cUp_x2->Draw("histsame");
    hc_SF_cDown_x2->Draw("histsame");

    btagSF->Draw("same");
    btagSFUp_x2->Draw("same");
    btagSFDown_x2->Draw("same");

    TLegend *legend6 = new TLegend(0.10,0.91,0.90,0.99);

    legend6->SetFillColor(kWhite);
    legend6->SetLineColor(kWhite);
    legend6->SetShadowColor(kWhite);
    legend6->SetTextFont(42);
    legend6->SetTextSize(0.035);

    legend6->SetNColumns(4);

    legend6->AddEntry(hc_SF,"c SF (ttH)","l");
    legend6->AddEntry(hc_SF_cUp_x2,"c SF Unc (ttH)","l");
    legend6->AddEntry(btagSF,"c SF (BTV)","l");
    legend6->AddEntry(btagSFUp_x2,"c SF Unc (BTV)","l");

    legend6->Draw();


    img = dirprefix + "cSF_totalErr_integralCompare_" + hist_name[iHist] + "_fit_All.png";
    if( makePlots ) c1->Print(img);
    img = dirprefix + "cSF_totalErr_integralCompare_" + hist_name[iHist] + "_fit_All.pdf";
    if( makePlots ) c1->Print(img);



    c_csv_ratio_final[iHist]->Write(Form("c_%s_final",hist_name[iHist].Data()));
    h_cErr1Up->Write(Form("c_%s_final_cErr1Up",hist_name[iHist].Data()));
    h_cErr1Down->Write(Form("c_%s_final_cErr1Down",hist_name[iHist].Data()));
    h_cErr2Up->Write(Form("c_%s_final_cErr2Up",hist_name[iHist].Data()));
    h_cErr2Down->Write(Form("c_%s_final_cErr2Down",hist_name[iHist].Data()));



    delete m_i;
    delete m_i_JESUp;
    delete m_i_JESDown;
    delete m_i_LFUp;
    delete m_i_LFDown;
    delete m_i_Stats1Up;
    delete m_i_Stats1Down;
    delete m_i_Stats2Up;
    delete m_i_Stats2Down;

    delete legend;
  }


  std::cout << " Done! " << std::endl;

  histofile.Write();
  histofile.Close();


}
