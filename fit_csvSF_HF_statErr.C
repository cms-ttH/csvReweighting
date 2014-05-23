#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"

#include "TRandom3.h"

#include "TCanvas.h"

#include "Math/MinimizerOptions.h"
#include <TSpline.h>
#include "Math/Interpolator.h"
#define MAXPOINTS 200

//______________________________________________________________________________
void fit_csvSF_HF_statErr( bool verbose_=false ){


  TH1::SetDefaultSumw2();

  TString dirprefix = "Images/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


  TFile *file = TFile::Open("HistoFiles/csv_rwt_hf_IT_v2.root");

  TString histofilename = "HistoFiles/fit_csvSF_HF_statErr.root";
  TFile histofile(histofilename,"recreate");
  histofile.cd();


  TRandom3 r;

  // //CSV original
  // int ncsvbins_CSV = 18;
  // double csvbins_CSV[] = { -0.04, 0.0, 0.122, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.723, 0.767, 0.810, 0.854, 0.898, 0.918, 0.939, 0.959, 0.980, 1.01 };

  //CSV new
  int ncsvbins = 18;
  double csvbins_new[] = { -0.04, 0.0, 0.122, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.723, 0.767, 0.810, 0.854, 0.898, 0.918, 0.939, 0.959, 0.980, 1.01 };


  std::vector<TString> hist_name;
  std::vector<TString> label;
  std::vector<double> x_center;
  std::vector<int> color;
  std::vector<TString> data_hist_name;
  std::vector<TString> mc_b_hist_name;
  std::vector<TString> mc_nonb_hist_name;

  int maxPt  = 5;
  int maxEta = 1;

  for( int iPt=0; iPt<maxPt; iPt++ ){
    for( int iEta=0; iEta<maxEta; iEta++ ){
      hist_name.push_back( Form("csv_ratio_Pt%d_Eta%d", iPt, iEta) );
      label.push_back( Form("Pt%d", iPt) );
      color.push_back( iPt+1 );
      x_center.push_back( iPt );

      data_hist_name.push_back( Form("csv_Data_Pt%d_Eta%d", iPt, iEta) );
      mc_b_hist_name.push_back( Form("csv_MC_bjets_Pt%d_Eta%d", iPt, iEta) );
      mc_nonb_hist_name.push_back( Form("csv_MC_nonbjets_Pt%d_Eta%d", iPt, iEta) );
    }
  }

  int NumHists = int( hist_name.size() );

  TH1D* h_csv_ratio[NumHists];

  TH1D* h_csv_data[NumHists];
  TH1D* h_csv_mc_b[NumHists];
  TH1D* h_csv_mc_nonb[NumHists];


  for( int iHist=0; iHist<NumHists; iHist++ ){
    TH1D* h_csv_data_temp1 = (TH1D*)file->Get( data_hist_name[iHist] )->Clone( Form("h_%s_temp1",data_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_b_temp1 = (TH1D*)file->Get( mc_b_hist_name[iHist] )->Clone( Form("h_%s_temp1",mc_b_hist_name[iHist].Data()) );
    TH1D* h_csv_mc_nonb_temp1 = (TH1D*)file->Get( mc_nonb_hist_name[iHist] )->Clone( Form("h_%s_temp1",mc_nonb_hist_name[iHist].Data()) );

    // rebin
    TH1D* h_csv_data_temp0 = (TH1D*)h_csv_data_temp1->Rebin( ncsvbins, Form("h_%s_temp1",data_hist_name[iHist].Data()), csvbins_new );
    TH1D* h_csv_mc_b_temp0 = (TH1D*)h_csv_mc_b_temp1->Rebin( ncsvbins, Form("h_%s_temp1",mc_b_hist_name[iHist].Data()), csvbins_new );
    TH1D* h_csv_mc_nonb_temp0 = (TH1D*)h_csv_mc_nonb_temp1->Rebin( ncsvbins, Form("h_%s_temp1",mc_nonb_hist_name[iHist].Data()), csvbins_new );

    h_csv_data[iHist] = new TH1D( Form("h_%s",data_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_b[iHist] = new TH1D( Form("h_%s",mc_b_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_nonb[iHist] = new TH1D( Form("h_%s",mc_nonb_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );

    for( int iBin=0; iBin<ncsvbins; iBin++ ){
      if( iBin==0 ){
	h_csv_data[iHist]->SetBinContent(iBin+1, h_csv_data_temp1->GetBinContent(iBin+1));
	h_csv_data[iHist]->SetBinError(iBin+1, h_csv_data_temp1->GetBinError(iBin+1));

	h_csv_mc_b[iHist]->SetBinContent(iBin+1, h_csv_mc_b_temp1->GetBinContent(iBin+1));
	h_csv_mc_b[iHist]->SetBinError(iBin+1, h_csv_mc_b_temp1->GetBinError(iBin+1));

	h_csv_mc_nonb[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_temp1->GetBinContent(iBin+1));
	h_csv_mc_nonb[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_temp1->GetBinError(iBin+1));
      }
      else {
	h_csv_data[iHist]->SetBinContent(iBin+1, h_csv_data_temp0->GetBinContent(iBin+1));
	h_csv_data[iHist]->SetBinError(iBin+1, h_csv_data_temp0->GetBinError(iBin+1));

	h_csv_mc_b[iHist]->SetBinContent(iBin+1, h_csv_mc_b_temp0->GetBinContent(iBin+1));
	h_csv_mc_b[iHist]->SetBinError(iBin+1, h_csv_mc_b_temp0->GetBinError(iBin+1));

	h_csv_mc_nonb[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_temp0->GetBinContent(iBin+1));
	h_csv_mc_nonb[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_temp0->GetBinError(iBin+1));
      }
    }


    h_csv_ratio[iHist]        = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_%d",iHist));

    h_csv_ratio[iHist]->Add(h_csv_mc_nonb[iHist],-1);
    h_csv_ratio[iHist]->Divide(h_csv_mc_b[iHist]);
  }

  int nBins = h_csv_ratio[0]->GetNbinsX();




  TH1D* h_csv_ratio_StatsUp[NumHists][nBins];
  TH1D* h_csv_ratio_StatsDown[NumHists][nBins];

  double threshold = 0.001;

  for( int iHist=0; iHist<NumHists; iHist++ ){

    for( int iBin=0; iBin<nBins; iBin++ ){
      double data = h_csv_data[iHist]->GetBinContent(iBin+1);
      double data_err = h_csv_data[iHist]->GetBinError(iBin+1);
	
      double mc_b = h_csv_mc_b[iHist]->GetBinContent(iBin+1);
      double mc_b_err = h_csv_mc_b[iHist]->GetBinError(iBin+1);
	
      double mc_nonb = h_csv_mc_nonb[iHist]->GetBinContent(iBin+1);
      double mc_nonb_err = h_csv_mc_nonb[iHist]->GetBinError(iBin+1);

      int iToy=0;
      bool found_match_StatsUp = false;
      bool found_match_StatsDown = false;
      while( !( found_match_StatsUp && found_match_StatsDown ) ){
	TH1D* h_csv_toys_data = (TH1D*)h_csv_data[iHist]->Clone(Form("h_%s_CSVbin%d_toy%d",data_hist_name[iHist].Data(),iBin,iToy));
	TH1D* h_csv_toys_mc_b = (TH1D*)h_csv_mc_b[iHist]->Clone(Form("h_%s_CSVbin%d_toy%d",mc_b_hist_name[iHist].Data(),iBin,iToy));
	TH1D* h_csv_toys_mc_nonb = (TH1D*)h_csv_mc_nonb[iHist]->Clone(Form("h_%s_CSVbin%d_toy%d",mc_nonb_hist_name[iHist].Data(),iBin,iToy));

	double new_data = r.Gaus(data,data_err);
	double new_mc_b = r.Gaus(mc_b,mc_b_err);
	double new_mc_nonb = r.Gaus(mc_nonb,mc_nonb_err);

	h_csv_toys_data->SetBinContent(iBin+1,new_data);
	h_csv_toys_mc_b->SetBinContent(iBin+1,new_mc_b);
	h_csv_toys_mc_nonb->SetBinContent(iBin+1,new_mc_nonb);

	double sum_data_toy = h_csv_toys_data->Integral();
	double sum_mc_toy = h_csv_toys_mc_b->Integral() + h_csv_toys_mc_nonb->Integral();

	h_csv_toys_mc_b->Scale( sum_data_toy / sum_mc_toy );
	h_csv_toys_mc_nonb->Scale( sum_data_toy / sum_mc_toy );

	TH1D* h_csv_toys_ratio_temp = (TH1D*)h_csv_toys_data->Clone(Form("h_csv_ratio_temp_ptbin%d_toy%d",iHist,iToy));
	h_csv_toys_ratio_temp->Add( h_csv_toys_mc_nonb, -1 );
	h_csv_toys_ratio_temp->Divide( h_csv_toys_mc_b );

	double sf = h_csv_ratio[iHist]->GetBinContent(iBin+1);
	double sf_err = h_csv_ratio[iHist]->GetBinError(iBin+1);
	double new_sf = h_csv_toys_ratio_temp->GetBinContent(iBin+1);

	if( fabs( (sf+sf_err-new_sf)/sf )<threshold && !found_match_StatsUp ){
	  h_csv_ratio_StatsUp[iHist][iBin] = (TH1D*)h_csv_toys_ratio_temp->Clone(Form("h_%s_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1));
	  found_match_StatsUp = true;
	  if( verbose_ ){
	    printf("\t StatsUp\t   iHist = %d,\t iBin = %d,\t SF = %.4f,\t SF_err = %.4f,\t new_SF = %.4f,\t (SF+SF_err-new_SF)/SF = %.4f \n",
		   iHist, iBin+1, sf, sf_err, new_sf, fabs( (sf+sf_err-new_sf)/sf ) );
	  }
	}
	else if( fabs( (sf-sf_err-new_sf)/sf )<threshold && !found_match_StatsDown ){
	  h_csv_ratio_StatsDown[iHist][iBin] = (TH1D*)h_csv_toys_ratio_temp->Clone(Form("h_%s_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1));
	  found_match_StatsDown = true;
	  if( verbose_ ){
	    printf("\t StatsDown\t iHist = %d,\t iBin = %d,\t SF = %.4f,\t SF_err = %.4f,\t new_SF = %.4f,\t (SF-SF_err-new_SF)/SF = %.4f \n",
		   iHist, iBin+1, sf, sf_err, new_sf, fabs( (sf-sf_err-new_sf)/sf ) );
	  }
	}

	iToy++;
      } // end loop over toys

      if( verbose_ ) std::cout << "\t Whew! It took " << iToy << " toys to find a match for everyone! " << std::endl;
    } // end loop over CSV bins
    // some verbosity to make sure everything isn't stuck
    std::cout << "\t On histogram " << iHist+1 << " / " << NumHists << std::endl;
  } // end loop over histograms


  TCanvas *c1 = new TCanvas("c1");

  for( int iHist=0; iHist<NumHists; iHist++ ){
    h_csv_ratio[iHist]->SetStats(0);
    h_csv_ratio[iHist]->SetLineWidth(2);
    h_csv_ratio[iHist]->SetMarkerStyle(20);

    for( int iBin=0; iBin<nBins; iBin++ ){
      h_csv_ratio[iHist]->Draw("pe1");
      h_csv_ratio_StatsUp[iHist][iBin]->Draw("histsame");
      h_csv_ratio_StatsDown[iHist][iBin]->Draw("histsame");

      c1->Print(Form("%shfSF_statsErrCompare_%s_CSVbin%d_Stats.png",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));
    }
  }



  histofile.Write();
  histofile.Close();

  std::cout << " Done! " << std::endl;

}

