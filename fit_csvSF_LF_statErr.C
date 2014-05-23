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
void fit_csvSF_LF_statErr( bool verbose_=false ){


  TH1::SetDefaultSumw2();

  TString dirprefix = "Images/";

  struct stat st;
  if( stat(dirprefix.Data(),&st) != 0 )  mkdir(dirprefix.Data(),0777);


  TFile *file = TFile::Open("HistoFiles/csv_rwt_lf_IT_v2.root");

  TString histofilename = "HistoFiles/fit_csvSF_LF_statErr.root";
  TFile histofile(histofilename,"recreate");
  histofile.cd();


  TRandom3 r;

  //CSV new
  int ncsvbins = 17;
  double csvbins[] = { -10.0, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };
  double csvbins_new[] = { -0.04, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };


  std::vector<TString> hist_name;
  std::vector<TString> label;
  std::vector<double> x_center;
  std::vector<int> color;
  std::vector<TString> data_hist_name;
  std::vector<TString> mc_b_hist_name;
  std::vector<TString> mc_nonb_hist_name;

  int maxPt  = 3;
  int maxEta = 3;

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
    TH1D* h_csv_data_temp0 = (TH1D*)h_csv_data_temp1->Rebin( ncsvbins, Form("h_%s_temp1",data_hist_name[iHist].Data()), csvbins );
    TH1D* h_csv_mc_b_temp0 = (TH1D*)h_csv_mc_b_temp1->Rebin( ncsvbins, Form("h_%s_temp1",mc_b_hist_name[iHist].Data()), csvbins );
    TH1D* h_csv_mc_nonb_temp0 = (TH1D*)h_csv_mc_nonb_temp1->Rebin( ncsvbins, Form("h_%s_temp1",mc_nonb_hist_name[iHist].Data()), csvbins );

    h_csv_data[iHist] = new TH1D( Form("h_%s",data_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_b[iHist] = new TH1D( Form("h_%s",mc_b_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );
    h_csv_mc_nonb[iHist] = new TH1D( Form("h_%s",mc_nonb_hist_name[iHist].Data()), ";CSV", ncsvbins, csvbins_new );

    for( int iBin=0; iBin<ncsvbins; iBin++ ){
      h_csv_data[iHist]->SetBinContent(iBin+1, h_csv_data_temp0->GetBinContent(iBin+1));
      h_csv_data[iHist]->SetBinError(iBin+1, h_csv_data_temp0->GetBinError(iBin+1));

      h_csv_mc_b[iHist]->SetBinContent(iBin+1, h_csv_mc_b_temp0->GetBinContent(iBin+1));
      h_csv_mc_b[iHist]->SetBinError(iBin+1, h_csv_mc_b_temp0->GetBinError(iBin+1));

      h_csv_mc_nonb[iHist]->SetBinContent(iBin+1, h_csv_mc_nonb_temp0->GetBinContent(iBin+1));
      h_csv_mc_nonb[iHist]->SetBinError(iBin+1, h_csv_mc_nonb_temp0->GetBinError(iBin+1));
    }


    h_csv_ratio[iHist] = (TH1D*)h_csv_data[iHist]->Clone(Form("h_csv_ratio_%d",iHist));

    h_csv_ratio[iHist]->Add(h_csv_mc_b[iHist],-1);
    h_csv_ratio[iHist]->Divide(h_csv_mc_nonb[iHist]);
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
	h_csv_toys_ratio_temp->Add( h_csv_toys_mc_b, -1 );
	h_csv_toys_ratio_temp->Divide( h_csv_toys_mc_nonb );

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

      if( verbose_ ) std::cout << "\t Stats Test, iHist = " << iHist << ",\t iBin = " << iBin+1 << ",\t\t It took " << iToy << " toys to find a match for everyone!" << std::endl;
      //printf("\t Stats Test, iHist = %d,\t iBin = %d,\t\t It took %d toys to find a match for everyone!\n",iHist, iBin+1, iToy );
    } // end loop over CSV bins
    // some verbosity to make sure everything isn't stuck
    std::cout << "\t On histogram " << iHist+1 << " / " << NumHists << std::endl;

  } // end loop over histograms


  TCanvas *c1 = new TCanvas("c1");

  //// pol6
  TF1* f0 = new TF1("f0","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",0,0.95 );
  TF1* f0_StatsUp = new TF1("f0_StatsUp","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",0,0.95 );
  TF1* f0_StatsDown = new TF1("f0_StatsDown","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x",0,0.95 );

  int nPars = f0->GetNpar();
  for( int iPar=0; iPar<nPars; iPar++ ){
    f0->SetParameter(iPar,1.);
    f0_StatsUp->SetParameter(iPar,1.);
    f0_StatsDown->SetParameter(iPar,1.);
  }

  f0->SetLineColor(kBlack);
  f0_StatsUp->SetLineColor(kRed);
  f0_StatsDown->SetLineColor(kBlack);

  f0->SetLineWidth(2);
  f0_StatsUp->SetLineWidth(2);
  f0_StatsDown->SetLineWidth(2);


  TH1D* h_csv_ratio_final_StatsUp[NumHists][nBins];
  TH1D* h_csv_ratio_final_StatsDown[NumHists][nBins];

  int NumFinalBins = 1000;
  double firstPoint = 0;
  double lastPoint  = 0.95;

  for( int iHist=0; iHist<NumHists; iHist++ ){
    h_csv_ratio[iHist]->SetStats(0);
    h_csv_ratio[iHist]->SetLineWidth(2);
    h_csv_ratio[iHist]->SetMarkerStyle(20);

    TFitResultPtr res = h_csv_ratio[iHist]->Fit(f0,"+mrN0S");

    for( int iBin=0; iBin<nBins; iBin++ ){
      h_csv_ratio_final_StatsUp[iHist][iBin] = new TH1D( Form("h_%s_final_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1), ";CSV", NumFinalBins, -0.04, 1.01 );
      h_csv_ratio_final_StatsDown[iHist][iBin] = new TH1D( Form("h_%s_final_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1), ";CSV", NumFinalBins, -0.04, 1.01 );

      h_csv_ratio_StatsUp[iHist][iBin]->SetLineWidth(2);
      h_csv_ratio_StatsDown[iHist][iBin]->SetLineWidth(2);

      h_csv_ratio_StatsUp[iHist][iBin]->SetLineColor(kRed);
      h_csv_ratio_StatsDown[iHist][iBin]->SetLineColor(kBlack);

      h_csv_ratio[iHist]->GetXaxis()->SetRangeUser(-0.041,0.9489);
      h_csv_ratio[iHist]->Draw("pe1");
      h_csv_ratio_StatsUp[iHist][iBin]->Draw("histsame");
      h_csv_ratio_StatsDown[iHist][iBin]->Draw("histsame");

      c1->Print(Form("%slfSF_statsErrCompare_%s_CSVbin%d_Stats.png",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));

      TFitResultPtr res_StatsUp = h_csv_ratio_StatsUp[iHist][iBin]->Fit(f0_StatsUp,"+mrN0S");
      printf("\t StatsUp:   iHist = %s,\t iBin = %d,\t Chisquare = %.2f,\t NDF = %d,\t Prob = %.4f \n",
	     hist_name[iHist].Data(), iBin, f0_StatsUp->GetChisquare(), f0_StatsUp->GetNDF(), f0_StatsUp->GetProb() );

      TFitResultPtr res_StatsDown = h_csv_ratio_StatsDown[iHist][iBin]->Fit(f0_StatsDown,"+mrN0S");
      printf("\t StatsDown: iHist = %s,\t iBin = %d,\t Chisquare = %.2f,\t NDF = %d,\t Prob = %.4f \n",
	     hist_name[iHist].Data(), iBin, f0_StatsDown->GetChisquare(), f0_StatsDown->GetNDF(), f0_StatsDown->GetProb() );


      for( int jBin=0; jBin<NumFinalBins; jBin++ ){
	double center = h_csv_ratio_final_StatsUp[iHist][iBin]->GetBinCenter(jBin+1);
	if( center<firstPoint ){
	  if( iBin==0 ){
	    h_csv_ratio_final_StatsUp[iHist][iBin]->SetBinContent(jBin+1,h_csv_ratio[iHist]->GetBinContent(1) + h_csv_ratio[iHist]->GetBinError(1));
	    h_csv_ratio_final_StatsDown[iHist][iBin]->SetBinContent(jBin+1,h_csv_ratio[iHist]->GetBinContent(1) - h_csv_ratio[iHist]->GetBinError(1));
	  }
	  else {
	    h_csv_ratio_final_StatsUp[iHist][iBin]->SetBinContent(jBin+1,h_csv_ratio_StatsUp[iHist][iBin]->GetBinContent(1));
	    h_csv_ratio_final_StatsDown[iHist][iBin]->SetBinContent(jBin+1,h_csv_ratio_StatsDown[iHist][iBin]->GetBinContent(1));
	  }
	}
	else if( center>lastPoint ){
	  h_csv_ratio_final_StatsUp[iHist][iBin]->SetBinContent(jBin+1,f0_StatsUp->Eval( lastPoint ));
	  h_csv_ratio_final_StatsDown[iHist][iBin]->SetBinContent(jBin+1,f0_StatsDown->Eval( lastPoint ));
	}
	else {
	  h_csv_ratio_final_StatsUp[iHist][iBin]->SetBinContent(jBin+1,f0_StatsUp->Eval( center ));
	  h_csv_ratio_final_StatsDown[iHist][iBin]->SetBinContent(jBin+1,f0_StatsDown->Eval( center ));
	}
      }

      h_csv_ratio_final_StatsUp[iHist][iBin]->SetLineColor(kRed);
      h_csv_ratio_final_StatsDown[iHist][iBin]->SetLineColor(kBlue);

      h_csv_ratio_final_StatsUp[iHist][iBin]->SetLineWidth(2);
      h_csv_ratio_final_StatsDown[iHist][iBin]->SetLineWidth(2);


      h_csv_ratio[iHist]->Draw("pe1");
      h_csv_ratio_final_StatsUp[iHist][iBin]->Draw("histsame");
      h_csv_ratio_final_StatsDown[iHist][iBin]->Draw("histsame");

      c1->Print(Form("%slfSF_statsErrCompareFunctions_%s_CSVbin%d_Stats.png",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));
      c1->Print(Form("%slfSF_statsErrCompareFunctions_%s_CSVbin%d_Stats.pdf",dirprefix.Data(),hist_name[iHist].Data(),iBin+1));

      h_csv_ratio_final_StatsUp[iHist][iBin]->Write(Form("%s_final_CSVbin%d_StatsUp",hist_name[iHist].Data(),iBin+1));
      h_csv_ratio_final_StatsDown[iHist][iBin]->Write(Form("%s_final_CSVbin%d_StatsDown",hist_name[iHist].Data(),iBin+1));
    }
  }

  std::cout << " Writing and closing the histogram file. " << std::endl;

  histofile.Write();
  histofile.Close();

  std::cout << " Done! " << std::endl;

}
