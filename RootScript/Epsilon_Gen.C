#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TGraph.h"

class TVirtualPad;

void Epsilon(char* FileName)
{

  // TH1F *EpHist = new TH1F ("h1",120, 0, 60.0);

  TFile * TreeFile = TFile::Open(FileName);

  if (!TreeFile || TreeFile->IsZombie())
    { 
      return; 
    }

  Int_t Ne[120];
  Int_t Np[120];
  for (Int_t i = 0; i<120; i++)
    {
      Ne[i] = 0;
      Np[i] = 0;
    }
  Float_t E = 0.0;
  //  Float_t Ep = 0;

  TTreeReader Reader ("tree", TreeFile);
  TTreeReaderArray<Float_t> Charge(Reader, "monitors.charge");
  TTreeReaderArray<Float_t>  Energy(Reader, "monitoirs.energy_m");


    if ((!Charge.IsEmpty()) && (!Energe.IsEmpty()))
    {
      Int_t nDel = Charge.GetSize()
      for (Int_t iDel = 0; iDel < nDel; ++iDel)
	{
	  E = Energy[iDel];

	  if( Charge == -1)
	    {
	      for(Int_t i = 0; i<120 ;i++)
		{
		  if ((Float_t)i/2.0 <= E && (Float_t)(i+1)/2.0 > E)
		    {
		      Ne[i]++;
		      break;
		    }
		  else
		    continue;
		}
	    }
	  else
	    {
	      for(Int_t i = 0; i<120 ;i++)
		{
		  if ((Float_t)i/2.0 <= E && (Float_t)(i+1)/2.0 > E)
		  {
		      Np[i]++;
		      break;
		    }
		  else
		    continue;
		}
	    }
	}

    }

  Float_t Epsilon([120];
		  for (Int_t i = 0; i<120; i++) Epsilon[i] = 0.0;
		  for (Int_t i = 0; i<120; i++) 
		    {
		      Epsilon[i] = ((Ne[60+i]-Np[60-i])/(Ne[60+i]+Np[60-i]));
		    }


		  Int_t n = 120;
		  Float_t Del[n];
		  for (Int_t i = 0; i < n; i++)
		    {
		      Del[i] = -30 + 0.5*i;
		    }

		  TGraph *Gr = new TGraph (n, Del, Epsilon);
		  
		  TCanvas *Ca = new TCanvas("c1","Ep - Del Graph",200,10,600,400);
		  Gr->Draw("AC*");
}

		      


