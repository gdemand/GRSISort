//g++ MakeMatrices.C -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries -lGRSIFormat -lAnalysisTreeBuilder -lGriffin -lSceptar -lDescant -lPaces -lGRSIDetector -lTigress -lSharc -lCSM -lTriFoil -lTGRSIint -lGRSILoop -lMidasFormat -lGRSIRootIO -lDataParser -lMidasFormat -lXMLParser -lXMLIO -lProof -lGuiHtml `root-config --cflags --libs` -lTreePlayer -o MakeMatrices

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>
#include <string>
#include <iterator>
#include <algorithm>
#include <climits>

#include "TROOT.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TVirtualIndex.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TGRSIRunInfo.h"

#ifndef __CINT__
#include "TGriffin.h"
#include "TSceptar.h"
#include "TPaces.h"
#include "TDescant.h"
#endif

#include "TDescantHit.h"

int checkSceptarHit(TSceptarHit* hit);
int checkDescantHit(TDescantHit* hit);

size_t descantlength = 0;

std::vector<unsigned int> cfddelayrange = {47, 47};
std::vector<unsigned int> cfdattenuationrange = {11, 11};

unsigned int cfddelayrangeentries = cfddelayrange.at(1) - cfddelayrange.at(0) + 1;
unsigned int cfdattenuationrangeentries = cfdattenuationrange.at(1) - cfdattenuationrange.at(0) + 1;

#ifdef partialsumcfd

//std::vector<unsigned int> zcdelayrange = {2, 30};
//std::vector<unsigned int> zcattenuationrange = {2, 30};

std::vector<unsigned int> zcdelayrange = {50, 200};
std::vector<unsigned int> zcattenuationrange = {30, 50};

unsigned int zcdelayrangeentries = zcdelayrange.at(1) - zcdelayrange.at(0) + 1;
unsigned int zcattenuationrangeentries = zcattenuationrange.at(1) - zcattenuationrange.at(0) + 1;


#else
std::vector<unsigned int> zcdiffrcrange = {20, 220};
std::vector<unsigned int> zcintrc1range = {50, 250};
std::vector<unsigned int> zcintrc2range = {50, 250};

unsigned int zcdiffrcstep = 50;
unsigned int zcintrc1step = 50;
unsigned int zcintrc2step = 50;

//std::vector<unsigned int> zcdiffrcrange = {30, 100};
//std::vector<unsigned int> zcintrc1range = {100, 300};
//std::vector<unsigned int> zcintrc2range = {100, 300};

//unsigned int zcdiffrcstep = 20;
//unsigned int zcintrc1step = 20;
//unsigned int zcintrc2step = 20;


unsigned int zcdiffrcrangeentries = (zcdiffrcrange.at(1) - zcdiffrcrange.at(0))/zcdiffrcstep + 1;
unsigned int zcintrc1rangeentries = (zcintrc1range.at(1) - zcintrc1range.at(0))/zcintrc1step + 1;
unsigned int zcintrc2rangeentries = (zcintrc2range.at(1) - zcintrc2range.at(0))/zcintrc2step + 1;

#endif

#ifdef __CINT__
void MakeMatrices() {
   if(!AnalysisTree) {
      printf("No analysis tree found!\n");
      return;
   }
   
   //coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default is 4000)
   TList *list = MakeMatrices(AnalysisTree, 0., 30., 80., 6000, 0., 6000.);
   
   std::string fileName = gFile->GetName();
   if(fileName.find_last_of("/") != std::string::npos) {
      fileName.insert(fileName.find_last_of("/")+1,"matrices_");
   } else {
      fileName.insert(0,"matrices_");
   }
   TFile *outfile = new TFile(file.c_str(),"recreate");
   list->Write();
}
#endif


TList *MakeMatrices(TTree* tree, std::vector<Double_t> &cfdStdDev, std::vector<Double_t> &zcMean, std::vector<Double_t> &zcStdDev, int coincLow = 0, int coincHigh = 10, int bg = 100, int nofBins = 4000, double low = 0., double high = 4000., long maxEntries = 0, TStopwatch* w = NULL) {
   if(w == NULL) {
      w = new TStopwatch;
      w->Start();
   }
   TList* list = new TList;
   
   const size_t MEM_SIZE = (size_t)1024*(size_t)1024*(size_t)1024*(size_t)8; // 8 GB
   
   unsigned int interpolationsteps = 256;
   
#ifndef partialsums
   Double_t zc;
#endif

#ifdef partialsumcfd
   std::vector<TH1F*> descantCfdArray(cfddelayrangeentries*cfdattenuationrangeentries*zcdelayrangeentries*zcattenuationrangeentries,0);
   std::vector<TH1F*> descantZcArray(cfddelayrangeentries*cfdattenuationrangeentries*zcdelayrangeentries*zcattenuationrangeentries,0);

   zcMean.resize(cfddelayrangeentries*cfdattenuationrangeentries*zcdelayrangeentries*zcattenuationrangeentries,0.);
   zcStdDev.resize(cfddelayrangeentries*cfdattenuationrangeentries*zcdelayrangeentries*zcattenuationrangeentries,0.);
#else
   std::vector<TH1F*> descantCfdArray(cfddelayrangeentries*cfdattenuationrangeentries*zcdiffrcrangeentries*zcintrc1rangeentries*zcintrc2rangeentries,0);
   std::vector<TH1F*> descantZcArray(cfddelayrangeentries*cfdattenuationrangeentries*zcdiffrcrangeentries*zcintrc1rangeentries*zcintrc2rangeentries,0);

   zcMean.resize(cfddelayrangeentries*cfdattenuationrangeentries*zcdiffrcrangeentries*zcintrc1rangeentries*zcintrc2rangeentries,0.);
   zcStdDev.resize(cfddelayrangeentries*cfdattenuationrangeentries*zcdiffrcrangeentries*zcintrc1rangeentries*zcintrc2rangeentries,0.);
#endif

   cfdStdDev.resize(cfddelayrangeentries*cfdattenuationrangeentries,0.);

   
   unsigned int currententry = 0;
   unsigned int currentcfdentry = 0;
   for(unsigned int cfddelay = cfddelayrange[0]; cfddelay <= cfddelayrange[1]; ++cfddelay) {
      for(unsigned int cfdattenuation = cfdattenuationrange[0]; cfdattenuation <= cfdattenuationrange[1]; ++cfdattenuation) {
         std::string title = "descantCfd"  + std::to_string(cfddelay) + "_" + std::to_string(cfdattenuation);
         std::string labels = "descant cfd:, cfddelay: " + std::to_string(cfddelay) + ", cfdattenuation: " + std::to_string(cfdattenuation) + ";time[0.2 ns]";
         descantCfdArray.at(currentcfdentry) = new TH1F(title.c_str(), labels.c_str(), 500, -350, 450); list->Add(descantCfdArray.at(currentcfdentry));
//         descantCfdArray.at(currentcfdentry) = new TH1F(title.c_str(), labels.c_str(), 500, 1600, 2100); list->Add(descantCfdArray.at(currentcfdentry));
         ++currentcfdentry;
#ifdef partialsumcfd
         for(unsigned int zcdelay = zcdelayrange[0]; zcdelay <= zcdelayrange[1]; ++zcdelay) {
            for(unsigned int zcattenuation = zcattenuationrange[0]; zcattenuation <= zcattenuationrange[1]; ++zcattenuation) {
               std::string title = "descantZc"  + std::to_string(cfddelay) + "_" + std::to_string(cfdattenuation) + "_" + std::to_string(zcdelay) + "_" + std::to_string(zcattenuation);
               std::string labels = "descant zc:, cfddelay: " + std::to_string(cfddelay) + ", cfdattenuation: " + std::to_string(cfdattenuation) + "zcdelay: " + std::to_string(zcdelay) + ", zcattenuation: " + std::to_string(zcattenuation) + ";time[0.2 ns]";
               descantZcArray.at(currententry) = new TH1F(title.c_str(), labels.c_str(), 500, 100, 1200); list->Add(descantZcArray.at(currententry));
               ++currententry;
            }
         }
#else
         for(unsigned int zcdiffrc = zcdiffrcrange[0]; zcdiffrc <= zcdiffrcrange[1]; zcdiffrc += zcdiffrcstep) {
            for(unsigned int zcintrc1 = zcintrc1range[0]; zcintrc1 <= zcintrc1range[1]; zcintrc1 += zcintrc1step) {
               for(unsigned int zcintrc2 = zcintrc2range[0]; zcintrc2 <= zcintrc2range[1]; zcintrc2 += zcintrc2step) {
                  std::string title = "descantZc"  + std::to_string(cfddelay) + "_" + std::to_string(cfdattenuation) + "_" + std::to_string(zcdiffrc) + "_" + std::to_string(zcintrc1) + "_" + std::to_string(zcintrc2);
                  std::string labels = "descant zc:, cfddelay: " + std::to_string(cfddelay) + ", cfdattenuation: " + std::to_string(cfdattenuation) + "zcdiffrc: " + std::to_string(zcdiffrc) + ", zcintrc1: " + std::to_string(zcintrc1) + ", zcintrc2: " + std::to_string(zcintrc2) + ";time[0.2 ns]";
                  descantZcArray.at(currententry) = new TH1F(title.c_str(), labels.c_str(), 500, 100, 1200); list->Add(descantZcArray.at(currententry));
                  ++currententry;
               }
            }
         }
         
#endif

      }
   }
   
   //even spaced binning in x and y, 1 keV bin width
   int nof3Dbins = 1500;
   double* xBins = new double[nof3Dbins+1];
   double* yBins = new double[nof3Dbins+1];
   for(int i = 0; i <= nof3Dbins; ++i) {
      xBins[i] = (double) i;
      yBins[i] = (double) i;
   }
   
   
   //set up branches
   TDescant* desc = 0;
   
   bool gotDescant;
   if(tree->FindBranch("TDescant") == 0) {
      gotDescant = false;
   } else {
      tree->SetBranchAddress("TDescant", &desc);
      gotDescant = true;
   }
   
   long entries = tree->GetEntries();
   
   std::vector<Short_t> temp_monitor;
   Double_t tempcfd;
   
   int entry;

   if(maxEntries == 0 || maxEntries > tree->GetEntries()) {
      maxEntries = tree->GetEntries();
   }
   for(entry = 0; entry < maxEntries; ++entry) {
      
      tree->GetEntry(entry);
      
      if(gotDescant && desc->GetMultiplicity() >= 1) {
         for(int b = 0; b < desc->GetMultiplicity(); ++b) {
            if(checkDescantHit(desc->GetDescantHit(b)) == -1) {
               continue;
            }
            //We have a good event!
            currententry = 0;
            currentcfdentry = 0;
            for(unsigned int cfddelay = cfddelayrange[0]; cfddelay <= cfddelayrange[1]; ++cfddelay) {
               for(unsigned int cfdattenuation = cfdattenuationrange[0]; cfdattenuation <= cfdattenuationrange[1]; ++cfdattenuation) {
                  
                  tempcfd = desc->GetDescantHit(b)->CalculateCfd(cfdattenuation/64., cfddelay, 0, interpolationsteps)/256.;
                  descantCfdArray.at(currentcfdentry)->Fill(tempcfd);
#ifdef partialsumcfd
                  for(unsigned int zcdelay = zcdelayrange[0]; zcdelay <= zcdelayrange[1]; ++zcdelay) {
                     for(unsigned int zcattenuation = zcattenuationrange[0]; zcattenuation <= zcattenuationrange[1]; ++zcattenuation) {
                        //Make zc plot
                        descantZcArray.at(currententry)->Fill(desc->GetDescantHit(b)->CalculateCfdForPartialSums(zcattenuation/64., zcdelay, interpolationsteps)/256.-tempcfd);
                        currententry++;
                     }
                  }
#else
                  for(unsigned int zcdiffrc = zcdiffrcrange[0]; zcdiffrc <= zcdiffrcrange[1]; zcdiffrc += zcdiffrcstep) {
                     for(unsigned int zcintrc1 = zcintrc1range[0]; zcintrc1 <= zcintrc1range[1]; zcintrc1 += zcintrc1step) {
                        for(unsigned int zcintrc2 = zcintrc2range[0]; zcintrc2 <= zcintrc2range[1]; zcintrc2 += zcintrc2step) {
                        //Make zc plot
                           
                           zc = TDescantHit::CalculateZeroCrossing(TDescantHit::Integrator(TDescantHit::Integrator(TDescantHit::Differentiator(TDescantHit::ShortVectorToDouble(*(desc->GetDescantHit(b)->GetWaveform())), zcdiffrc), zcintrc1), zcintrc2), interpolationsteps)/256.;

//                           std::cout << "Zc: " << zc << std::endl;
                           descantZcArray.at(currententry)->Fill(zc-tempcfd);
                           currententry++;
                        }
                     }
                  }
#endif
                  currentcfdentry++;
               }
            }
            
         }
      }
      if(entry%250 == 0) {
         std::cout << "\t" << entry << " / " << entries << " = "<< (float)entry/entries*100.0 << "%. " << w->RealTime() << " seconds" << "\r" << std::flush;
         w->Continue();
      }
   }
   
   currententry = 0;
   currentcfdentry = 0;
   //Store means and standard deviations - note GetRMS() returns a standard deviation
   for(unsigned int cfddelay = cfddelayrange[0]; cfddelay <= cfddelayrange[1]; ++cfddelay) {
      for(unsigned int cfdattenuation = cfdattenuationrange[0]; cfdattenuation <= cfdattenuationrange[1]; ++cfdattenuation) {
//         std::cout << "Cfddelay: " << cfddelay << ", cfdattenuation: " << cfdattenuation << std::endl;
         
         cfdStdDev.at(currentcfdentry)=descantCfdArray.at(currentcfdentry)->GetRMS();
         currentcfdentry++;
#ifdef partialsumcfd
         for(unsigned int zcdelay = zcdelayrange[0]; zcdelay <= zcdelayrange[1]; ++zcdelay) {
            for(unsigned int zcattenuation = zcattenuationrange[0]; zcattenuation <= zcattenuationrange[1]; ++zcattenuation) {
               //Make zc plot
               zcMean.at(currententry) = descantZcArray.at(currententry)->GetMean();
               zcStdDev.at(currententry) = descantZcArray.at(currententry)->GetRMS();
               currententry++;
            }
         }
#else
         for(unsigned int zcdiffrc = zcdiffrcrange[0]; zcdiffrc <= zcdiffrcrange[1]; zcdiffrc += zcdiffrcstep) {
            for(unsigned int zcintrc1 = zcintrc1range[0]; zcintrc1 <= zcintrc1range[1]; zcintrc1 += zcintrc1step) {
               for(unsigned int zcintrc2 = zcintrc2range[0]; zcintrc2 <= zcintrc2range[1]; zcintrc2 += zcintrc2step) {
                  //Make zc plot
                  zcMean.at(currententry) = descantZcArray.at(currententry)->GetMean();
                  zcStdDev.at(currententry) = descantZcArray.at(currententry)->GetRMS();
                  currententry++;
               }
            }
         }
#endif
      }
   }
   
   for(size_t i = 0; i < cfdStdDev.size() ; i++) {
      std::cout << "cfdStdDev[" <<i <<"]: " << cfdStdDev.at(i) << std::endl;
   }

   
   //Find best figure of merit from matrix and report settings
   //Generate slices at various parameters
   
   std::cout << "\t" << entry << " / " << entries << " = "<< (float)entry/entries*100.0 << "%. " << w->RealTime() << " seconds" << std::endl;
   w->Continue();
   
   //create all background corrected matrices
   
#ifdef __CINT__
   TCanvas* c = new TCanvas;
   c->Divide(2,2);
   c->cd(1);
   matrix->Draw("colz");
   c->cd(2);
   matrix_bgcorr->Draw("colz");
   c->cd(3);
   matrix_coinc->Draw("colz");
   c->cd(4);
   matrix_bg->Draw("colz");
   TCanvas* c2 = new TCanvas;
   c2->Divide(2);
   c2->cd(1);
   timeDiff->Draw();
   timeEGate->SetLineColor(2);
   timeEGate->Draw("same");
   timeModuleOne->SetLineColor(4);
   timeModuleOne->Draw("same");
   c2->cd(2);
   grifMult->Draw();
   grifMultCut->SetLineColor(2);
   grifMultCut->Draw("same");
#endif
   
   list->Sort();
   std::cout << "creating histograms done after " << w->RealTime() << " seconds" << std::endl;
   w->Continue();
   return list;
}


#ifndef __CINT__

int main(int argc, char **argv) {
   if(argc != 4 && argc != 3) {
      printf("try again (usage: %s <neutron analysis tree file> <gamma analysis tree file> <max entries>).\n",argv[0]);
      return 0;
   }
   
   std::vector<Double_t> neutronCfdStdDev;
   std::vector<Double_t> neutronZcMean;
   std::vector<Double_t> neutronZcStdDev;
   
   std::vector<Double_t> gammaCfdStdDev;
   std::vector<Double_t> gammaZcMean;
   std::vector<Double_t> gammaZcStdDev;
   
   TStopwatch w;
   w.Start();
   
   std::string neutronFileName;
   neutronFileName = argv[1];
   if(neutronFileName.find_last_of("/") != std::string::npos) {
      neutronFileName.insert(neutronFileName.find_last_of("/")+1,"optimize_");
   } else {
      neutronFileName.insert(0,"optimize_");
   }
   
   TFile* neutronFile = new TFile(argv[1]);
   if(neutronFile == NULL || !neutronFile->IsOpen()) {
      printf("Failed to open file '%s'!\n",argv[1]);
      return 1;
   }
   
   TTree* neutronTree = (TTree*) neutronFile->Get("AnalysisTree");
   
   if(neutronTree == NULL) {
      printf("Failed to find analysis tree in file '%s'!\n",argv[1]);
      return 1;
   }
   
   std::cout << argv[0] << ": starting MakeMatrices after " << w.RealTime() << " seconds" << std::endl;
   w.Continue();
   
   TChannel::ReadCalFromTree(neutronTree);
   
   //coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default is 4000)
   TList *neutronList;
   if(argc < 4) {
      neutronList = MakeMatrices(neutronTree, neutronCfdStdDev, neutronZcMean, neutronZcStdDev , -20., 40., 80., 6000, 0., 6000., 0, &w);
   } else {
      int entries = atoi(argv[4]);
      std::cout<<"Limiting processing of analysis tree to "<<entries<<" entries!"<<std::endl;
      neutronList = MakeMatrices(neutronTree, neutronCfdStdDev, neutronZcMean, neutronZcStdDev, 0., 30., 80., 6000, 0., 6000., entries, &w);
   }
   
   TFile *neutronOutfile = new TFile(neutronFileName.c_str(),"recreate");
   neutronList->Write();
   neutronOutfile->Close();
   neutronFile->Close();
   
   std::string gammaFileName;
   gammaFileName = argv[2];
   if(gammaFileName.find_last_of("/") != std::string::npos) {
      gammaFileName.insert(gammaFileName.find_last_of("/")+1,"optimize_");
   } else {
      gammaFileName.insert(0,"optimize_");
   }

   TFile* gammaFile = new TFile(argv[2]);
   if(gammaFile == NULL || !gammaFile->IsOpen()) {
      printf("Failed to open file '%s'!\n",argv[1]);
      return 1;
   }
   
   TTree* gammaTree = (TTree*) gammaFile->Get("AnalysisTree");
   
   if(gammaTree == NULL) {
      printf("Failed to find analysis tree in file '%s'!\n",argv[1]);
      return 1;
   }
   
   std::cout << argv[0] << ": starting MakeMatrices after " << w.RealTime() << " seconds" << std::endl;
   w.Continue();
   
   TChannel::ReadCalFromTree(gammaTree);
   
   //coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default is 4000)
   TList *gammaList;
   if(argc < 4) {
      gammaList = MakeMatrices(gammaTree, gammaCfdStdDev, gammaZcMean, gammaZcStdDev, -20., 40., 80., 6000, 0., 6000., 0, &w);
   } else {
      int entries = atoi(argv[4]);
      std::cout<<"Limiting processing of analysis tree to "<<entries<<" entries!"<<std::endl;
      gammaList = MakeMatrices(gammaTree, gammaCfdStdDev, gammaZcMean, gammaZcStdDev, 0., 30., 80., 6000, 0., 6000., entries, &w);
   }
   
   TFile *gammaOutfile = new TFile(gammaFileName.c_str(),"recreate");
   gammaList->Write();
   gammaOutfile->Close();
   gammaFile->Close();
   
   unsigned int currententry = 0;
   unsigned int currentcfdentry = 0;
   
   double tempfom = 0.;
   
   double bestneutroncfdstddev = INT_MAX;
   double bestgammacfdstddev = INT_MAX;
   
   unsigned int bestneutroncfddelay = 0;
   unsigned int bestneutroncfdattenuation=0;
   unsigned int bestgammacfddelay = 0;
   unsigned int bestgammacfdattenuation=0;
   
   double bestfom = 0.;
   unsigned int bestzccfddelay = 0;
   unsigned int bestzccfdattenuation = 0;

#ifdef partialsumcfd
   unsigned int bestzcdelay = 0;
   unsigned int bestzcattenuation = 0;
#else
   unsigned int bestzcdiffrc = 0;
   unsigned int bestzcintrc1 = 0;
   unsigned int bestzcintrc2 = 0;
#endif
   for(size_t i = 0; i < gammaCfdStdDev.size() ; i++) {
      std::cout << "gammaCfdStdDev[" <<i <<"]: " <<gammaCfdStdDev.at(i) << std::endl;
   }
   
   for(unsigned int cfddelay = cfddelayrange[0]; cfddelay <= cfddelayrange[1]; ++cfddelay) {
      for(unsigned int cfdattenuation = cfdattenuationrange[0]; cfdattenuation <= cfdattenuationrange[1]; ++cfdattenuation) {
         
         std::cout << "For cfddelay : " << cfddelay << ", and cfdattenuation: " << cfdattenuation << ", neutron std = " << neutronCfdStdDev.at(currentcfdentry) << " and gamma std = " << gammaCfdStdDev.at(currentcfdentry) << std::endl;
         
         if(neutronCfdStdDev.at(currentcfdentry) < bestneutroncfdstddev) {
            bestneutroncfdstddev = neutronCfdStdDev.at(currentcfdentry);
            bestneutroncfddelay = cfddelay;
            bestneutroncfdattenuation = cfdattenuation;
         }
         
         if(gammaCfdStdDev.at(currentcfdentry) < bestgammacfdstddev) {
            bestgammacfdstddev = gammaCfdStdDev.at(currentcfdentry);
            bestgammacfddelay = cfddelay;
            bestgammacfdattenuation = cfdattenuation;
         }
         ++currentcfdentry;
#ifdef partialsumcfd
         for(unsigned int zcdelay = zcdelayrange[0]; zcdelay <= zcdelayrange[1]; ++zcdelay) {
            for(unsigned int zcattenuation = zcattenuationrange[0]; zcattenuation <= zcattenuationrange[1]; ++zcattenuation) {
               tempfom = (neutronZcMean.at(currententry) - gammaZcMean.at(currententry))/(neutronZcStdDev.at(currententry) + gammaZcStdDev.at(currententry));
               
               std::cout << "For cfddelay : " << cfddelay << ", cfdattenuation: " << cfdattenuation << ", zcdelay: " << zcdelay << ", zcattenuation: " << zcattenuation << ", zc FOM = " << tempfom << std::endl;
               
               if(tempfom > bestfom) {
                  bestfom = tempfom;
                  bestzccfddelay = cfddelay;
                  bestzccfdattenuation = cfdattenuation;
                  bestzcdelay = zcdelay;
                  bestzcattenuation = zcattenuation;
               }
               
               ++currententry;
               
            }
         }
#else
         for(unsigned int zcdiffrc = zcdiffrcrange[0]; zcdiffrc <= zcdiffrcrange[1]; zcdiffrc += zcdiffrcstep) {
            for(unsigned int zcintrc1 = zcintrc1range[0]; zcintrc1 <= zcintrc1range[1]; zcintrc1 += zcintrc1step) {
               for(unsigned int zcintrc2 = zcintrc2range[0]; zcintrc2 <= zcintrc2range[1]; zcintrc2 += zcintrc2step) {
                  
                  tempfom = (neutronZcMean.at(currententry) - gammaZcMean.at(currententry))/(neutronZcStdDev.at(currententry) + gammaZcStdDev.at(currententry));
                  
                  std::cout << "For cfddelay : " << cfddelay << ", cfdattenuation: " << cfdattenuation << ", zcdiffrc: " << zcdiffrc << ", zcintrc1: " << zcintrc1 << ", zcintrc2: " << zcintrc2 << ", zc FOM = " << tempfom << std::endl;
                  
                  if(tempfom > bestfom) {
                     bestfom = tempfom;
                     bestzccfddelay = cfddelay;
                     bestzccfdattenuation = cfdattenuation;
                     bestzcdiffrc = zcdiffrc;
                     bestzcintrc1 = zcintrc1;
                     bestzcintrc2 = zcintrc2;
                  }
                  
                  ++currententry;
                  
               }
            }
         }
#endif
      }
   }

   std::cout << "Best neutron cfd standard deviation: " << bestneutroncfdstddev << ", at delay: " << bestneutroncfddelay << ", and attenuation: " << bestneutroncfdattenuation <<std::endl;
   std::cout << "Best gamma cfd standard deviation: " << bestgammacfdstddev << ", at delay: " << bestgammacfddelay << ", and attenuation: " << bestgammacfdattenuation <<std::endl;
#ifdef partialsumcfd
   std::cout << "Best figure of merit:" << bestfom << ", at zc delay: " <<  bestzcdelay << ", zc attenuation: " << bestzcattenuation << ", delay: " << bestzccfddelay << ", and attenuation: " << bestzccfdattenuation << std::endl;
#else
   std::cout << "Best figure of merit:" << bestfom << ", at cfd delay: " <<  bestzccfddelay << ", cfd attenuation: " << bestzccfdattenuation << ", bestzcdiffrc: " << bestzcdiffrc << ", and zcintrc1: " << bestzcintrc1  << ", and zcintrc2: " << bestzcintrc2 << std::endl;
#endif
   std::cout << argv[0] << " done after " << w.RealTime() << " seconds" << std::endl;
   
   return 0;
}

#endif

int checkDescantHit(TDescantHit* hit) {
   
   std::vector<Short_t>* waveform;
   std::vector<Short_t>::iterator maxelement;
   if(hit == 0) {
      return -1;
   }
   
   waveform = hit->GetWaveform();
   
   maxelement = std::max_element(waveform->begin(), waveform->end());
   
   if(std::distance(maxelement, waveform->end()) > 2 && maxelement[0] == maxelement[1] && maxelement[0] == maxelement[1]) {
      //Overflow condition satisified (max value is repeated 3 times in a row)
      return -1;
   }
   
   //TODO: Write an overflow check.
   
   //Scope Data
   //Max waveform value is between 500 and 2500, the max element does not occur in the first 100 or last 200 elements of the waveform.
   //The max element is within 10 samples of 420, the first and last elements are under 200.  The waveform does not go under -500
   //The cfd time is within 100 of 280.
   if(*maxelement > 2500 || *maxelement < 20 || std::distance(waveform->begin(),maxelement)<100 || std::distance(maxelement, waveform->end())<200 /*|| std::abs(std::distance(waveform->begin(), maxelement) - 420) > 10 */|| *std::min_element(waveform->begin(), waveform->end()) < -500 || waveform->front()>200 || waveform->back()>200 /*|| std::abs(hit->GetCfd()/256+280) < 100*/) {
      return -1;
   }
   
   return 0;
   
}




