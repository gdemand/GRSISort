
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>

#include "TTree.h"
#include "TSystem.h"

#include "TFile.h"
#include "TTree.h"
#include "TDescant.h"
#include "TGriffin.h"
#include "TPaces.h"
#include "TSceptar.h"

#include "/Users/gdemand/DESCANT/Programs/LecroyScope2Root/Oscilloscope.hh"

int CreateAndAddDescantHitFromScope(TDescant &Descant, int eventnumber, std::vector<double> *time, std::vector<double> *voltage, double minvalue);

int main(int argc, char* argv[]) {
   
   if(argc < 2) {
      std::cout << "Usage is ConvertScopeData inputfile.root" << std::endl;
      return 0;
   }
   std::cout << argc << std::endl;
   
   std::string inputfilename(argv[1]);
   std::string outputfilename(argv[1]);

   if(outputfilename.find_last_of("/") != std::string::npos) {
      outputfilename.insert(outputfilename.find_last_of("/")+1,"analysis_");
   } else {
      outputfilename.insert(0,"analysis_");
   }

   if(outputfilename.find_last_of("/") != std::string::npos) {
      outputfilename.erase(0, outputfilename.find_last_of("/")+1);
   }

   std::cout << "outputfilename is now:" << outputfilename<< std::endl;
   
   gSystem->Load("/Users/gdemand/lib/libLecroyScope2Root.so");
   
   TFile inputfile(inputfilename.c_str());
   
   TFile f(outputfilename.c_str(), "recreate");
   
   TDescant Descant;
   TPaces Paces;
   TSceptar Sceptar;
   TGriffin Griffin;
   
   TChannel *channel = 0;
   
   int eventnumber = 0;
   
   std::string inputdirectory(argv[1]);
   
   std::vector<double> *inputVoltage;
   std::vector<double> *inputTime;
   
   TTree* inputTree = (TTree*) inputfile.Get("tree");

   TTree AnalysisTree("AnalysisTree", "AnalysisTree");

   
   inputTree->SetBranchAddress("time", &inputTime);
   inputTree->SetBranchAddress("voltage", &inputVoltage);

   AnalysisTree.Branch("TDescant","TDescant",&Descant,128000,99);
   AnalysisTree.Branch("TPaces","TPaces",&Paces,128000,99);
   AnalysisTree.Branch("TSceptar","TSceptar",&Sceptar,128000,99);
   AnalysisTree.Branch("TGriffin","TGriffin",&Griffin,128000,99);

   channel = new TChannel("");
   channel->SetChannelName("DSC01XN00X");
   channel->SetAddress(0);
   channel->SetUseCalFileIntegration(false);
   TChannel::AddChannel(channel);
   channel->WriteToRoot();

   LecroyHeader* scopeHeader = (LecroyHeader*) inputfile.Get("header");
   
   std::cout << "Min value: " << scopeHeader->GetMinVoltage() << std::endl;
   double minvalue = scopeHeader->GetMinVoltage();
   
   for(long long entry = 1; entry < inputTree->GetEntries(); ++entry) {
      inputTree->GetEntry(entry);
      
//      std::cout << "Time: " << inputTime->at(1000) << ", voltage: " << inputVoltage->at(1000) <<std::endl;
      Descant.Clear();
      if(CreateAndAddDescantHitFromScope(Descant, eventnumber, inputTime, inputVoltage, minvalue) == 0) {
         AnalysisTree.Fill();
      }
      //return 0;
      
   }
   
   AnalysisTree.Write();
   
   return 0;
   
}


int CreateAndAddDescantHitFromScope(TDescant &Descant, int eventnumber, std::vector<double> *time, std::vector<double> *voltage, double minvalue)	{
   //Builds the GRIFFIN Hits from the "data" structure. Basically, loops through the data for and event and sets observables.
   //This is done for both GRIFFIN and it's suppressors.
   
   std::string waveform_string;
   
   std::vector<Short_t> waveform;
   
   waveform.clear();
   for(size_t i = 0; i < voltage->size(); i++) {
      if((*voltage)[i] < minvalue*0.99) {
         return 1;
      }
      waveform.push_back(round((*voltage)[i]*(-200)));
   }
   
   TDescantHit dethit;
   
   dethit.SetAddress(0);//gdata->GetDetAddress(i));
   dethit.SetCharge(0);//gdata->GetDetCharge(i));
   
   dethit.SetTime(0);//gdata->GetDetTime(i)); Should likely be SetTimeStamp
   dethit.SetCfd(0);//gdata->GetDetCFD(i));
   dethit.SetZc(0);//gdata->GetDetZc(i));
   dethit.SetCcShort(0);//gdata->GetDetCcShort(i));
   dethit.SetCcLong(0);//gdata->GetDetCcLong(i));
   
   dethit.SetWaveform(waveform);//gdata->GetDetWave(i));

   dethit.fTimes = *time;
   
   if(dethit.GetWaveform()->size() > 0) {
      dethit.AnalyzeScopeWaveform();
   }
   
   Descant.AddHit(&dethit);

   
   return 0;
   
}