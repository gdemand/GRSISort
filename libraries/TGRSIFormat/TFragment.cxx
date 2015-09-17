#include"TFragment.h"
#include"TChannel.h"
#include <iostream>

#include <TClass.h>

ClassImp(TFragment)

////////////////////////////////////////////////////////////////
//                                                            //
// TFragment                                                  //
//                                                            //
// This Class contains all of the information in an event     //
// fragment                                                   //
//                                                            //
////////////////////////////////////////////////////////////////

TFragment::TFragment() {
   // Default Constructor
#if MAJOR_ROOT_VERSION < 6
   Class()->IgnoreTObjectStreamer(kTRUE);
#endif
   Clear();
}

TFragment::TFragment(const TFragment& rhs, int hit) {
  //copy constructor that copies only the requested hit (if hit is in range 0 - Cfd.size()), if hit is negative, it act's as a normal copy constructor
  
  //first copy all "normal" data members
  MidasTimeStamp = rhs.MidasTimeStamp;
  MidasId = rhs.MidasId;
  TriggerId = rhs.TriggerId;
  FragmentId = rhs.FragmentId;
  TriggerBitPattern = rhs.TriggerBitPattern;

  NetworkPacketNumber = rhs.NetworkPacketNumber;

  ChannelNumber = rhs.ChannelNumber;
  ChannelAddress = rhs.ChannelAddress;

  TimeStampLow = rhs.TimeStampLow;
  TimeStampHigh = rhs.TimeStampHigh;

  PPG = rhs.PPG;
  DeadTime = rhs.DeadTime;
  NumberOfFilters = rhs.NumberOfFilters;
  NumberOfPileups = rhs.NumberOfPileups;
  DataType = rhs.DataType;
  DetectorType = rhs.DetectorType;
  ChannelId = rhs.ChannelId;

  KValue = rhs.KValue;

  wavebuffer = rhs.wavebuffer;

  if(hit < 0 || hit >= Cfd.size()) {
	  Cfd = rhs.Cfd;
	  Zc = rhs.Zc;
	  ccShort = rhs.ccShort;
	  ccLong = rhs.ccLong;
	  Led = rhs.Led;
	  Charge = rhs.Charge;
  } else {
	  Cfd.push_back(rhs.Cfd[hit]);
	  Zc.push_back(rhs.Zc[hit]);
	  ccShort.push_back(rhs.ccShort[hit]);
	  ccLong.push_back(rhs.ccLong[hit]);
	  Led.push_back(rhs.Led[hit]);
	  Charge.push_back(rhs.Charge[hit]);
  }

  NumberOfHits = Cfd.size();
  HitIndex = hit;
}

TFragment::~TFragment() {
  // Default destructor does nothing right now
  //Clear();
}

void TFragment::Clear(Option_t *opt) {  
   // Clears all fields of the TFragment
   MidasTimeStamp    = 0;
   MidasId           = 0;

   TriggerId         = 0;
   FragmentId        = 0;
   TriggerBitPattern = 0;

   NetworkPacketNumber = 0;

   ChannelAddress    = -1;
   Cfd.clear();//              = -1;
   Zc.clear();//              = -1;
   ccShort.clear();//
   ccLong.clear();//
   Led.clear();//               = -1;
   Charge.clear();//            = -1;
  
   //TimeStamp         = -1;
   TimeStampLow  = 0;
   TimeStampHigh = 0;


   PPG                = 0; 
   DeadTime           = 0; 
   NumberOfFilters    = 0; 
   NumberOfPileups    = 0; 
   DataType           = 0; 
   DetectorType       = 0; 
   ChannelId          = 0; 



   if(!wavebuffer.empty())
      wavebuffer.clear();

   if(!KValue.empty())     //->
   	KValue.clear();      //->

	NumberOfHits = -1;
	HitIndex = -1;
}

long TFragment::GetTimeStamp() const {
   long time = TimeStampHigh;
   time  = time << 28;
   time |= TimeStampLow & 0x0fffffff;
   return time;
}


double TFragment::GetTime() const {
   double dtime = (double)(GetTimeStamp())+ gRandom->Uniform();
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan )//|| Charge.size()<1)
      return dtime;
   return dtime - chan->GetTZero(GetEnergy());
}



double TFragment::GetTZero() const {
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan )//|| Charge.size()<1)
      return 0.000;
   return chan->GetTZero(GetEnergy());
}


long TFragment::GetCfdTime() {
   long ns = 0;
	if(DataType < 0 && index < Cfd.size()) {//TGriffin
	  //cfd replaces the lowest 26 bits of the timestamp
	  //so ns = CFD - 26 LSB of timestamp
	  ns = Cfd.at(index) - TimeStampLow & 0x03ffffff
	} else if(DataType == 2 && index < Cfd.size()) {//TSceptar
	  //cfd is a difference from the leading edge time stamp in 1/256th of a ns
	  //but also contains an offset between the 100Mhz timestamp clock and the 125MHz data clock
	  //so we add/subtract (?) the proper cfd divided by 256 to the offset
	  ns = (Cfd.at(index) >> 21) & 0xf - (Cfd.at(index) >> 8) & 0x1fff;
   }
   return 10*GetTimeStamp() + ns;  
}



long TFragment::GetTimeStamp_ns() {
   long ns = 0;
   if(DataType==2 && Cfd.size()>0) {
     ns = (Cfd.at(0) >> 21) & 0xf;
   }
   return 10*GetTimeStamp() + ns;  
}

Int_t TFragment::Get4GCfd(int i) { // return a 4G cfd in terms 
  if(Cfd.size()==0)                // of 1/256 ns since the trigger
     return -1;
  if(Cfd.size()<i)
     i = Cfd.size()-1;
  return  Cfd.at(i)&0x001fffff;

}


const char *TFragment::GetName() const {
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan)
      return "";
   return chan->GetChannelName();
}

/*
double TFragment::GetEnergy() const {
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan || Charge.size()<1)
      return 0.00;
   return chan->CalibrateENG((int)(Charge.at(0)));
}
*/

double TFragment::GetEnergy(int i) const {
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan || !(Charge.size()>i))
      return 0.00;
   if(chan->UseCalFileIntegration()) {
      //printf("I am here\n");
     return chan->CalibrateENG((int)(Charge.at(i)),0);  // this will use the integration value
                                                        // in the TChannel if it exists.
   }
   if(KValue.size()>i && KValue.at(i)>0)
     return chan->CalibrateENG((int)(Charge.at(i)),(int)KValue.at(i));
   return chan->CalibrateENG((int)(Charge.at(i)));
}

Float_t TFragment::GetCharge(int i) const {
   TChannel *chan = TChannel::GetChannel(ChannelAddress);
   if(!chan || !(Charge.size()>i))
      return 0.00;
   if(chan->UseCalFileIntegration()) {
      return ((Float_t)Charge.at(i)+gRandom->Uniform())/((Float_t)chan->GetIntegration());// this will use the integration value
   }                                                                                      // in the TChannel if it exists.
   if(KValue.size()>i && KValue.at(i)>0){
      return ((Float_t)Charge.at(i)+gRandom->Uniform())/((Float_t)KValue.at(i));// this will use the integration value
   }
   return ((Float_t)Charge.at(i)+gRandom->Uniform());// this will use no integration value
}

void TFragment::Print(Option_t *opt) const {
   //Prints out all fields of the TFragment


   TChannel *chan = TChannel::GetChannel(this->ChannelAddress);
   //printf("%s Event at	%i:\n", chan->GetDigitizerType().c_str(), MidasId);
   char buff[20];
   ctime(&MidasTimeStamp);
   struct tm * timeinfo = localtime(&MidasTimeStamp);
   strftime(buff,20,"%b %d %H:%M:%S",timeinfo);
   printf("MidasTimeStamp: %s\n",buff);
   printf("MidasId    	%i\n", MidasId);
   printf("TriggerId: 	%lu\n", TriggerId);
   printf("FragmentId:   %i\n", FragmentId);
   printf("TriggerBit:	0x%08x\n", TriggerBitPattern);
   printf("NetworkPacketNumber: %i\n", NetworkPacketNumber);
   if(chan)
	   printf("Channel: %i\tName: %s\n", chan->GetNumber(), chan->GetChannelName());
   printf("\tChannel Address: 0x%08x\n", ChannelAddress);
   printf("\tChannel Num:      %i\n", ChannelNumber);
   printf("\tCharge[%lu]	  ",Charge.size());   for(int x=0;x<Charge.size();x++){printf( "     0x%08x", Charge.at(x));} printf("\n");
   printf("\tCFD[%lu]		  ",Cfd.size());      for(int x=0;x<Cfd.size();x++)   {printf( "     0x%08x", Cfd.at(x));} printf("\n");
   printf("\tZC[%lu]		  ",Zc.size());      for(int x=0;x<Zc.size();x++)   {printf( "     0x%08x", Zc.at(x));} printf("\n");
   printf("\tLED[%lu]		  ",Led.size());      for(int x=0;x<Led.size();x++)   {printf( "     0x%08x", Led.at(x));} printf("\n");
   printf("\tTimeStamp High: 0x%08x\n", TimeStampHigh);
   printf("\tTimeStamp Low:    0x%08x\n", TimeStampLow);
   //unsigned short temptime = (TimeStampLow & 0x0000ffff) - ((Cfd >> 4) & 0x0000ffff);   //TimeStampLow&0x0000ffff; 
   //printf("\ttime from timestamp(to the nearest 10ns):    0x%04x\t%ins\n", temptime, temptime * 10);
   if (!wavebuffer.empty())
      printf("Has a wave form stored.\n");
   else
      printf("Does Not have a wave form stored.\n");

}




bool TFragment::IsDetector(const char * prefix, Option_t *opt) const {
   //Checks to see if the current fragment constains the same "prefix", for example "GRG"
   //The option determines whether the channel should be:
   // - C : The core of a segmented detector
   // - S : The segments of a segmented detector
   // - A : Low gain output of a detector
   // - B : High gain output of a detectora
   // If C or S are not given, default to C
   // If A or B are not given, default to A
   //Note that multiple options add to the output, so "CAB" would return the core with both high and low gain
   //One should eventually add N,P,T options as well.
   std::string pre = prefix;
   TString option = opt;
   std::string channame = this->GetName();
   if(channame.length()<9)
      return false;

   option.ToUpper();
   //Could also do everything below with MNEMONIC Struct. This limits the amount of string processing that needs to be done
   //Because it returns false after every potential failure while the mnemonic class sets all of the strings, and then checks
   //for conditions.
   if(channame.compare(0,pre.length(),pre)) {     //channame.BeginsWith(pre)){
      if(option.Length()<1) //no option.
         return true;
      if(channame.length()>8) {
        if(option.Contains("B") && (std::toupper(channame[9])==std::toupper('B')))
          return true;
        else if(option.Contains("A") && (std::toupper(channame[9])==std::toupper('A')))
          return true;
      }
      if(option.Contains("C") && !channame.compare(7,2,"00"))
        return true;
      else if(option.Contains("S") && channame.compare(7,2,"00"))
         return true;
   } else 
     return false;
   
   return false;
}
