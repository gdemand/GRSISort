#ifndef TGRSIOPTIONS_H
#define TGRSIOPTIONS_H

/** \addtogroup Sorting
 *  @{
 */

#include <cstdio>
#include <string>

#include "TObject.h"

namespace TGRSIOptions {
    namespace priv{
      extern std::string fHostName;
      extern std::string fExptName;
      
      extern std::vector<std::string> fInputRootFile;
      extern std::vector<std::string> fInputMidasFile;
      extern std::vector<std::string> fInputCalFile;
      extern std::vector<std::string> fInputOdbFile;
    	extern std::vector<std::string> fExternalRunInfo;
    	extern std::vector<std::string> fMacroFile;

      extern bool fCloseAfterSort;
      extern bool fLogErrors;
      extern bool fUseMidFileOdb;
      extern bool fMakeAnalysisTree;
      extern bool fProgressDialog;
      extern bool fWorkHarder;
      extern bool fReadingMaterial;
      extern bool fIgnoreFileOdb;
		extern bool fIgnoreScaler;
		extern bool fIgnoreEpics;
      extern bool fWriteBadFrags;
		extern bool fWriteDiagnostics;
      }
      std::string GetHostName();
      std::string GetExptName();
     
      std::vector<std::string> GetInputRoot();  
      std::vector<std::string> GetInputMidas(); 
      std::vector<std::string> GetInputCal();   
      std::vector<std::string> GetInputOdb();   
      std::vector<std::string> GetMacroFile();   

		const char *GetXMLODBFile(int runNumber=0,int subRunNumber=-1);
      const char *GetCalFile(int runNumber=0,int subRunNumber=-1);

      void AddExternalRunInfo(std::string);
      void SetExternalRunInfo();
      bool ExternalRunInfo();

      void SetCloseAfterSort(bool flag=true); 
      bool CloseAfterSort();                  

      void SetIgnoreFileOdb(bool flag=true);
      bool IgnoreFileOdb();

      void SetIgnoreScaler(bool flag=true);
      bool IgnoreScaler();

      void SetIgnoreEpics(bool flag=true);
      bool IgnoreEpics();
      
      void SetIgnoreSCLR(bool flag=true);
      bool IgnoreSCLR();

      void SetLogErrors(bool flag=true);      
      bool LogErrors();			
	
      void SetProgressDialog(bool flag=true); 
      bool ProgressDialog();                  
	
      void SetUseMidFileOdb(bool flag=true);  
      bool UseMidFileOdb();                   
      
      void SetMakeAnalysisTree(bool flag=true);
      bool MakeAnalysisTree();                

      void SetWorkHarder(bool flag=true);
      bool WorkHarder();      

      void SetReadingMaterial(bool flag=true);
      bool ReadingMaterial();

      void SetWriteBadFrags(bool flag=true);
      bool WriteBadFrags();

      void SetWriteDiagnostics(bool flag=true); 
      bool WriteDiagnostics();                  

      void SetHostName(std::string &host);
      void SetExptName(std::string &expt); 
      
      void AddInputRootFile(std::string &input);  
      void AddInputMidasFile(std::string &input); 
      void AddInputCalFile(std::string &input);   
      void AddInputOdbFile(std::string &input);   
      void AddMacroFile(std::string &input);   
}
/*! @} */
#endif
