#ifndef GMCTWaveformAnalysis_H
#define GMCTWaveformAnalysis_H

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// GMCTWaveformAnalysis                                                       //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#include "generated/GMCTWaveformAnalysis_Base.h"

class GMCTWaveformAnalysis : public GMCTWaveformAnalysis_Base
{

private:
   GMCTWaveformAnalysis(const GMCTWaveformAnalysis &c); // not implemented
   GMCTWaveformAnalysis &operator=(const GMCTWaveformAnalysis &c); // not implemented

public:
   GMCTWaveformAnalysis(const char *name = 0, const char *title = 0, int level = 0, const char *taskSuffix = 0, TFolder *histoFolder = 0)
   :GMCTWaveformAnalysis_Base(name,title,level,taskSuffix,histoFolder) {}
   virtual ~GMCTWaveformAnalysis() {}

protected:
   // Event Methods
   void Init();
   void BeginOfRun();
   void Event();
   void EndOfRun();
   void Terminate();


   ClassDef(GMCTWaveformAnalysis,0)
};

#endif   // GMCTWaveformAnalysis_H
