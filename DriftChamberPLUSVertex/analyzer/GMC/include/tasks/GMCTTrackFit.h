#ifndef GMCTTrackFit_H
#define GMCTTrackFit_H

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// GMCTTrackFit                                                            //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#include "generated/GMCTTrackFit_Base.h"

class GMCTTrackFit : public GMCTTrackFit_Base
{

private:
   GMCTTrackFit(const GMCTTrackFit &c); // not implemented
   GMCTTrackFit &operator=(const GMCTTrackFit &c); // not implemented

public:
   GMCTTrackFit(const char *name = 0, const char *title = 0, int level = 0, const char *taskSuffix = 0, TFolder *histoFolder = 0)
   :GMCTTrackFit_Base(name,title,level,taskSuffix,histoFolder) {}
   virtual ~GMCTTrackFit() {}

protected:

   // Event Methods
   void Init();
   void BeginOfRun();
   void Event();
   void EndOfRun();
   void Terminate();
   Int_t MinuitFit();
   Int_t LeftRightMinuitFit();

   ClassDef(GMCTTrackFit,0)
};

#endif   // GMCTTrackFit_H
