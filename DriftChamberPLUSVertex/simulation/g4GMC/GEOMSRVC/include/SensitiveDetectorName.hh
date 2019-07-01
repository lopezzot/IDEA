/*
 * SensitiveDetectorName.hh
 *
 *  Created on: Apr 3, 2017
 *      Author: tassiell
 */

#ifndef SENSITIVEDETECTORNAME_HH_
#define SENSITIVEDETECTORNAME_HH_

//namespace cdch {

class SensitiveDetectorName {
//static const char* name[4] = { "DCHtrackerSD", "DCHSWiresSD",  "DCHFWiresSD", "DCHWallsSD" };

public:
  enum SDtype { DCHtracker=0, DCHSWires,  DCHFWires, DCHWalls, SVXtracker, PSHWtracker };

  static char const * TrackerGas(){
    return name[DCHtracker];
  }

  static char const * TrackerSWires(){
    return name[DCHSWires];
  }

  static char const * TrackerFWires(){
    return name[DCHFWires];
  }

  static char const * TrackerWalls(){
    return name[DCHWalls];
  }

  static char const * SVXTrackerRO(){
    return name[SVXtracker];
  }

  static char const * PSHWTrackerRO(){
    return name[PSHWtracker];
  }

private:
  static const char* name[6];

};

//const char* SensitiveDetectorName::name[4] = { "DCHtrackerSD", "DCHSWiresSD",  "DCHFWiresSD", "DCHWallsSD" };

//} // namespace cdch

#endif /* SENSITIVEDETECTORNAME_HH_ */
