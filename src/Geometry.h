#ifndef Geometry_h
#define Geometry_h

//-----------------
//Geometry class
#include "vector"
#include "map"
#include <TGeoManager.h>
#include "TGeoTube.h"
#include "TGeoTrd2.h"
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TClass.h>
class Geometry {
public:
  Geometry(const char *geoFN, int verbose = 0);
  bool is2DStripMPGD(std::string coll_name);
  bool getGeoMats(int mpgd, bool &isFullyMSV);
  bool WorldToLocal(int idet, unsigned long cellID,
		    double *gpos,  // In mm
		    double *lpos); // In mm
  bool WorldToLocal(int idet, unsigned long cellID,
		    double *gpos, double *gmom, // In mm
		    double *lpos, double *lmom);
  bool LocalToWorld(int idet, unsigned long cellID,
		    double *lpos,  // In mm
		    double *gpos); // In mm
  static constexpr int N_MPGDs = 4;
  unsigned int MPGD2DStrips;
  std::vector<std::vector<const TGeoVolume*>> geoDetVols[4];
  std::vector<std::vector<TGeoHMatrix>>       geoDetMats[N_MPGDs];
  std::vector<std::vector<double>>            geoDetGaps[N_MPGDs];
private:
  int verboseLevel;
  // Tags allowing to find MPGD layers and modules
  static constexpr const char *tag1s[N_MPGDs] =
    {"Inner","Outer","Endcap",  "Endcap"};
  static constexpr const char *tag2s[N_MPGDs] =
    {0,      0,      "Backward","Forward"};
};
#endif // #ifdef Geometry_cxx
