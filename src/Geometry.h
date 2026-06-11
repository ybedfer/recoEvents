#ifndef Geometry_h
#define Geometry_h

//-----------------
//Geometry class
#include "vector"
#include <TGeoManager.h>
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
  // ***** SINGLETON CLASS
  static inline Geometry* address;
  static Geometry* Ptr()
  {
    if (address) {
      return(address);
    }
    else {
      printf("** Geometry::Ptr(): Inconsistency: The object of Geometry class is not yet created or already destructed\n");
      return 0;
    }
  };
  static constexpr int N_MPGDs = 4;
  unsigned int MPGD2DStrips;
private:
  std::vector<std::vector<TGeoHMatrix>> geoDetMats[N_MPGDs];
  int verboseLevel;
  // Tags allowing to find MPGD layers and modules
  static constexpr const char *tag1s[N_MPGDs] =
    {"Inner","Outer","Endcap",  "Endcap"};
  static constexpr const char *tag2s[N_MPGDs] =
    {0,      0,      "Backward","Forward"};
};
#endif // #ifdef Geometry_cxx
