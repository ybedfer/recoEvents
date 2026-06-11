#include "recoEvents.h"
#include <TVector3.h>
#include <TRotation.h>

//#define USE_OLDER_GEOM // From before commit #896b85e80 by Matt

// ********** GEOMETRY
void recoEvents::initGeometry(int idet, bool hasStrips)
{
  radiatorThicknesses[idet] = 0;
  if (isMPGD(idet)) {
    //<constant name="BackwardMPGDDriftGap_thickness"           value="3.0*mm" />
    volumeThicknesses[idet] = 3;
    if (hasStrips) {
      const double radThickness = (3-3*.01)/2; // 3mm thickness-3*HELPER SUBVOLUMES
      radiatorThicknesses[idet] = radThickness;
    }
  }
  else if (idet==4) {
    // <constant name="VertexBarrelMod_thickness" value="0.2*mm" />
    volumeThicknesses[idet] = .04;
  }
  else if (idet==5) {
    // <constant name="SiVertexSensor_thickness" value="40*um" />
    volumeThicknesses[idet] = .04;
  }
  if (idet==0) {                         // ***** CyMBaL = 2 radii, 4 sections
    const int nSs = 4; nSections[0] = nSs;
#ifdef USE_OLDER_GEOM
    radii[0].push_back(556.755); radii[0].push_back(578.755);
#else
    radii[0].push_back(561.755); radii[0].push_back(583.755);
#endif
    // <constant name="InnerMPGDBarrel_zmin"            value="1025*mm"/> <comment> negative z </comment>
    // <constant name="InnerMPGDBarrel_zmax"            value="1450*mm"/> <comment> positive z </comment>
    // <constant name="MMOutwardFrameWidth"                    value="5.0*cm"/>
    ZAbscissae[0].push_back(-1025+50); ZAbscissae[0].push_back(1450-50);
    double dZs[nSs] = {670.0,103.5,-528.5,-1095.0}; // mm
    for (int section = 0; section<nSs; section++)
      sectionDZs[0].push_back(dZs[section]);
    // <constant name="MMModuleWidth"                          value="46.0*cm"/>
    // Width is converted in angle
    // <constant name="MMOuterSector_R"                        value="57.7*cm"/>
    // <constant name="MMInnerSector_R"                        value="55.5*cm"/>
    double hwidth = 230, rmins[2] = {555,577};
    for (int io = 0; io<2; io++) hWidths[0].push_back(hwidth/rmins[io]);
    // <constant name="MMModuleLength"                         value="61.0*cm"/>
    ZHLengths[0] = 305;
    //<constant name="MMnStripsPhi"    value = "512" /> 
    //<constant name="MMnStripsZ"      value = "512" />
    int nStripsPhi = 512; pitches[0].push_back(2*hwidth/nStripsPhi);
    int nStripsZ =   512; pitches[0].push_back(2*ZHLengths[0]/nStripsZ);
    nChannels[0].push_back(nStripsPhi); nChannels[0].push_back(nStripsZ);
  }
  else if (idet==1) {                    // ***** OUTER
#ifdef USE_OLDER_GEOM
    radii[1].push_back(737.4650);
    // <constant name="MPGDOuterBarrelModule_zmin1"     value="1795*mm"/>
    // <constant name="MPGDOuterBarrelModule_zmin2"     value="1845*mm"/>
    // <constant name="MPGDOuterBarrelModule_PCB_offset"              value="110*mm"/>
    // <constant name="MPGDOuterBarrelFrame_width"            value="15*mm"/>
    ZAbscissae[1].push_back(-1795+110+15); ZAbscissae[1].push_back(1845-110-15);
    sectionDZs[1].push_back(880); sectionDZs[1].push_back(-830);
#else
    // Radius = Distance to REFERENCE subVolume (cellID>>28&0xf)=0x0
    // Obtained by MasterToLocal transforming (0,0,0) in MPGDTrackerDigi.
    radii[1].push_back(732.4650);
    // <constant name="MPGDOuterBarrelModule_zmin1"     value="1925*mm"/>
    // <constant name="MPGDOuterBarrelModule_zmin2"     value="1675*mm"/>
    // <constant name="MPGDOuterBarrelModule_PCB_offset"              value="110*mm"/>
    // <constant name="MPGDOuterBarrelFrame_width"            value="15*mm"/>
    ZAbscissae[1].push_back(-1925+110+15); ZAbscissae[1].push_back(1675-110-15);
    sectionDZs[1].push_back(720); sectionDZs[1].push_back(-970);
#endif
    // <constant name="MPGDOuterBarrelModule_width"                   value="360*mm"/>
    hWidths[1].push_back(180-15);
    ZHLengths[1] = 840;
    //<constant name="MPGDOuterBarrelPitch"                        value = "800*um" />
    double outerPitch = 800; pitches[1].push_back(outerPitch/1000);
    // <constant name="MPGDOuterBarrelnStrips"                      value = "1792" />
    nChannels[1].push_back(1792); nChannels[1].push_back(1792);
  }
  else if (idet==2 || idet==3) {         // ***** Endcaps
    // <constant name="BackwardMPGDMod1_rmin"         value="70.0*mm" />
    // <constant name="BackwardMPGDMod1_rmax"         value="400.0*mm" />
    radii[idet].push_back(70); radii[idet].push_back(400);
    if (idet==2) {
      // <constant name="BackwardMPGD_zmin"             value="1075.0*mm"/>
      // <constant name="BackwardMPGDMod_offset"        value="125.0*mm"/>
      ZAbscissae[2].push_back(-1075); ZAbscissae[2].push_back(-1075-125);
    }
    else {
      // <constant name="ForwardMPGD_zmin"             value="1500.0*mm"/>
      // <constant name="ForwardMPGDMod_offset"        value="125.0*mm"/>
      ZAbscissae[3].push_back( 1500); ZAbscissae[3].push_back( 1500+125);
    }
    double f_rmin = 70, f_rmax = 400;
    radii[3].push_back(f_rmin); radii[3].push_back(f_rmax);
  }
  else if (idet==4) {
    // <constant name="VertexBarrel_rmin"               value="3.6*cm"/>
    // <constant name="VertexBarrel_rmax"               value="12.6*cm"/>
    radii[4].push_back(36); radii[4].push_back(126);
    // <constant name="VertexBarrel_length"             value="26.6*cm"/>
    // <constant name="RSU_length"  value="21.666*mm" />
    // <constant name="VertexBarrelMod_length"   value="RSU_length/2" />
    // <constant name="VertexBarrelLayer_nz"     value="12*2" />
    // <constant name="VertexBarrelLayer_length" value="VertexBarrelMod_length*VertexBarrelLayer_nz" />
    double hLength = 65/3.*12/2;
    ZAbscissae[4].push_back(-hLength); ZAbscissae[4].push_back(+hLength);
  }
  else if (idet==5) {
    // <constant name="SiBarrelMod1_rc" value="26.5*cm" /> <!-- 26.5 cm is the average radius for inner sub-layer of 262mm and outer of 267mm-->
    // <constant name="SiBarrelMod2_rc" value="42*cm" /> <!-- 42 cm, inner/outer 417mm, 423mm -->
    radii[5].push_back(265); radii[5].push_back(420);
    // <constant name="SiBarrelMod2_length" value="84*cm - 4.7*cm" /> <!--UPDATED from 84*cm to 84*cm - 4.7*cm = 79.3cm-->
    ZAbscissae[5].push_back(-396.5); ZAbscissae[5].push_back(+396.5);    
  }
  // ***** GAINS, THRESHOLDS
  if (isMPGD(idet)) {
    if (idet==0 || idet==1) {
      //digi_cfg.gain                = 10000;
      gains[idet] = 10000;
    }
    else
      gains[idet] = 1;
    // Thresholds on eDep
    //.threshold      = 100 * dd4hep::eV, in "MPGD.cc"
    eDThresholds[idet] = .1;
    //digi_cfg.stripResolutions[0] = digi_cfg.stripResolutions[1] = 150 * dd4hep::um;
    resolutions[idet] = 150; // in µm
  }
  else {
    gains[idet] = 1;
    //.threshold = 0.54 * dd4hep::keV, in "BVTX.cc","BTRK.cc"
    eDThresholds[idet] = .54;
    resolutions[idet] = 0; // Not yet set...
  }
}
