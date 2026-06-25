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
    //double dZs[nSs] = {670.0,103.5,-528.5,-1095.0}; // mm
    //for (int section = 0; section<nSs; section++)
    //sectionDZs[0].push_back(dZs[section]);
    // <constant name="MMModuleWidth"                          value="46.0*cm"/>
    // Width is converted in angle
    /*
    // <constant name="MMOuterSector_R"                        value="57.7*cm"/>
    // <constant name="MMInnerSector_R"                        value="55.5*cm"/>
    */
    // <constant name="InnerMPGDBarrel_rmin"            value="555*mm"/>
    // <constant name="MMInnerSection_R"                       value="InnerMPGDBarrel_rmin + 0.5*cm"/>
    // <constant name="MMOuterSection_R"                       value="MMInnerSection_R + 2.2*cm"/>
    double hwidth = 230, rmins[2] = {555+5,577+5};
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
    //sectionDZs[1].push_back(880); sectionDZs[1].push_back(-830);
#else
    // Radius = Distance to REFERENCE subVolume (cellID>>28&0xf)=0x0
    // Obtained by MasterToLocal transforming (0,0,0) in MPGDTrackerDigi.
    radii[1].push_back(732.4650);
    // <constant name="MPGDOuterBarrelModule_zmin1"     value="1925*mm"/>
    // <constant name="MPGDOuterBarrelModule_zmin2"     value="1675*mm"/>
    // <constant name="MPGDOuterBarrelModule_PCB_offset"              value="110*mm"/>
    // <constant name="MPGDOuterBarrelFrame_width"            value="15*mm"/>
    ZAbscissae[1].push_back(-1925+110+15); ZAbscissae[1].push_back(1675-110-15);
    //sectionDZs[1].push_back(720); sectionDZs[1].push_back(-970);
#endif
    // <constant name="MPGDOuterBarrelModule_width"                   value="360*mm"/>
    hWidths[1].push_back(180-15);
    // <constant name="MPGDOuterBarrelModule_Inset_length"            value="MPGDOuterBarrelModule_length - MPGDOuterBarrelModule_PCB_offset"/>
    ZHLengths[1] = ((ZAbscissae[1][1]-ZAbscissae[1][0])/2-15/*InnerFrame*/)/2;
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
      // <constant name="ForwardMPGD_zmin"             value="1285.0*mm"/>
      // <constant name="ForwardMPGDMod_offset"        value="125.0*mm"/>
      ZAbscissae[3].push_back( 1285); ZAbscissae[3].push_back( 1285+125);
    }
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
bool recoEvents::parseGeometry()
{
  Geometry *geo = Geometry::Ptr();
  if (!geo) {
    printf("** parseGeometry: No access to geometry\n");
    return false;
  }

  const char *shapes[Geometry::N_MPGDs] = {"TGeoTubeSeg","TGeoBBox","TGeoTrd2","TGeoTrd2"};

  // LOOP ON MPGDs
  printf("+++++++++++++++++ Parsing %d MPGDs\n",Geometry::N_MPGDs);
  int mpgd; unsigned int error; for (mpgd = 0, error = 0; mpgd<Geometry::N_MPGDs; mpgd++) {
    if (!(0x1<<mpgd&processedDetectors)) continue;
    {
      vector<vector<const TGeoVolume*>> &geoVols = geo->geoDetVols[mpgd];
      vector<vector<TGeoHMatrix>>       &geoMats = geo->geoDetMats[mpgd];
      vector<vector<double>>            &geoGaps = geo->geoDetGaps[mpgd];
      using Module2StaveType = map<int,int>;
      Module2StaveType &module2StaveType = module2StaveTypes[mpgd];
      // Init mpgd
      using VolumeMap = map<const char*,pair<const TGeoVolume*,int>>;
      VolumeMap volumeMap;
      NLayers[mpgd] = geoVols.size(); if ((int)geoMats.size()!=NLayers[mpgd]) {
	printf("** parseGeometry: MPGD 0x%x Inconsistency: geoMats/geoVols differ in size: %d/%d\n",
	       0x1<<mpgd,(int)geoMats.size(),NLayers[mpgd]);
	error = 0x1; break;
      }
      const char *shape = shapes[mpgd];
      // Init mpgd
      ZExtrema[mpgd][0] = 1e6; ZExtrema[mpgd][1] = -1e6;
      NModules[mpgd] = 0;
      int modulesPerSection = 0, firstSection = 1;
      int prvNModules = -1;
      // ***** LOOP ON LAYERS,modules
      for (int layer = 0; layer<NLayers[mpgd]; layer++) {
	vector<const TGeoVolume*> gVols = geoVols[layer];
	vector<TGeoHMatrix> &gMats = geoMats[layer];
	vector<double> &gGaps = geoGaps[layer];
	int nLayerModules = gVols.size(); if ((int)gMats.size()!=nLayerModules) {
	  printf("** parseGeometry: MPGD 0x%x Inconsistency: gMats/gVols[%d] differ in size: %d/%d\n",
		 0x1<<mpgd,layer,(int)gMats.size(),nLayerModules);
	  error |= 0x2; break;
	}
	if (prvNModules>=0 && nLayerModules!=prvNModules) {
	  printf("** parseGeometry: Unforeseen configuration: Varying #modules per Layer: %d in layer %d != %d in layer %d\n",
		 prvNModules,layer-1,nLayerModules,layer);
	  error |= 0x10; break;
	}
	NModules[mpgd] = nLayerModules;
	for (int module = 0; module<nLayerModules; module++) {
	  const TGeoVolume *v = gVols[module]; const TGeoShape *s = v->GetShape();
	  if (strcmp(s->ClassName(),shape)) {
	    printf("** parseGeometry: MPGD 0x%x Inconsistency: Volume of Layer,Module %d,%d not a \"%s\"\n",
		   0x1<<mpgd,layer,module,shapes[mpgd]);
	    error = 0x4; break;
	  }
	  // ***** STAVETYPE (i.e. distinct volume shape)
	  //       SECTION
	  // - New staveType/section? Then map module2StaveType.
	  int newStaveType, staveType;
	  const char *sN = s->GetName();
	  VolumeMap::const_iterator iv = volumeMap.find(sN);
	  if (iv==volumeMap.end()) {
	    newStaveType = 1;
	    staveType = volumeMap.size(); volumeMap[sN] = {v,staveType};
	  }
	  else {
	    newStaveType = 0;
	    staveType = iv->second.second;
	  }
	  module2StaveType[module] = staveType;
	  // # of modules per Section
	  if (firstSection) modulesPerSection++;
	  if (newStaveType) firstSection = 0; // Finalise # of modules par section
	  //*****  mpgd SPECIFICS
	  const TGeoHMatrix &m = gMats[module]; double gap = gGaps[module];
	  if (mpgd==0) {
	    if (newStaveType) {
	      // CyMBaL: Parse TGeoTubeseg
	      const TGeoTubeSeg *tube = dynamic_cast<const TGeoTubeSeg*>(s);
	      if (!tube) { error |= 0x4+16*module; break; }
	      double startPhi = tube->GetPhi1(), endPhi = tube->GetPhi2();
	      // In "https://root.cern.ch/root/html534/guides/users-guide/Geometry.html"
	      // TGeoTubeSeg: "phi1 is converted to [0,360] and phi2 > phi1."
	      // => Convert it to [-pi,+pi].
	      startPhi *= TMath::Pi()/180; endPhi *= TMath::Pi()/180;
	      startPhi -= 2 * TMath::Pi(); endPhi -= 2 * TMath::Pi();
	      if (fabs(startPhi+endPhi)>1e-6) {
		printf("** Traverse: CyMBaL Module %d Inconsistency: |startPhi|(%.6f) != endPhi(%.6f)\n",
		       module,startPhi,endPhi);
		error = 0x8; 
	      }
	      HWidths[mpgd].push_back(-startPhi);
	      double rMin = tube->GetRmin()*10, rMax = tube->GetRmax()*10; // in mm
	      double R = (rMin+rMax)/2; Radii[mpgd].push_back(R);
	      double dZ = tube->GetDz()*10; /* in mm */ zHLengths[mpgd] = dZ;
	    }
	    double R = Radii[mpgd][staveType], hW = HWidths[mpgd][staveType];
	    // CyMBaL: Parse TGeoHMatrix
	    double dZ = zHLengths[mpgd];
	    double Z = m.GetTranslation()[2]; Z += Z>0 ? +dZ : -dZ;
	    if (Z<ZExtrema[mpgd][0]) ZExtrema[mpgd][0] = Z;
	    if (Z>ZExtrema[mpgd][1]) ZExtrema[mpgd][1] = Z;
	    /*
	    printf("CyMBal: L,M %d,%-2d: hW %.4f rad. R %.3f Z [%.1f,%.1f] dZ %.1f G %.2f mm Stave %d/%d MpS %d\n",
	           layer,module,hW,R,
		   ZExtrema[mpgd][0],ZExtrema[mpgd][1],dZ,gap,
		   staveType,(int)volumeMap.size(),modulesPerSection);
	    */
	  }
	  else if (mpgd==1) {
	    if (newStaveType) {
	      // Outer: Parse TGeoBBox
	      const TGeoBBox *box = dynamic_cast<const TGeoBBox*>(s);
	      if (!box) { error |= 0x4+16*module; break; }
	      double dX = box->GetDX()*10; HWidths[mpgd].push_back(dX);
	      double dY = box->GetDY()*10; zHLengths[mpgd] = dY;
	      double dZ = box->GetDZ()*10;
	    }
	    // Outer: Parse TGeoHMatrix
	    double dY = zHLengths[mpgd];
	    double Z = m.GetTranslation()[2]; Z += Z>0 ? +dY : -dY;
	    double X = m.GetTranslation()[0], Y = m.GetTranslation()[1];
	    double R = sqrt(X*X+Y*Y);
	    double hW = HWidths[mpgd][staveType];
	    if (newStaveType) Radii[mpgd].push_back(R);
	    if (Z<ZExtrema[mpgd][0]) ZExtrema[mpgd][0] = Z;
	    if (Z>ZExtrema[mpgd][1]) ZExtrema[mpgd][1] = Z;
	    /*
	    printf("Outer: L,M %d,%-2d: hW %.3f R %.3f Z [%.1f,%.1f] dZ %.1f G %.2f mm Stave %d/%d\n",
		   layer,module,hW,R,
		   ZExtrema[mpgd][0],ZExtrema[mpgd][1],dY,gap,staveType,(int)volumeMap.size());
	    */
	  }
	  else {
	    // Endcap: Parse TGeoTrd2
	    const TGeoTrd2 *trd2 = dynamic_cast<const TGeoTrd2*>(s);
	    if (!trd2) { error |= 0x4+16*module; break; }
	    double dX1 = trd2->GetDx1()*10, dX2 = trd2->GetDx2()*10;
	    //double dY1 = trd2->GetDy1()*10, dY2 = trd2->GetDy2()*10; already known: its the gap
	    double dZ = trd2->GetDZ()*10;	    
	    // Endcap: Parse TGeoTrd2
	    double X = m.GetTranslation()[0], Y = m.GetTranslation()[1];
	    double R = 2*dZ+sqrt(X*X+Y*Y);
	    double Z = m.GetTranslation()[2];
	    if (Z<ZExtrema[mpgd][0]) ZExtrema[mpgd][0] = Z;
	    if (Z>ZExtrema[mpgd][1]) ZExtrema[mpgd][1] = Z;
	    if (newStaveType) {
	      double rMin = R-dZ, rMax = R+dZ;
	      Radii[mpgd].push_back(rMin); Radii[mpgd].push_back(rMax);
	    }
	    /*
	    printf("Endcap: L,M %d,%-2d %.3f/%.3f %.3f R %.3f [%.3f,%.3f] Z [%.3f,%.3f] G %.2f mm Stave %d/%d <%s>\n",
		   layer,module,dX1,dX2,dZ,R,R-dZ,R+dZ,
		   ZExtrema[mpgd][0],ZExtrema[mpgd][1],gap,staveType,(int)volumeMap.size(),sN);
	    */
	  }
	}
	if (error&0x4) {
	  int module = error>>4;
	  printf("** parseGeometry: MPGD 0x%x Inconsistency: Volume of Layer,Module %d,%d can't be cast into a \"%s\"\n",
		 0x1<<mpgd,layer,module,shapes[mpgd]);
	}
      }
      if (error) return false;
    }
  } // End loop on mpgds
  // ***** COMPARISON
  double diff;
  for (mpgd = 0; mpgd<Geometry::N_MPGDs; mpgd++) {
    if (!(0x1<<mpgd&processedDetectors)) continue;
    printf("===+==+==== %d\n",mpgd);
    using Module2StaveType = map<int,int>;
    Module2StaveType &module2StaveType = module2StaveTypes[mpgd];
    int verbose = mpgd==5;
    for (int ud = 0; ud<2; ud++) {
      double ZAbscissa = ZAbscissae[mpgd][ud];
      double ZExtremum = ZExtrema  [mpgd][ud];
      diff = fabs(ZAbscissa-ZExtremum); if (diff>1.e-6 || verbose) {
	printf("%d ZExtremum %.6f ZAbscissa %.6f %.6f\n",
	       ud,ZExtremum,ZAbscissa,diff);
      }
    }
    // # of layers
    int NLayer = NLayers[mpgd],  nLayer = nLayers[mpgd];
    if (NLayer!=nLayer) {
      printf("#Layers(%d) != #layers(%d)\n",NLayer,nLayer);
    }
    // # of modules
    int NModule = NModules[mpgd],  nModule = nModules[mpgd];
    if (NModule!=nModule) {
      printf("#Modules(%d) != #modules(%d)\n",NModule,nModule);
    }
    // # of radii
    int nrs = radii[mpgd].size(), nRs = Radii[mpgd].size();
    if (nRs!=nrs) {
      printf("#Radii(%d) != #radii(%d)\n",nRs,nrs);
    }
    // radii
    for (int module = 0; module<(NModule<nModule?NModule:nModule); module++) {
      unsigned long cellID = (mpgd==2 || mpgd==3) ? (module<<10) : (module<<12);
      int Type = GetStaveType(mpgd,cellID), type = getStaveType(mpgd,cellID);
      double Radius = Radii[mpgd][Type];
      double radius = radii[mpgd][type];
      diff = fabs(Radius-radius); if (diff>1.e-3 || verbose) {
	printf("%d:%d/%d Radius %.6f radius %.6f %.6f\n",
	       module,Type,type,Radius,radius,diff);
	break;
      }
    }
    // # of hWidths
    int nHWs = HWidths[mpgd].size(), nhWs = hWidths[mpgd].size();
    if (nHWs!=nhWs) {
      printf("#HWidths(%d) != #hWidths(%d)\n",nHWs,nhWs);
    }
    if (mpgd>=2) continue;
    // hWidths
    for (int module = 0; module<(NModule<nModule?NModule:nModule); module++) {
      unsigned long cellID = (mpgd==2 || mpgd==3) ? (module<<10) : (module<<12);
      int Type = GetStaveType(mpgd,cellID), type = getStaveType(mpgd,cellID);
      double HWidth = HWidths[mpgd][Type];
      double hWidth = hWidths[mpgd][type];
      diff = fabs(HWidth-hWidth); if (diff>1.e-6 || verbose) {
	printf("%d:%d/%d HWidth %.6f hWidth %.6f %.6f\n",
	       module,Type,type,HWidth,hWidth,diff);
	break;
      }
    }
    double zHLength = zHLengths[mpgd], ZHLength = ZHLengths[mpgd];
    diff = fabs(zHLength-ZHLength); if (diff>1.e-6 || verbose) {
      printf("zHLength %.6f ZHLength %.6f %.6f\n",
	     zHLength,ZHLength,diff);
    }
  }
  return true;
}
