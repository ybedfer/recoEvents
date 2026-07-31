#include "recoEvents.h"
#include <TVector3.h>
#include <TRotation.h>

// ********** GEOMETRY
bool recoEvents::extraGeometry(int idet, bool hasStrips)
{
  // Set geometry of non MPGDs (which latter get their geometry from Geometry class
  // Set extra parameters (#channels and pitch, resolution, gain and pitch)

  if (!isMPGD(idet)) {
    // ***** #LAYERS, #MODULES
    int N_MPGDs = Geometry::N_MPGDs;
    if (N_DETs-N_MPGDs!=2) {
      printf("** extraGeometry: Inconsistency: N_DETs(=%d) != N_MPGDs(=%d)+2\n",
	     N_DETs,N_MPGDs);
      return false;
    }
    const int NLayers[2]  = {   3, 2};
    const int NModules[2] = {1920,80};
    int jdet = idet-N_MPGDs;
    nLayers[idet] = NLayers[jdet]; nModules[idet] = NModules[jdet];
  }
  // ***** THICKNESSES
  radiatorThicknesses[idet] = 0;
  if (isMPGD(idet)) {
    if (hasStrips) {
      const double radThickness = (3-3*.01)/2; // 3mm thickness-3*HELPER SUBVOLUMES
      radiatorThicknesses[idet] = radThickness;
    }
  }
  else if (idet==4) {
    // <constant name="VertexBarrelMod_thickness" value="0.2*mm" />
    volumeThicknesses[idet] = .20;
  }
  else if (idet==5) {
    // <constant name="SiVertexSensor_thickness" value="40*um" />
    volumeThicknesses[idet] = .04;
  }
  // ***** #CHANNELS, PITCH
  if (idet==0) {
    if (isCyMBaL_8S) {
      //<constant name="MMnStripsPhi"    value = "512" /> 
      //<constant name="MMnStripsZ"      value = "512" />
      int nStripsPhi = 512; pitches[0].push_back(2*hWidths[0][0]/nStripsPhi);
      int nStripsZ =   512; pitches[0].push_back(2*ZHLengths[0]/nStripsZ);
      nChannels[0].push_back(nStripsPhi); nChannels[0].push_back(nStripsZ);
    }
    else {
      //<constant name="MMnStripsX"        value="192"/>
      //<constant name="MMPitchX"          value="1.766*mm"/>
      //<constant name="MMnStripsY"        value="576"/>
      //<constant name="MMPitchY"          value="1.021"/>
      int nStripsPhi = 192, nStripsZ = 576;
      nChannels[0].push_back(192); nChannels[0].push_back(576);
      pitches[0].push_back(1.766); pitches[0].push_back(1.021);
    }
  }
  else if (idet==1) {                    // ***** OUTER
    //<constant name="MPGDOuterBarrelPitch"                        value = "800*um" />
    double outerPitch = 800; pitches[1].push_back(outerPitch/1000);
    //<constant name="MPGDOuterBarrelnStrips"                      value = "1792" />
    nChannels[1].push_back(1792); nChannels[1].push_back(1792);
  }
  else if (idet==2 || idet==3) {         // ***** Endcaps
    // No #channels, not pitch
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
    ZAbscissae[4][0] = -hLength; ZAbscissae[4][1] = +hLength;
  }
  else if (idet==5) {
    // <constant name="SiBarrelMod1_rc" value="26.5*cm" /> <!-- 26.5 cm is the average radius for inner sub-layer of 262mm and outer of 267mm-->
    // <constant name="SiBarrelMod2_rc" value="42*cm" /> <!-- 42 cm, inner/outer 417mm, 423mm -->
    radii[5].push_back(265); radii[5].push_back(420);
    // <constant name="SiBarrelMod2_length" value="84*cm - 4.7*cm" /> <!--UPDATED from 84*cm to 84*cm - 4.7*cm = 79.3cm-->
    ZAbscissae[5][0] = -396.5;   ZAbscissae[5][0] = +396.5;
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
    if (idet==0) {
      if (isCyMBaL_8S) {
	//digi_cfg.stripResolutions[0] = digi_cfg.stripResolutions[1] = 150 * dd4hep::um;
	resolutions[0].push_back(150); resolutions[0].push_back(150); // in µm
      }
      else {
	//<constant name="MMumResolutionPhi"     value="260"/>
	//<constant name="MMumResolutionZ"       value="150"/>
	resolutions[0].push_back(260); resolutions[0].push_back(150); // in µm
      }
    }
    else {
      resolutions[idet].push_back(150); // in µm
    }
  }
  else {
    gains[idet] = 1;
    //.threshold = 0.54 * dd4hep::keV, in "BVTX.cc","BTRK.cc"
    eDThresholds[idet] = .54;
    resolutions[idet].push_back(0); // Not yet set...
  }
  return true;
}

bool recoEvents::parseGeometry()
{
  const char *shapes[Geometry::N_MPGDs] = {"TGeoTubeSeg","TGeoBBox","TGeoTrd2","TGeoTrd2"};

  // LOOP ON MPGDs
  printf("+++++++++++++++++ Parsing %d MPGDs\n",Geometry::N_MPGDs);
  int mpgd; unsigned int error; for (mpgd = 0, error = 0; mpgd<Geometry::N_MPGDs; mpgd++) {
    if (!(0x1<<mpgd&processedDetectors)) continue;
    vector<vector<const TGeoVolume*>> &geoVols = geo->geoDetVols[mpgd];
    vector<vector<TGeoHMatrix>>       &geoMats = geo->geoDetMats[mpgd];
    vector<vector<double>>            &geoGaps = geo->geoDetGaps[mpgd];
    using Module2StaveType = map<int,int>;
    Module2StaveType &module2StaveType = module2StaveTypes[mpgd];
    // Init mpgd
    using VolumeMap = map<const char*,pair<const TGeoVolume*,int>>;
    VolumeMap volumeMap;
    nLayers[mpgd] = geoVols.size(); if ((int)geoMats.size()!=nLayers[mpgd]) {
      printf("** parseGeometry: MPGD 0x%x Inconsistency: geoMats/geoVols differ in size: %d/%d\n",
	     0x1<<mpgd,(int)geoMats.size(),nLayers[mpgd]);
      error = 0x1; break;
    }
    const char *shape = shapes[mpgd];
    // Init mpgd
    ZAbscissae[mpgd][0] = 1e6; ZAbscissae[mpgd][1] = -1e6;
    nModules[mpgd] = 0;
    volumeThicknesses[mpgd] = 0;
    int modulesPerSection = 0, firstSection = 1;
    int prvnModules = -1;
    // ***** LOOP ON LAYERS,modules
    for (int layer = 0; layer<nLayers[mpgd]; layer++) {
      vector<const TGeoVolume*> gVols = geoVols[layer];
      vector<TGeoHMatrix> &gMats = geoMats[layer];
      vector<double> &gGaps = geoGaps[layer];
      int nLayerModules = gVols.size(); if ((int)gMats.size()!=nLayerModules) {
	printf("** parseGeometry: MPGD 0x%x Inconsistency: gMats/gVols[%d] differ in size: %d/%d\n",
	       0x1<<mpgd,layer,(int)gMats.size(),nLayerModules);
	error |= 0x2; break;
      }
      if (prvnModules>=0 && nLayerModules!=prvnModules) {
	printf("** parseGeometry: Unforeseen configuration: Varying #modules per Layer: %d in layer %d != %d in layer %d\n",
	       prvnModules,layer-1,nLayerModules,layer);
	error |= 0x10; break;
      }
      nModules[mpgd] = nLayerModules;
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
	double prvThickness = volumeThicknesses[mpgd];
	if (prvThickness && gap!=prvThickness) {
	  printf("** parseGeometry: MPGD %d Module %d Inconsistency: Gap(=%.2f) != previous(%.2f)\n",
		 mpgd,module,gap,prvThickness);
	  error |= 0x20;
	}
	volumeThicknesses[mpgd] = gap;
	if (mpgd==0) {
	  if (newStaveType) {
	    // CyMBaL: Parse TGeoTubeseg
	    const TGeoTubeSeg *tube = dynamic_cast<const TGeoTubeSeg*>(s);
	    if (!tube) { error |= 0x4+256*module; break; }
	    double startPhi = tube->GetPhi1(), endPhi = tube->GetPhi2();
	    // In "https://root.cern.ch/root/html534/guides/users-guide/Geometry.html"
	    // TGeoTubeSeg: "phi1 is converted to [0,360] and phi2 > phi1."
	    // => Convert it to [-pi,+pi].
	    startPhi *= TMath::Pi()/180; endPhi *= TMath::Pi()/180;
	    startPhi -= 2 * TMath::Pi(); endPhi -= 2 * TMath::Pi();
	    if (fabs(startPhi+endPhi)>1e-6) {
	      printf("** parseGeometry: CyMBaL Module %d Inconsistency: |startPhi|(%.6f) != endPhi(%.6f)\n",
		     module,startPhi,endPhi);
	      error = 0x8;
	    }
	    hWidths[mpgd].push_back(-startPhi);
	    double rMin = tube->GetRmin()*10, rMax = tube->GetRmax()*10; // in mm
	    double R = (rMin+rMax)/2; radii[mpgd].push_back(R);
	    double dZ = tube->GetDz()*10; /* in mm */ ZHLengths[mpgd] = dZ;
	  }
	  double R = radii[mpgd][staveType], hW = hWidths[mpgd][staveType];
	  // CyMBaL: Parse TGeoHMatrix
	  double dZ = ZHLengths[mpgd];
	  double Z = m.GetTranslation()[2]; Z += Z>0 ? +dZ : -dZ;
	  if (Z<ZAbscissae[mpgd][0]) ZAbscissae[mpgd][0] = Z;
	  if (Z>ZAbscissae[mpgd][1]) ZAbscissae[mpgd][1] = Z;
	  /*
	    printf("CyMBal: L,M %d,%-2d: hW %.4f rad. R %.3f Z [%.1f,%.1f] dZ %.1f G %.2f mm Stave %d/%d MpS %d\n",
	    layer,module,hW,R,
	    ZAbscissae[mpgd][0],ZAbscissae[mpgd][1],dZ,gap,
	    staveType,(int)volumeMap.size(),modulesPerSection);
	  */
	}
	else if (mpgd==1) {
	  if (newStaveType) {
	    // Outer: Parse TGeoBBox
	    const TGeoBBox *box = dynamic_cast<const TGeoBBox*>(s);
	    if (!box) { error |= 0x4+256*module; break; }
	    double dX = box->GetDX()*10; hWidths[mpgd].push_back(dX);
	    double dY = box->GetDY()*10; ZHLengths[mpgd] = dY;
	    double dZ = box->GetDZ()*10;
	  }
	  // Outer: Parse TGeoHMatrix
	  double dY = ZHLengths[mpgd];
	  double Z = m.GetTranslation()[2]; Z += Z>0 ? +dY : -dY;
	  double X = m.GetTranslation()[0], Y = m.GetTranslation()[1];
	  double R = sqrt(X*X+Y*Y);
	  double hW = hWidths[mpgd][staveType];
	  if (newStaveType) radii[mpgd].push_back(R);
	  if (Z<ZAbscissae[mpgd][0]) ZAbscissae[mpgd][0] = Z;
	  if (Z>ZAbscissae[mpgd][1]) ZAbscissae[mpgd][1] = Z;
	  /*
	    printf("Outer: L,M %d,%-2d: hW %.3f R %.3f Z [%.1f,%.1f] dZ %.1f G %.2f mm Stave %d/%d\n",
	    layer,module,hW,R,
	    ZAbscissae[mpgd][0],ZAbscissae[mpgd][1],dY,gap,staveType,(int)volumeMap.size());
	  */
	}
	else {
	  // Endcap: Parse TGeoTrd2
	  const TGeoTrd2 *trd2 = dynamic_cast<const TGeoTrd2*>(s);
	  if (!trd2) { error |= 0x4+256*module; break; }
	  double dX1 = trd2->GetDx1()*10, dX2 = trd2->GetDx2()*10;
	  //double dY1 = trd2->GetDy1()*10, dY2 = trd2->GetDy2()*10; already known: its the gap
	  double dZ = trd2->GetDZ()*10;	    
	  // Endcap: Parse TGeoTrd2
	  double X = m.GetTranslation()[0], Y = m.GetTranslation()[1];
	  double R = sqrt(X*X+Y*Y);
	  double Z = m.GetTranslation()[2];
	  if (Z<ZAbscissae[mpgd][0]) ZAbscissae[mpgd][0] = Z;
	  if (Z>ZAbscissae[mpgd][1]) ZAbscissae[mpgd][1] = Z;
	  if (newStaveType) {
	    double rMin = R-dZ, rMax = R+dZ;
	    radii[mpgd].push_back(rMin); radii[mpgd].push_back(rMax);
	  }
	  /*
	    printf("Endcap: L,M %d,%-2d %.3f/%.3f %.3f R %.3f [%.3f,%.3f] Z [%.3f,%.3f] G %.2f mm Stave %d/%d <%s>\n",
	    layer,module,dX1,dX2,dZ,R,R-dZ,R+dZ,
	    ZAbscissae[mpgd][0],ZAbscissae[mpgd][1],gap,staveType,(int)volumeMap.size(),sN);
	  */
	}
      }
      if (error&0x4) {
	int module = error>>8;
	printf("** parseGeometry: MPGD 0x%x Inconsistency: Volume of Layer,Module %d,%d can't be cast into a \"%s\"\n",
	       0x1<<mpgd,layer,module,shapes[mpgd]);
      }
    }
    if (error) return false;
  } // End loop on mpgds
  return true;
}
