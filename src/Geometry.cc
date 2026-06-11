// ********** GEOMETRY

#include "Geometry.h"
using namespace std;

Geometry::Geometry(const char *geoFN, int verbose)
{
  // ***** CONSTRUCTOR
  // - Explore the TGeometry to:
  //  i) Determine which of the MPGDs are 2DStrip (in fact what's what's rather
  //    determined is whether they have MultipleSensitiveVolume, which as of
  //    2026/04, is the only way MPGDs can be 2DStrip).
  // => Stored in data member "MPGD2DStrips"
  // ii) Retrieve the World-To-Local transforms.
  // => Stored in "geoDetMatrix[detector][layer][module]".
  // - Assumptions (some effort is made to limit their number):
  //   - All MPGD modules can be reached via the TGeometry tree using a minimal
  //    set of tags (see "findLayersModules").
  //     Layer's node name includes "(L|l)ayer". Module's include "(M|m)odule".
  //     Sensitive surface is on a node named...
  //    ..."ReferenceThinGap" for 2DStrip (in fact MSV) MPGDs,
  //    ..."DriftGap" otherwise.
  //   - Encoding of layer# and module# numbering in IDDescriptor starts @ 0 (
  //    as opposed to @ 1).
  //   - Possibly others...
  
  // Init
  MPGD2DStrips = 0;
  if(address == 0) //!< Protection against multiple instances
    address = this;
  verboseLevel = verbose;

  int error = 0;
  // Get the TGeometry.
  // - Let's first check that the file exists in PWD, in order to avoid the
  //  possible warning message printed out by "TGeoManager::Import" when the
  //  file is missing.
  //const char geoFN[] = "hit_matching.geometry.root";
  TSystemDirectory *sys = new TSystemDirectory(".",".");
  if (!sys->GetListOfFiles()->FindObject(geoFN)) {
    printf("** Geometry: No \"%s\" file in $PWD.\n",geoFN);
    error = 1;
  }
  if (!error) {
  // - Import the TGeometry
    gGeoManager->SetVerboseLevel(verbose);
    if (!gGeoManager->Import(geoFN)) {
      printf("Error Importing \"%s\" into gGeoManager. => All MPGDs to be processed as pixels\n",
	     geoFN);
      error = 1;
    }
  }

  if (!error) {
    // Set "MPGD2DStrips" and "geoDetMats".
    for (int mpgd = 0; mpgd<N_MPGDs; mpgd++) {
      bool isMSV = false;
      if (getGeoMats(mpgd,isMSV)) {
	if (isMSV) MPGD2DStrips |= 0x1<<mpgd;
      }
    }
  }
  else {
    address = 0;
    throw std::runtime_error("Failing to init Geometry. => No World<->Local available.");
  }
};
bool findLayersModules(int LM,       // 0: Layer, 1: Module 
		       string &path, // Initial path
		       const char *tag1, const char *tag2,
		       int verbose,
		       vector<string>& targetNodes);
bool findDriftGap(const char *modulePath, int verbose,
		  const TGeoHMatrix *&m, bool &isMSV);
bool Geometry::getGeoMats(int mpgd, bool &isFullyMSV)
{
  // 

  isFullyMSV = true;

  vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[mpgd];

  const char *tag1 = tag1s[mpgd], *tag2 = tag2s[mpgd];
  vector<string> layerNodes;
  string path("/"); gGeoManager->CdTop();
  path += string(gGeoManager->GetCurrentNode()->GetName());
  bool ok = findLayersModules(0,path,tag1,tag2,verboseLevel,layerNodes);
  if (ok) {
    int nLayers = layerNodes.size();
    if (verboseLevel>1)
      printf("\"%s%c%s\": Found %d Layer(s): \n",
	     tag1,tag2?',':'\0',tag2?tag2:"\0",nLayers);
    for (int layer = 0; layer<nLayers; layer++) {
      const char *layerPath = layerNodes[layer].c_str();
      if (verboseLevel>1)
	printf("%2d/%2d: \"%s\"\n",layer,nLayers,layerPath);
      vector<string> moduleNodes;
      path = layerNodes[layer]; gGeoManager->cd(path.c_str());
      ok = findLayersModules(1,path,tag1,tag2,verboseLevel,moduleNodes);
      if (!ok) break;
      int nModules = moduleNodes.size();
      if (verboseLevel>1)
	printf("\"%s%c%s\", Layer %d: Found %d Modules: \n",
	       tag1,tag2?',':'\0',tag2?tag2:"\0",layer,nModules);
      vector<TGeoHMatrix> gMats;
      for (int module = 0; module<nModules; module++) {
	const char *modulePath = moduleNodes[module].c_str();
	if (verboseLevel>1)
	  printf("%2d/%2d: \"%s\"\n",module,nModules,modulePath);
	const TGeoHMatrix *m = 0; bool isMSV = false;
	if (findDriftGap(modulePath,verboseLevel,m,isMSV)) {
	  // cm -> mm
	  TGeoHMatrix mmm(*m);
	  const double *trcm = mmm.GetTranslation(); double trmm[3];
	  for (int xyz = 0; xyz<3; xyz++) trmm[xyz] = trcm[xyz]*10;
	  mmm.SetTranslation(trmm);
	  gMats.push_back(mmm);
	  isFullyMSV &= isMSV;
	}
	else {
	  if (verboseLevel) 
	    printf("** getGeoMats: Could not get TGeoHMatrix of Drift Gap in node \"%s\"\n",
		   modulePath);
	}
      }
      geoMats.push_back(gMats);
    }
  }
  if (!ok) {
    printf("** getGeoMats: No Module found for \"%s%c%s\"\n",
	   tag1,tag2?',':'\0',tag2?tag2:"\0");
    isFullyMSV = false;
  }
  if (verboseLevel) {
    int nLayers = geoMats.size();
    if (verboseLevel)
      printf("\"%s%c%s\": Found %d Layers: \n",
	     tag1,tag2?',':'\0',tag2?tag2:"\0",nLayers);
    for (int layer = 0; layer<nLayers; layer++) {
      vector<TGeoHMatrix>& gMats = geoMats[layer];
      int nModules = gMats.size();
      if (verboseLevel)
	printf("\"%s%c%s\", Layer %d: Found %d Modules: \n",
	       tag1,tag2?',':'\0',tag2?tag2:"\0",layer,nModules);
      for (int module = 0; module<nModules; module++) {
	TGeoHMatrix& geoMat = gMats[module];
	const double *tr = geoMat.GetTranslation();
	if (verboseLevel>1)
	  printf("%2d/%d: %8.4f,%8.4f,%8.4f cm\n",
		 module,nModules,tr[0],tr[1],tr[2]);
      }
    }
  }
  return true;
}
bool findLayersModules(int LM, string &path,
		       const char *tag1, const char *tag2,
		       int verbose,
		       vector<string>& targetNodes)
{
  // Find all TGeoNodes w/ name including:
  //  i) input tags ("tag2" can be =0, it's needed for Endcaps which share a
  //    common "assembly" node),
  //  AND if LM=0, the word "MPGD" (it's needed for the earlier levels, to
  //    raise ambiguities),
  // ii) LM=0: the word "Layer" (or "layer").
  //     LM=1: the word "Module" (or "module").
  // Collect their path in "targetNodes".
  // - Starting point is input "path"
  // - Then node fulfilling (i) is recursively searched...
  // - ...until nodes fulfilling (ii) are found (in the mean time, only one
  //  node is allowed, else we wouldn't know which path to then follow).
  // This strategy is presumed to work in all kinds of circumstance, possibly
  // differing from the situation as of 2026/04, where we have, e.g.:
  //  "/world_volume_1
  //    /OuterBarrelMPGDSubAssembly_10
  //      /MPGDOuterBarrel_0
  //        /MPGDOuterBarrel_layer0_0
  //          /MPGDOuterBarrelModule_0,2,...,23"
  // A maximum depth of 10 is still imposed in addition.
  const TGeoNode *node = gGeoManager->GetCurrentNode(), *found;
  for (int depth = 0; depth<10; depth++) {
    const TObjArray* subNodeArray = node->GetNodes();
    TObject *o = subNodeArray->First();
    found = 0; while (o) {
      if (!o->IsA()->InheritsFrom(TGeoNode::Class())) {
	printf("** findLayersModules: Object \"%s\" of Class \"%s\" among daughter of TGeoNode \"%s\"\n",
	       o->GetName(),o->ClassName(),node->GetName());
	return false;
      }
      const TGeoNode *subNode = (TGeoNode*)o;
      const char *name = subNode->GetName();
      if (verbose>2) {
	const TGeoMatrix *m = subNode->GetMatrix();
	const double *tr = m->GetTranslation();
	printf("  \"%s\": %8.4f,%8.4f,%8.4f cm\n",name,tr[0],tr[1],tr[2]);
      }
      // Apply requirement (i)
      if ((LM==1 || strstr(name,"MPGD")) &&
	  (strstr(name,tag1) || tag2 && strstr(name,tag2))) {
	bool isTarget = LM==0 ? strstr(name,"Layer")  || strstr(name,"layer")
	  /* */               : strstr(name,"Module") || strstr(name,"module");
	if (!isTarget && found) {
	  printf("** findLayersModules: Ambiguity: TGeoNode \"%s\"(=> not a %s) has more than one subNode: \"%s\" and \"%s\" ",
		 node->GetName(),LM==0?"Layer":"Module",found->GetName(),name);
	  return false;
	}
	if (isTarget) // Requirement (ii): Is Add path to the list
	  targetNodes.push_back(path+string("/")+string(subNode->GetName()));
	found = subNode;
      }
      o = subNodeArray->After(o);
    }
    if (!found) return false;
    if (targetNodes.empty()) { // Increment "path" and keep searching
      path += string("/") + string(found->GetName());
      node = found;
      if (verbose>1) {
	const TGeoMatrix *m = found->GetMatrix();
	const double *tr = m->GetTranslation();
	printf("\"%s\": %8.4f,%8.4f,%8.4f cm\n",path.c_str(),tr[0],tr[1],tr[2]);
      }
    }
    else
      break;
  }
  return !targetNodes.empty();
}
bool findDriftGap(const char *modulePath, int verbose,
		  const TGeoHMatrix *&m, bool &isMSV)
{
  // Return TGeoHMatrix of World-to-Local of Reference SubVolume, which name
  // is assumed to be "ReferenceThinGap".
  if (!gGeoManager->cd(modulePath)) {
    printf("** findDriftGap: Bad input path \"%s\"\n",
	   modulePath);
    return false;
  }
  TGeoNode *node = gGeoManager->GetCurrentNode();
  const TObjArray* subNodeArray = node->GetNodes();
  TObject *o = subNodeArray->First();
  const TGeoNode* found = 0; while (o) {
    if (!o->IsA()->InheritsFrom(TGeoNode::Class())) {
      printf("** findDriftGap: Object \"%s\" of Class \"%s\" among daughter of TGeoNode \"%s\"\n",
	     o->GetName(),o->ClassName(),node->GetName());
      return false;
    }
    const TGeoNode *subNode = (TGeoNode*)o;
    const char *name = subNode->GetName();
    if (verbose>2) {
      const TGeoMatrix *m = subNode->GetMatrix();
      const double *tr = m->GetTranslation();
      printf("  \"%s\": %8.4f,%8.4f,%8.4f cm\n",name,tr[0],tr[1],tr[2]);
    }
    isMSV = strstr(name,"ReferenceThinGap");
    if (isMSV || strstr(name,"DriftGap")) {
      found = subNode; break;
    }
    o = subNodeArray->After(o);
  }
  if (!found) return false;
  string nodePath = string(modulePath)+string("/")+string(found->GetName());
  if (!gGeoManager->cd(nodePath.c_str())) {
    printf("** findDriftGap: Inconsistency: gGeoManager cannot cd(\"%s\")\n",
	   nodePath.c_str());
    return false;
  }
  m = gGeoManager->GetCurrentMatrix();
  if (verbose>1) {
    const double *tr = m->GetTranslation();
    printf("\"%s\": %8.4f,%8.4f,%8.4f cm\n",found->GetName(),tr[0],tr[1],tr[2]);
  }
  return true;
}
void cellID2LayerModule(int idet, unsigned long cellID,
			unsigned int &layer, unsigned int &module)
{
  // Assuming IDDescriptor to be:
  if (idet<2) { // Barrel:
    // <id>system:8,layer:4,module:12,         strip:28:4,u:-16,v:-16</id>
    layer =  cellID>>8&0xf, module = cellID>>12&0xfff;
  }
  else {        // ECT
    // <id>system:8,layer:2,module:6,sensor:16,strip:28:4,x:32:-16,z:-16</id>
    layer =  cellID>>8&0x3, module = cellID>>10&0x3f;
  }
}
bool Geometry::WorldToLocal(int idet, unsigned long cellID,
			    double *gpos, // In mm
			    double *lpos) // In mm
{
  unsigned int layer, module; cellID2LayerModule(idet,cellID,layer,module);
  bool ok = false; 
  if (idet<N_MPGDs) {
    vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[idet];
    size_t nLayers = geoMats.size(), nModules = 0; if (layer<nLayers) {
      if (module<geoMats[layer].size()) ok = true;
    }
    if (ok) {
      const TGeoHMatrix &geoMat = geoMats[layer][module];
      geoMat.MasterToLocal(gpos,lpos);
    }
  }
  if (!ok) {
    printf("** Geometry::WorldToLocal: No Geometry found for mpgd=%d/%d, cellID = 0x%08lx => (layer,module) = (%u,%u)\n",idet,N_MPGDs,cellID&0xffffffff,layer,module);
    return false;
  }
  return true;
}
bool Geometry::WorldToLocal(int idet, unsigned long cellID,
			    double *gpos, double *gmom, // In mm
			    double *lpos, double *lmom) // In mm
{
  unsigned int layer, module; cellID2LayerModule(idet,cellID,layer,module);
  bool ok = false; 
  if (idet<N_MPGDs) {
    vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[idet];
    size_t nLayers = geoMats.size(), nModules = 0; if (layer<nLayers) {
      if (module<geoMats[layer].size()) ok = true;
    }
    if (ok) {
      const TGeoHMatrix &geoMat = geoMats[layer][module];
      geoMat.MasterToLocal(gpos,lpos);
      geoMat.MasterToLocalVect(gmom,lmom);
    }
  }
  if (!ok) {
    printf("** Geometry::WorldToLocal: No Geometry found for mpgd=%d/%d, cellID = 0x%08lx => (layer,module) = (%u,%u)\n",idet,N_MPGDs,cellID&0xffffffff,layer,module);
    return false;
  }
  return true;
}
bool Geometry::LocalToWorld(int idet, unsigned long cellID,
			    double *lpos, // in mm
			    double *gpos) // in mm
{
  unsigned int layer, module; cellID2LayerModule(idet,cellID,layer,module);
  bool ok = false; 
  if (idet<N_MPGDs) {
    vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[idet];
    size_t nLayers = geoMats.size(), nModules = 0; if (layer<nLayers) {
      if (module<geoMats[layer].size()) ok = true;
    }
    if (ok) {
      const TGeoHMatrix &geoMat = geoMats[layer][module];
      geoMat.LocalToMaster(lpos,gpos);
    }
  }
  if (!ok) {
    printf("** Geometry::LocalToWorld: No Geometry found for mpgd=%d/%d, cellID = 0x%08lx => (layer,module) = (%u,%u)\n",idet,N_MPGDs,cellID&0xffffffff,layer,module);
    return false;
  }
  return true;
}
