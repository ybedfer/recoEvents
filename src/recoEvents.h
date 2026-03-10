//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  5 23:53:48 2024 by ROOT version 6.30/02
// from TChain events/
//////////////////////////////////////////////////////////

/*
// ********** USAGE
// ***** INPUT = TChain
const char tag[7] = "struv", nEvtsTag[3] = "20"; 
TChain *events = new TChain("events"); char fName[] = "podio.sensor.1234.20,root"; size_t lN = strlen(fName)+1; int seeds[] = {1234,4567,8910,1112}; int nSeeds = sizeof(seeds)/sizeof(int);
for (int i = 0; i<nSeeds; i++) { snprintf(fName,lN,"podio.%s.%s.%4d.root",tag,nEvtsTag,seeds[i]); events->Add(fName); }
// ***** INPUT = TTree
TTree *events = (TTree*)gDirectory->Get("events");
// ***** INSTANTIATION
.L ../recoEvents/install/librecoEvents.so
recoEvents ana(events,0xf); // Instantiate for all of (0x1:CyMBaL,0x2:Outer,0x4:Vertex,0x8:Si)
// ***** EVENT CONTROL, DEBUGGING
ana.select = new TTreeFormula("select", "@MCParticles.size()==1", events);// Add rejection cut
ana.requirePDG = 13;        // Require MCParticle = mu-
ana.requireQuality = 2;     // Require primary (!=0) and reject primary w/ interfering secondary in same module (>1)
ana.verbose = 0x1111;       // Debugging printout (1 for CyMBaL, 2 for Outer...)
// ***** 5-SUBVOLUME: is default. SimHits are coalesced alla MPGDTrackerDigi
recoEvents ana(events,0xf,0x0); // Instantiate w/ no strips (i.e. pixels) in CyMBaL and Outer
ana.SetNSensitiveSurfaces(1);   // Overwrite default.
// ***** LOOP
ana.Loop();                 // For debugging: "Loop(<nEvents>,<firstEvent>)"
// ***** DRAW
ana.DrawphithZR(0,0xf,true); // Draw SimHits, w/ if true, colour highlighting of module type.
ana.DrawphithZR(1,0xf,true); // Draw RecHits (1: 1st coord, 2: 2nd coord), w/ if true, colour highlighting of module type.
ana.DrawResiduals(1,0x1);    // Draw residuals RecHit-SimHit for phi (iRec=1) of CyMBaL (=0x1)
ana.DrawSimHit(3328,0x1,0);  // Draw 0th SimHit of detector 0x1 for evt #3328
ana.DrawSimHit(3328,0x1,4,1);// Superimpose 4th SimHit
// Direct access to histograms
ana.recHs[0][0].ZR[0]; h2->SetTitle("CyMBaL Rec"); h2->Draw(); SetPaveText(h2,0);
// ***** SEVERAL recoEvents OBJECTS
recoEvents ana2(events2,0xf);
// [...]
ana.DrawResiduals(1,0x1,0x6)
ana2.DrawResiduals(1,0x1,0x6,cCyMBaLRes,3)
// ***** PDF
cCyMBaLRes->Print("cCyMBaLRes.pdf","EmbedFonts")
// ***** Scan
t->Scan("EventHeader.eventNumber:@MPGDBarrelHits.size():MPGDBarrelHits.cellID&0xffffffff:MPGDBarrelHits.quality:MPGDBarrelHits.cellID>>32:MPGDBarrelHits.momentum.z:@MPGDBarrelRawHits.size():MPGDBarrelRawHits.cellID&0xffffffff:MPGDBarrelRawHits.cellID>>32:@MPGDBarrelRecHits.size():MPGDBarrelRecHits.cellID&0xffffffff","@MPGDBarrelHits.size()","col=3d:2d:8llx:2d:8llx:6.2f:2d:8llx:8llx:2d:8llx",1,1555);
t->Scan("EventHeader.eventNumber:@OuterMPGDBarrelHits.size():OuterMPGDBarrelHits.cellID&0xffffffff:OuterMPGDBarrelHits.quality:OuterMPGDBarrelHits.cellID>>32:OuterMPGDBarrelHits.momentum.x:OuterMPGDBarrelHits.momentum.y:OuterMPGDBarrelHits.momentum.z:@MCParticles.size():MCParticles[2].MCParticles.PDG","@OuterMPGDBarrelHits.size()","col=4d:2d:8llx:2d:8llx:6.2f:6.2f:6.2f:2d:2d",1,14724);
*/

#ifndef recoEvents_h
#define recoEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeFormula.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "edm4eic/MCRecoClusterParticleAssociationData.h"
#include "podio/ObjectID.h"
#include "edm4eic/ClusterData.h"
#include "edm4hep/RawCalorimeterHitData.h"
#include "edm4eic/CalorimeterHitData.h"
#include "edm4eic/MCRecoTrackerHitAssociationData.h"
#include "edm4hep/SimTrackerHitData.h"
#include "edm4eic/RawTrackerHitData.h"
#include "edm4eic/TrackerHitData.h"
#include "edm4eic/TrackSegmentData.h"
#include "edm4eic/TrackParametersData.h"
#include "edm4eic/TrackData.h"
#include "edm4eic/TrajectoryData.h"
#include "edm4eic/Measurement2DData.h"
#include "edm4eic/VertexData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4eic/CherenkovParticleIDData.h"
#include "edm4hep/EventHeaderData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4eic/HadronicFinalStateData.h"
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4eic/MCRecoParticleAssociationData.h"
//#include "edm4hep/TrackerHitData.h"

using namespace std;
using namespace edm4hep;

typedef struct{ TDirectory *dir; TH2D *X, *Y, *Z, *R, *phi, *phir, *th, *mod, *thphi, *XY, *ZR, *Rr, *xyr, *Ur, *Vr, *eDep; }
  Histos;
typedef struct{ TDirectory *dir; TH1D *X, *Y, *Z, *phi, *R, *Rphir, *Ur, *Vr; /* TProfile2D *xyr; */}
  Resids;

void SetPaveText(TH1 *h, int mode = 0, int opt = 0);

class LayerModules {
  // Structure encoding layer# and module#
  // (Note: substructures ("modulePat" and index into it) do not match exactly
  // layer and module entities defined in IDDescriptor.) 
public:
  LayerModules() {
    init();
  };
  LayerModules(int idet, unsigned long cellID) {
    int num, index;
    if        (idet==2|| idet==3) {
      index = (cellID>>8&0x3)-1; num = cellID>>10&0x3f;
    } else if (idet==0 || idet==1) {
      index = 0;                 num = cellID>>12&0xfff;
    } else if (idet==4) { // Vertex: 1920 modules (at most) => divide by 32
      index = (cellID>>8&0xf);   num = cellID>>12&0xfff; num /= 32;
    }
    else if (idet==5) { // Si: 2 systemIDs: 0x3b and 0x3c; sublayer = module parity
      /* */                      num = cellID>>12&0xfff;
      index = (1-(cellID%2))*2 + num%2; num /= 2;
    }
    if (index>3 || num>63) {
      printf("** LayerModules: Invalid args: idet,cellID = %d,0x%08lx,0x%08lx => index,num = %d,%d\n",
	     idet,cellID&0xffffffff,cellID>>32,index,num);
    }
    for (int idx = 0; idx<4; idx++) {
      if (idx==index) patterns[idx] = ((unsigned long)0x1)<<num;
      else            patterns[idx] = 0;
    }
  }
  void operator=(LayerModules &lm) {
    for (int idx = 0; idx<4; idx++) patterns[idx] = lm.patterns[idx];
  };
  ~LayerModules() {};
  void init() {
    for (int idx = 0; idx<4; idx++) patterns[idx] = 0;
  };
  void add(LayerModules &lm) {
    for (int idx = 0; idx<4; idx++) patterns[idx] |= lm.patterns[idx];
  };
  void subtract(LayerModules &lm) {
    for (int idx = 0; idx<4; idx++) patterns[idx] &= ~lm.patterns[idx];
  };
  bool contains(LayerModules &lm) {
    for (int idx = 0; idx<4; idx++) {
      if (lm.patterns[idx]&patterns[idx]) return true;
    }
    return false;
  };
  unsigned long patterns[4];
};

class recoEvents: public TNamed {
public :
   recoEvents(TTree *tree,
	      // Detector pattern: 0x1=CyMBaL, 0x2=Outer, 0x4=Vertex, 0x8=Si
	      // mode: Processed detectors
	      // hasStrips: If not set, detector has pixels.
	      unsigned int mode = 0xf, unsigned int hasStrips = 0x3);
   virtual ~recoEvents();

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // "nObjCreated": counting # of objects in a STATIC MEMBER. "iObjCreated"
   // is specific to an instance and used to give distinct names to its TCanvas'
   // or TDirectory's.
   static unsigned int nObjCreated; unsigned int iObjCreated;
   unsigned int processedDetectors;
   unsigned int stripMode;  // Pattern of detectors w/ strip readout
   bool reconstruction;
   bool requireTraversing; // Conditions coalescing
   int getDetHit(int idet, int ih, double &X, double &Y, double &Z,
		 unsigned int &module, unsigned int &div);
   // ***** EVENT CONTROL, DEBUGGING
   int evtNum; // Current event# (used to document error messages).
  // (Layer,Module)
   LayerModules prvLayerModule;
   LayerModules allLayerModules; // Pattern of all (layer,module)'s where hit
   // (Layer,module)'s where only LONE PRIMARIES, i.e. primary hits w/o secondary offsprings
   LayerModules prvLayerModuleLP;
   LayerModules layerModulesLP;  // All (layer,module)'s where only LPs
   LayerModules layerModulesLP1; // where in addition single LP
   // Debugging
   unsigned int doDebug; int debuggedDet; // Debug control
   int evtToDebug; // "evtNum" of event to be debugged
   // Verbosity:
   // 0xh<<Idet, h = 0x1: SimHit/RecHit, 0x10: Association, 0x100: Coalescing/Extending, 0x1000: Large residuals
   // 0x10000: More info
   // 0x20000: Header
   unsigned int verbose;

   // ***** DETECTOR NAMES
#define N_DETs 6
   const char *detectorNames[N_DETs];
   // STRIP SEGMENTATION
   const char *coordNames[2][N_DETs];
   void SetNSensitiveSurfaces(int nSurfaces);
   // GEOMETRY/CONFIGURATION
   // CyMBal: 4 sections * 8 staves, numbering from 0
   unsigned int MPGDs, Barrels;
   bool isMPGD(int idet), isBarrel(int idet);
   static constexpr int nModules[N_DETs] =   {32,24,48,48,1920,70};
   static constexpr int moduleMns[N_DETs] = { 0, 0, 1, 1,  1, 1};
   double volumeThicknesses[N_DETs];   // Overall thickness
   double radiatorThicknesses[N_DETs]; // Thickness of RADIATOR SUBVOLUME
   vector<double> radii[N_DETs];      // in mm
   vector<double> ZAbscissae[N_DETs]; // in mm
   int nSections[N_DETs];
   vector<double> sectionDZs[N_DETs]; // Transform global -> local
   vector<double> hWidths[N_DETs];    // HalfWidths
   double ZHLengths[N_DETs];          // HalfLengths
   vector<double> pitches[N_DETs];    // pitch in mm
   double gains[N_DETs], eDThresholds[N_DETs], resolutions[N_DETs];

   // ***** SELECTION
   TTreeFormula *select; // Provides for specifying a rejection cut.
   // ***** EVENT SELECTION
   int requireNMCs;
   // ***** MODULE SELECTION
   bool requireModules;
   void AddRequiredLayerModules(int idet, int index, unsigned long pattern);
   // ***** HIT SELECTION
   // Requirements for filling "simHs" and residuals
   int requirePDG;     // Associated MCParticle 
   int requireQuality; // "SimTrackerHit::quality". =1: reject modules where any secondary hit, >1: reject modules where more than one hit
   int requireOffEdge; //  SimHit (or coalesced) not on edge
  
   // ********** HISTOS
   void BookHistos(Histos *Hs, const char *tag);
   void getxyArgs(bool isRec, int idet, double &xMx, double &yMx, string &sT);
   unsigned int getStatus(int idet, int ih);
   unsigned int getStatus(int idet, int ih, map<int,int> &sim2coa);
   MCParticleData& getMCParticle(int idet, int ih);
   bool getrec2coas(int idet, map<int,int> &sim2coa,
		    map<int,vector<int>,less<int>> &rec2coas);
   void fillHit(int iSimRec, int idet,
		double X, double Y, double Z, double E, unsigned long cellID);
   bool fillResids(int idet,
		   const Vector3f &pos, const Vector3d &psim, unsigned long cellID);
   unsigned int isOnEdge(int idet, unsigned long cellID,
			 double Xs, double Ys, double Zs);
   double getResCut(int idet, int module, int strip, bool isOnEdge);
   void parseCellID(int idet, unsigned long ID,
		    unsigned int &module, unsigned int &div, unsigned int &strip);
   bool parseStrip(int idet, int simOrRec, unsigned int &strip);
   void g2lCyMBaL(double X,   double Y,   double Z, unsigned int div,
		  double &Xr, double &Yr, double &Zr,
		  double &Rr, double &phir);
   void g2lCyMBaL(double Px, double Py, double Pz, unsigned int div,
		  double *lmom);
   void l2gCyMBaL(double *lpos, unsigned int div, double *gpos);
   void g2lOuter(double X, double Y, double Z, unsigned int module,
		 double &Rcphi, double &Xr, double &Yr, double &Zr,
		 double &Ur, double &Vr);
   void g2lOuter(double Px, double Py, double Pz, unsigned int module,
		 double *lmom);
   void l2gOuter(double *lpos, unsigned int module, double *gpos);
   bool samePMO(int idet, int ih, int jh);
   bool extrapolate(int idet, int ih, int jh);
   bool coalesce(int idet, vector<int> coalesced, SimTrackerHitData &hext);
   void extend(int idet, int ih, SimTrackerHitData &hext);
   bool checkTraversing(int idet, SimTrackerHitData &hit,
			double &pathDepth, double &depth);
   bool extendHit(int idet, int ih, int direction, double *lext);
   // ***** EVENT CONTROL, DEBUGGING
   void initDetEvent();
   void requestDebug(int idet, unsigned int level);
   bool debugIsOn(unsigned int level);
   void updateDetEvent(int idet, int ih);
   void finaliseDetEvent();
   bool moduleSelection(int idet, unsigned long cellID);
   void printHit(int idet,
		 double X, double Y, double Z, unsigned long cellID);
   void debugHit(int idet, int ih, int nHs, SimTrackerHitData &hit, unsigned int status);
   void debugHit(int idet, SimTrackerHitData &hit);
   void debugRec(int idet, int ir);
   void debugHitRec(int idet, int is, int ncoas, int cIndex, SimTrackerHitData &hit, int ir);
   void debugAssoc(int idet);
   void debugAssoc(int idet, map<int,int>raw2rec, map<int,int> sim2coa, map<int,vector<int>> rec2sims);
   static constexpr int violet = kMagenta+2, bleu = kBlue+1, vert = kGreen+2;
   static constexpr int jaune = kOrange+0, orange = kOrange+1, rouge  = kRed+0;
   static constexpr int marron = 50, gris   = 11;
   void DrawphithZR(int iSimRec = 0 /* 0: sim, 1: rec, 2: rec 2nd coord of STRIP */,
		    unsigned int detectorPattern = 0x3f, 
		    unsigned int histoPattern = 0xf, bool decompose = false,
		    TCanvas *cPrv = 0, int ipad = 1, int col = -1);
   void DrawModules(int iSimRec = 0, unsigned int detectorPattern = 0x1, bool decompose = false);
   void DrawResiduals(int iRec = 1,   unsigned int detectorPattern = 0x1,
		      unsigned int histoPattern = 0xf,
		      TCanvas *cPrv = 0, int ipad = 1, int col = -1);
   void SetMinima(double min);
   TDirectory *dSim, *dSimDets[N_DETs];
   TDirectory *dRec, *dRecDets[N_DETs];
   Histos simHs[N_DETs], recHs[2][N_DETs];
   Resids resHs[2][N_DETs];
   TH1D *hMult;
   void DrawSimHit(int jentry, unsigned int detectorPattern, int ih,
		   bool addToPreExisting = false);

   // ***** BRANCHES
   vector<SimTrackerHitData> *hits[N_DETs];
   TBranch *simBranches[N_DETs];
   vector<podio::ObjectID> *amcs[N_DETs];
   TBranch *amcBranches[N_DETs];
   vector<edm4eic::TrackerHitData> *recs[N_DETs];
   TBranch *recBranches[N_DETs];
   vector<podio::ObjectID> *arhs[N_DETs];
   TBranch *arhBranches[N_DETs];
   vector<podio::ObjectID> *ashs[N_DETs];
   TBranch *ashBranches[N_DETs];
   vector<podio::ObjectID> *aRhs[N_DETs]; // RECHIT -> RAWHIT
   TBranch *aRhBranches[N_DETs];
   vector<MCParticleData> *mcParticles;
   TBranch *MCParticles;
   vector<EventHeaderData> *eventHeader;
   TBranch *EventHeader;

  virtual void     Loop(int nEvents = 0, int firstEvent = 0);

 private:
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
  
   void initGeometry(int idet, bool hasStrips);

   // ***** Data members not to be modified from command line:
   // Depending upon setup, MPGDs can have more than one sensitive surfaces.
   // - Can be =1 or =5.
   // - But source code for =1 may not be up-to-date.
   int nSensitiveSurfaces;
   // ***** MODULE SELECTION
   LayerModules requiredLayerModules[N_DETs];

   ClassDef(recoEvents,1); // Must be the last item before the closing '};'
};

#endif

#ifdef recoEvents_cxx
recoEvents::recoEvents(TTree *tree, unsigned int detectors, unsigned int hasStrips) : fChain(0) 
{
  // Init global settings
  processedDetectors = detectors; verbose = 0; select = 0;
  // ***** SELECTION
  requireNMCs =    0; // Default: no requirement
  requireModules = 0; // Default: all modules allowed
  for (int idet = 0; idet<N_DETs; idet++) requiredLayerModules[idet].init();
  requirePDG =     0; // Default: do not require any ID
  requireQuality = 0; // Default: do not require "SimTrackerHit::quality"
  requireOffEdge = 0; // Default: do not require off edge
  stripMode = hasStrips;
  requireTraversing = false; // Default = coalesce all extrapolate-compatible hits
  // ***** GEOMETRY
  MPGDs = 0xf;    // CyMBaL, Outer and Endcaps
  // Check that arg. "hasStrips" fits the pattern of MPGDs
  if ((hasStrips&MPGDs)!=hasStrips) {
    printf("** recoEvents: Invalid arg. <hasStrips>(=0x%x): does not fit pattern of MPGDs(=0x%x)! => Aborting...\n",
	   hasStrips,MPGDs);
    return;
  }
  Barrels = 0x33; // CyMBaL, Outer and Vertex and Si
  nSensitiveSurfaces = 5; // Default = 5. Can be changed via "SetNSensitiveSurfaces"
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&processedDetectors)) continue;
    initGeometry(idet,hasStrips&0x1<<idet);
  }
  // ***** OBJECT ID
  iObjCreated = nObjCreated++;
  evtToDebug = -1;
  Init(tree);
}

recoEvents::~recoEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t recoEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t recoEvents::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void recoEvents::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) {
     printf("** recoEvents::Init: Arg. TTree* is null. => Abort!\n");
     return;
   }
   fChain = (TChain*)tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1); // Leads to run time crash as soon as one act on mpgdHs
   // Set margins, in order for recoEvents axis settings to work well
   gStyle->SetPadBottomMargin(.15); gStyle->SetPadLeftMargin(.14);

   // ********** BRANCHES
   // ***** SIMHITS
   for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "MPGDBarrel","OuterMPGDBarrel","BackwardMPGDEndcap","ForwardMPGDEndcap",
       "VertexBarrel","SiBarrel"};
     const char *branchName = branchNames[idet];
     string name(branchName); name += string("Hits");
     fChain->SetBranchAddress(name.c_str(),&hits[idet],&simBranches[idet]);
   }
   // ***** MC Particle
   for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "_MPGDBarrel","_OuterMPGDBarrel","_BackwardMPGDEndcap","_ForwardMPGDEndcap",
       "_VertexBarrel","_SiBarrel"};
     const char *branchName = branchNames[idet];
     //#define MCPARTICLE
#ifdef MCPARTICLE
     string name(branchName); name += string("Hits_MCParticle");
#else
     string name(branchName); name += string("Hits_particle");
#endif
     fChain->SetBranchAddress(name.c_str(),&amcs[idet],&amcBranches[idet]);
   }
   // ***** REC BRANCHES
   // Determine whether available for all requested detectors
   unsigned int recBranchPat = 0; for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "MPGDBarrel","OuterMPGDBarrel","BackwardMPGDEndcap","ForwardMPGDEndcap",
       "SiBarrelVertex","SiBarrelTracker"};
     const char *branchName = branchNames[idet];
     string name(branchName); name += string("RecHits");
     if (fChain->GetBranch(name.c_str())) recBranchPat |= 0x1<<idet;
   }
   // "reconstruction"? It's either all or nothing
   reconstruction = (recBranchPat&processedDetectors)==processedDetectors;
   if (!reconstruction && recBranchPat) {
     printf(" * Init: Some detectors (0x%x) have RecHits TBranches in input TTree, but not all of the 0x%x requested.\n => Only SimHits will be processed.\n",
	    recBranchPat,processedDetectors);
   }
   else if (reconstruction) {
     printf(" * Init: RecHits TBranches in input TTree = 0x%x, matching requested detectors 0x%x.\n => Both SimHits and RecHits will be processed.\n",
	    recBranchPat,processedDetectors);
   }
   else {
     printf(" * Init: No RecHits TBranch in input TTree matching requested detectors 0x%x.\n => Only SimHits will be processed.\n",
	    processedDetectors);
   }
   if (reconstruction) { // Simulation Only" => bypass RecHits and associations
     // ********** RECHITS
     for (int idet = 0; idet<N_DETs; idet++) {
       const char *branchNames[N_DETs] = {
	 "MPGDBarrel","OuterMPGDBarrel","BackwardMPGDEndcap","ForwardMPGDEndcap",
	 "SiBarrelVertex","SiBarrelTracker"};
       const char *branchName = branchNames[idet];
       string name(branchName); name += string("RecHits");
       fChain->SetBranchAddress(name.c_str(),&recs[idet],&recBranches[idet]);
     }
     // ***** RAWHITS <-> SIMHITS ASSOCIATION
     for (int idet = 0; idet<N_DETs; idet++) {
       const char *branchNames[N_DETs] = {
	 "_MPGDBarrel","_OuterMPGDBarrel","_BackwardMPGDEndcap","_ForwardMPGDEndcap",
	 "_SiBarrelVertex","_SiBarrel"};
       const char *branchName = branchNames[idet];
       string namr(branchName); namr += string("RawHitAssociations_rawHit");
       string nams(branchName); nams += string("RawHitAssociations_simHit");
       fChain->SetBranchAddress(namr.c_str(),&arhs[idet],&arhBranches[idet]);
       fChain->SetBranchAddress(nams.c_str(),&ashs[idet],&ashBranches[idet]);
     }
     // ***** RAWHIT   -> RECHIT  ASSOCIATION
     for (int idet = 0; idet<N_DETs; idet++) {
       const char *branchNames[N_DETs] = {
	 "_MPGDBarrel","_OuterMPGDBarrel","_BackwardMPGDEndcap","_ForwardMPGDEndcap",
	 "_SiBarrelVertex","_SiBarrelTracker"};
       const char *branchName = branchNames[idet];
       string namR(branchName); namR += string("RecHits_rawHit");
       fChain->SetBranchAddress(namR.c_str(),&aRhs[idet],&aRhBranches[idet]);
     }
   }
   fChain->SetBranchAddress("MCParticles",&mcParticles,&MCParticles);
   fChain->SetBranchAddress("EventHeader",&eventHeader,&EventHeader);

   Notify();

   // ***** DETECTOR NAMES
   const char *detNs[N_DETs] = {"CyMBaL","Outer","BECT","FECT","Vertex","Si"};
   const char *coordNs[2][2] = {{"#scale[1.2]{#varphi}","#font[32]{Z}"},
				{"#font[32]{U}","#font[32]{V}"}};
   // ***** STRIPS
   for (int idet = 0; idet<N_DETs; idet++) {
     detectorNames[idet] = detNs[idet];
     if (0x1<<idet&stripMode) {
       coordNames[0][idet] = coordNs[idet][0];
       coordNames[1][idet] = coordNs[idet][1];
     }
   }

   // ********** HISTOS
   // ***** BASE DIRECTORY = tag it w/ recoEvents instance#
   TDirectory *dSave = gDirectory;
   gDirectory->cd("/");
   char dN[] = "dSim00"; size_t lN = strlen(dN)+1;
   int iObj = iObjCreated;
   if (iObj) snprintf(dN,lN,"dSim%d",iObj);
   else      snprintf(dN,lN,"dSim");
   dSim = gDirectory->mkdir(dN); dSim->cd();
   printf(" * Init: TDirectory \"%s\"\n",gDirectory->GetName());
   // ***** INSTANTIATION
   BookHistos(simHs,"s");
   dSave->cd();
   if (reconstruction) {
     // ***** RECO HISTOS
     if (iObj) snprintf(dN,lN,"dRec%d",iObj);
     else      snprintf(dN,lN,"dRec");
     dRec = gDirectory->mkdir(dN); dRec->cd();
     printf(" * Init: TDirectory \"%s\"\n",gDirectory->GetName());
     // ***** INSTANTIATION
     BookHistos(recHs[0],"r0");
     BookHistos(recHs[1],"r1");
     dSave->cd();
   }
}
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
    radii[0].push_back(556.755); radii[0].push_back(578.755);
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
  }
  else if (idet==1) {                    // ***** OUTER
    radii[1].push_back(737.4650);
    // <constant name="MPGDOuterBarrelModule_zmin1"     value="1795*mm"/>
    // <constant name="MPGDOuterBarrelModule_zmin2"     value="1845*mm"/>
    // <constant name="MPGDOuterBarrelModule_PCB_offset"              value="110*mm"/>
    // <constant name="MPGDOuterBarrelFrame_width"            value="15*mm"/>
    ZAbscissae[1].push_back(-1795+110+15); ZAbscissae[1].push_back(1845-110-15);
    sectionDZs[1].push_back(880); sectionDZs[1].push_back(-830);
    // <constant name="MPGDOuterBarrelModule_width"                   value="360*mm"/>
    hWidths[1].push_back(180-15);
    ZHLengths[1] = 840;
    //<constant name="MPGDOuterBarrelPitch"                        value = "800*um" />
    double outerPitch = 800; pitches[1].push_back(outerPitch/1000);
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
// *************** HISTOS
void setAxes(TAxis *ax, int nDiv, int maxDigits)
{
  ax->SetNdivisions(nDiv); ax->SetLabelFont(62);
  ax->SetLabelSize(.055);  ax->SetLabelOffset(.006);
  ax->SetTitleSize(.065);  ax->SetTitleOffset(1.0);
  if (maxDigits) ax->SetMaxDigits(maxDigits);
}
void recoEvents::BookHistos(Histos *Hs, const char* tag)
{
  using namespace ROOT;
  using namespace std;

  // ***** PARSE <tag> ARG.
  // Check internal consistency
  if (strcmp(tag,"s") && strcmp(tag,"r0") && strcmp(tag,"r1")) {
    printf("** BookHistos: Inconsistency: Invalid <tag> arg.(=\"%s\"). Aborting...\n",tag);
    exit(1);
  }
  bool isRec = tag[0]=='r';
  // ***** STRIPS
  int iStrip = 0;
  if (!strcmp(tag,"r1")) iStrip = 2; // 2nd coordinate, for STRIPS
  else                   iStrip = 1;

  // ***** BASE DIRECTORY
  TDirectory *dSave = gDirectory;
  if (!isRec) // Only once
    hMult = new TH1D("hMult","@MCParticles.size()",512,0,512);

  // ********** LOOP ON DETECTORS
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&processedDetectors)) continue;
    if (iStrip==2 && !(0x1<<idet&stripMode)) // 2nd coordinate: skip if not STRIP
      continue;
    Histos &hs = Hs[idet];
    // ***** SUB-DIRECTORY
    string sdS = string("d")+string(detectorNames[idet]);
    if ((0x1<<idet&stripMode) && iStrip==2) sdS += string("2"); 
    const char *sdN = sdS.c_str();
    TDirectory *dSD = dSave->mkdir(sdN); dSD->cd(); hs.dir = dSD;
    printf(" * BookHistos: TDirectory \"%s\"",gDirectory->GetName());
    if ((0x1<<idet&stripMode) && iStrip==2) {
      string sdT("2nd coord."); gDirectory->SetTitle(sdT.c_str());
      printf(" - %s",gDirectory->GetTitle());
    }
    printf("\n");
    // *************** INSTANTIATE HISTOS
    // ***** SET HISTO RANGES 
    double dX, dY, dR, ZMn, ZMx, RMn, RMx;
    int nMods = nModules[idet], modMn = moduleMns[idet], nDivs, nLayers;
    double UMn, UMx; // Outer specific: U binning
    double dRr = 0; // Endcap specific: R binning
    if      (idet==0) {             // ********** CyMBaL
      dX=dY = 600;
      double RAve = (radii[0][0]+radii[0][1])/2, deltaR = 25;
      RMn = RAve-deltaR; RMx = RAve+deltaR;
      ZMn = ZAbscissae[0][0]; ZMx = ZAbscissae[0][1];
      nDivs = nMods;
      nLayers = 1;
    }
    else if (idet==1) {             // ********** µRWELL
      dX=dY = 800;
      double RAve = radii[1][0], deltaR = 20;
      RMn = RAve-deltaR; RMx = RAve+deltaR;
      ZMn = ZAbscissae[1][0]; ZMx = ZAbscissae[1][1];
      double modHL = 2*ZHLengths[1], modHW = 2*hWidths[1][0];
      // Ur: 224 for core part + 2*16 bins on the edge to get possible stray hits
      double dRangeUr = (modHL+modHW)*sqrt(2)/2;
      double deltaUr = dRangeUr/224; UMx = dRangeUr+16*deltaUr; UMn = -UMx;
      nDivs = 2;
      nLayers = nDivs;
    }
    else if (idet==2 || idet==3) {  // ********** Endcaps
      RMn = radii[idet][0]; RMx = radii[idet][1];
      dX=dY = RMx*1.1;
      double ZAve = (ZAbscissae[idet][0]+ZAbscissae[idet][1])/2;
      double deltaZ = fabs(ZAbscissae[idet][1]-ZAbscissae[idet][0])/2;
      deltaZ *= 1.2;
      ZMn = ZAve-deltaZ; ZMx = ZAve+deltaZ;
      nDivs = 2;
      nLayers = 2;
    }
    else if (idet==4 || idet==5) {  // ********** VERTEX, Si
      double RAve = (radii[idet][0]+radii[idet][1])/2;
      double deltaR = (radii[idet][1]-radii[idet][0])/2; deltaR *= 1.1;
      RMn = RAve-deltaR; RMx = RAve+deltaR;      
      dX=dY = RMx*1.1;
      ZMn = ZAbscissae[idet][0]; ZMx = ZAbscissae[idet][1];
      if (idet==4) nDivs = 3; // 3 layers
      else         nDivs = 2;
      nLayers = nDivs;
    }
    // Z|R: 224 for core part + 2*16 bins on the edge to get possible stray hits
    if (isBarrel(idet)) {
      double deltaZ = 16*(ZMx-ZMn)/224; ZMn -= deltaZ; ZMx += deltaZ;
    }
    else {
      double deltaR = 16*(RMx-RMn)/224; RMn -= deltaR; RMx += deltaR;
    }      
    double divMn = -.5, divMx = nDivs-.5; 
    const double pi = TMath::Pi();
    char hN[] = "hthphi"; size_t lN = strlen(hN)+1;
    // ***** HISTO TITLE = DETECTOR NAME [+ STRIP COORDINATE]
    char hT[] =
      "Outer#font[32]{U};d#font[32]{Rc#delta#scale[1.2]{#varphi}}  #font[22]{(#mum)}   "; size_t lT = strlen(hT)+1;
    string dS(detectorNames[idet]);
    if (isRec && // hence one histo per strip's coord., as opposed to isRec =0
	(0x1<<idet&stripMode)) // then add strip's coord.-name as a distinctive tag
      dS += string(coordNames[iStrip-1][idet]);
    const char *dN = dS.c_str();
    // ***** LOOP ON HISTOS
    snprintf(hN,lN,"%s%s",tag,"X");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'X');
    hs.X =   new TH2D(hN,hT,256,-dX,  dX,nDivs,divMn,divMx);
    snprintf(hN,lN,"%s%s",tag,"Y");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'Y');
    hs.Y =   new TH2D(hN,hT,256,-dY,  dY,nDivs,divMn,divMx);
    snprintf(hN,lN,"%s%s",tag,"R");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'R');
    hs.R =   new TH2D(hN,hT,256,RMn, RMx,nDivs,divMn,divMx);
    if      (idet==0) { // If CyMBaL
      snprintf(hN,lN,"%s%s",tag,"Rr");
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'R');
      hs.Rr =  new TH2D(hN,hT,1024,RMn, RMx,nDivs,divMn,divMx);
    }
    else if (idet==1) { // If Outer
      snprintf(hN,lN,"%s%s",tag,"Rr");
      snprintf(hT,lT,"%s;#font[32]{%s}  #font[22]{(mm)}   ",dN,
	       "Rc#delta#scale[1.2]{#varphi}");
      hs.Rr =  new TH2D(hN,hT,1024,RMn, RMx,nDivs,divMn,divMx);
      snprintf(hN,lN,"%s%s",tag,"Ur");
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'U');
      hs.Ur =  new TH2D(hN,hT,1024,UMn,UMx,nDivs,divMn,divMx);
      snprintf(hN,lN,"%s%s",tag,"Vr");
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'V');
      hs.Vr =  new TH2D(hN,hT,1024,UMn,UMx,nDivs,divMn,divMx);
    }
    else if (idet==2 || idet==3) {
      snprintf(hN,lN,"%s%s",tag,"Rr");
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'R');
      hs.Rr =  new TH2D(hN,hT,1024,RMn,RMx,nDivs,divMn,divMx);
    }
    snprintf(hN,lN,"%s%s",tag,"Z");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'Z');
    hs.Z =   new TH2D(hN,hT,256,ZMn, ZMx,nDivs,divMn,divMx);
    snprintf(hN,lN,"%s%s",tag,"phi");
    snprintf(hT,lT,"%s;#scale[1.2]{#varphi}  #font[22]{(rad)}   ",dN);
    hs.phi = new TH2D(hN,hT,512,-pi,  pi,nDivs,divMn,divMx); 
    snprintf(hN,lN,"%s%s",tag,"phir");
    snprintf(hT,lT,"%s;#scale[1.2]{#varphi}r  #font[22]{(rad)}   ",dN);
    hs.phir = new TH2D(hN,hT,512,-pi,  pi,nDivs,divMn,divMx); 
    snprintf(hN,lN,"%s%s",tag,"th");
    snprintf(hT,lT,"%s;#theta  #font[22]{(rad)}   ",dN);
    hs.th =  new TH2D(hN,hT,256,  0,  pi,nDivs,divMn,divMx);
    snprintf(hN,lN,"%s%s",tag,"mod");
    snprintf(hT,lT,"%s;module#",dN);
    hs.mod = new TH2D(hN,hT,nMods,modMn-.5,modMn+nMods-.5,nLayers,divMn,divMx);
    // ***** eDep
    // Log binning: set a bin edge at threshold
    double eDThr = eDThresholds[idet];
    double sq4_10 = sqrt(sqrt(10.)), sq16_10 = sqrt(sqrt(sq4_10)), sq64_10 = sqrt(sqrt(sq16_10));
#define N_eDBINS 256
    double eDBin; int bin; double eDBins[N_eDBINS+1];
    for (bin = 64, eDBin = eDThr; bin>=0; bin--) {
      eDBins[bin] = eDBin; eDBin /= sq64_10;
    }
    for (bin = 64, eDBin = eDThr; bin<=N_eDBINS; bin++) {
      eDBins[bin] = eDBin; eDBin *= sq64_10;
    }
    snprintf(hN,lN,"%s%s",tag,"eDep");
    snprintf(hT,lT,"%s;eDep (keV)   ",dN);
    hs.eDep = new TH2D(hN,hT,N_eDBINS,eDBins,nDivs,divMn,divMx);

    TH2D *h2s[] = {hs.X,hs.Y,hs.Z,hs.R,hs.phi,hs.th,hs.mod,
		   hs.Rr,        // CyMBaL/Outer specific
		   hs.Ur,hs.Vr,  // Outer specific
		   hs.eDep};
    unsigned int flags[] =
      /* */       {0x3f,0x3f,0x3f,0x3f,  0x3f, 0x3f,  0x3f,
		    0x3,
		    0x2, 0x2,
		    0x3f};
    int nh2s = sizeof(h2s)/sizeof(TH2D*);
    for (int ih = 0; ih<nh2s; ih++) {
      if (!(0x1<<idet&flags[ih])) continue;
      setAxes(h2s[ih]->GetXaxis(),505,0);
      setAxes(h2s[ih]->GetYaxis(),505,2);
    }
    // ***** 2D HISTOS: X vs. Y, R vs. Z, theta vs. phi
    string sT;
    snprintf(hN,lN,"%s%s",tag,"XY");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'X');
    sT = string(hT);
    snprintf(hT,lT,";#font[32]{%c}  #font[22]{(mm)}   ",'Y');
    sT += string(hT); const char *hTXY = sT.c_str();
    hs.XY = new TH2D(hN,hTXY,256,-dX, dX,256,-dY,dY);
    if (isMPGD(idet) && isBarrel(idet)) { // Barrel MPGD: xyr = 2D in local coordinates
      double xMx, yMx; string sT = string(dN); getxyArgs(isRec,idet,xMx,yMx,sT);
      hs.xyr = new TH2D("hxyr",sT.c_str(),256,-xMx,xMx,256,-yMx,yMx);
    }
    snprintf(hN,lN,"%s%s",tag,"ZR");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'Z');
    sT = string(hT);
    snprintf(hT,lT,";#font[32]{%c}  #font[22]{(mm)}   ",'R');
    sT += string(hT); const char *hTZR = sT.c_str();
    hs.ZR = new TH2D(hN,hTZR,256,ZMn,ZMx,256,RMn,RMx);
    snprintf(hN,lN,"%s%s",tag,"thphi");
    snprintf(hT,lT,"%s;#scale[1.2]{#varphi}  #font[22]{(rad)}   ",dN);
    sT = string(hT);
    snprintf(hT,lT,";#theta  #font[22]{(rad)}");
    sT += string(hT); const char *hTthphi = sT.c_str();
    hs.thphi = new TH2D(hN,hTthphi,128,-pi,pi,128,0,pi);
    TH2D *H2s[] = {hs.XY,hs.ZR,hs.thphi,
      hs.xyr}; // MPGD specific
    int nH2s = sizeof(H2s)/sizeof(TH2D*);
    if (!isMPGD(idet) || !isBarrel(idet)) nH2s -= 1; // If !BarrelMPGD, cancel hs.xyr
    for (int ih = 0; ih<nH2s; ih++) {
      setAxes(H2s[ih]->GetXaxis(),505,2); setAxes(H2s[ih]->GetYaxis(),505,2);
    }
    if (isRec) {
      // ******************** RESIDUALS ********************
      int i01 = tag[1]=='0' ? 0 : 1;
      Resids &rs = resHs[i01][idet]; rs.dir = dSD;
      // ***** RANGES
      double dx, dy, dz, dr, du, dv; // in µm
      double dphi = 1.2; // in mrad
      if (idet==0 || idet==1) {
	dx = 608;
	double thickness = volumeThicknesses[idet], deltaR = thickness/192;
	dr = thickness/2+32*deltaR; dr *= 1000;
	dz = dx; du = dx; dv = dx;
	if (0x1<<idet&stripMode) { // ***** STRIPS: UPDATE RANGES DEPENDING ON iStrip
	  if (idet==0) { // CyMBaL
	    // 256-32-32 bins for the core part, 2x32 bins outside 
	    double modL = ZHLengths[0]*2, deltaZ = modL/192, dz_phi = modL/2+32*deltaZ;
	    // For phi, X, etc... no way to have limit at bin edge, since we have
	    // 2 distinct sensitive surface radii.
	    double dphi_Z = 2*pi/8/2; dphi_Z *= 1.20;
	    double RsurfMx = radii[0][1], dx_Z = dphi_Z*RsurfMx;
	    dz_phi *= 1000; dx_Z *= 1000; dphi_Z *= 1000; // mm -> µm, rad -> mrad
	    if (iStrip==1) dz = dz_phi;
	    else { dphi = dphi_Z; dx = dx_Z; }
	  }
	  else if (idet==1) { // Outer
	    double modL = 2*ZHLengths[1];
	    double modW = 2*hWidths[1][0];
	    // 256-16-16 bins for the core part, 2x32 bins outside 
	    double stripRange = sqrt(2)*(modW+modL)/2;
	    double deltaUV = stripRange/192, dUV = stripRange/2+32*deltaUV;
	    dUV *= 1000;
	    if (iStrip==1) dv = dUV;
	    else           du = dUV;
	    double deltaX = (modW+modL)/192; dx = modW/2+32*deltaX; dx *= 1000;
	    dz = dx;
	  }
	}
      }
      else if (idet==2 || idet==3) {
	dx = 608;
	double thickness = volumeThicknesses[idet], deltaZ = thickness/192;
	dz = thickness/2+32*deltaZ; dz *= 1000;
	dr = dx;
	dphi = dx/radii[idet][0];
      }
      else if (idet==4 || idet==5) {
	dx = 25.6;
	dz = dx;
	double thickness = volumeThicknesses[idet], deltaR = thickness/192;
	dr = thickness/2+32*deltaR; dr *= 1000;
	// phi? empirical
	if (idet==4) dphi = 0.5;
	else         dphi = 0.08;
      }
      snprintf(hN,lN,"%c%s",'d',"X");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'X');
      rs.X =   new TH1D(hN,hT,256,-dx,dx);
      snprintf(hN,lN,"%c%s",'d',"Y");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'Y');
      rs.Y =   new TH1D(hN,hT,256,-dx,dx);
      if      (idet==0) { // If CyMBaL
	snprintf(hN,lN,"%c%s",'d',"Rr"); // Rr = Reduced R
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'R');
	rs.R =   new TH1D(hN,hT,256,-dr,dr);
      }
      else if (idet==1) { // If Outer
	snprintf(hN,lN,"%c%s",'d',"Rr"); // Rr = Reduced R
	snprintf(hT,lT,"%s;d#font[32]{%s}  #font[22]{(#mum)}   ",dN,
		 "Rc#delta#scale[1.2]{#varphi}");
	rs.R =   new TH1D(hN,hT,256,-dr,dr);
	snprintf(hN,lN,"%c%s",'d',"Ur");
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'U');
	rs.Ur =   new TH1D(hN,hT,256,-du,du);
	snprintf(hN,lN,"%c%s",'d',"Vr");
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'V');
	rs.Vr =   new TH1D(hN,hT,256,-dv,dv);
      }
      else {
	snprintf(hN,lN,"%c%s",'d',"R");
	snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'R');
	rs.R =   new TH1D(hN,hT,256,-dr,dr);
      }
      snprintf(hN,lN,"%c%s",'d',"Z");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'Z');
      rs.Z =   new TH1D(hN,hT,256,-dz,dz);
      snprintf(hN,lN,"%c%s",'d',"phi");
      snprintf(hT,lT,"%s;d#scale[1.2]{#varphi}  #font[22]{(mrad)}   ", dN);
      rs.phi = new TH1D(hN,hT,512,-dphi,dphi); 
      if (idet==0) { // CyMBaL specific
	snprintf(hN,lN,"%c%s",'d',"Rphir");
	snprintf(hT,lT,"%s;Rd#scale[1.2]{#varphi}r  #font[22]{(#mum)}   ",dN);
	rs.Rphir = new TH1D(hN,hT,512,-dx,dx); 
      }
      TH1D *r1s[] =          {rs.X,rs.Y,rs.Z,rs.phi,rs.R,
			      rs.Rphir,      // CyMBaL specific
			      rs.Ur,rs.Vr};  // Outer specific
      unsigned int flags[] = {0x3f,0x3f,0x3f,  0x3f,0x3f,
			       0x1,
			       0x2,  0x2};
      int nr1s = sizeof(r1s)/sizeof(TH1D*);
      for (int ih = 0; ih<nr1s; ih++) {
	if (!(0x1<<idet&flags[ih])) continue;
	setAxes(r1s[ih]->GetXaxis(),515,3); setAxes(r1s[ih]->GetYaxis(),505,3);
      }
      /*
      // ***** 2D HISTOS: residuals vs. (x,y)
      if (idet==1) { // Only if Outer: xyr = 2D in local coordinates
	// (The idea is to illustrate how the shape of the distribution of
	// residuals along the non-measurement axis (typically asymmeric) arises
	// from the dependence of the hit count upon Z.)
	char coord = iStrip==1?'V':'U';
	snprintf(hN,lN,"d%cxy",coord);
	char sN[] = " - dU"; sprintf(sN," - d%c",coord);
	string sT = string(dN); sT += string(sN);
	double xMx, yMx; getxyArgs(idet,xMx,yMx,sT);
	rs.xyr = new TProfile2D(hN,sT.c_str(),256,-xMx,xMx,256,-yMx,yMx);
	setAxes(rs.xyr->GetXaxis(),505,3); setAxes(rs.xyr->GetYaxis(),505,3);
      }
      */
    }
  }
  dSave->cd();
}
void recoEvents::getxyArgs(bool isRec, int idet, double &xMx, double &yMx, string &sT)
{
  char hT[] =
    ";#font[32]{Z}r  #font[22]{(mm);#font[32]{R#varphi}r  #font[22]{(mm)}}   ";
  size_t lT = strlen(hT)+1;;
  if (idet==0) {
    // Since TCanvas' tend to be wider along the horizontal awis, let's have Zr
    // as abscissa and Rphir as ordinate.
    snprintf(hT,lT,";#font[32]{%s}r  #font[22]{(mm)}   ","Z");
    sT += string(hT);
    snprintf(hT,lT,";R#font[32]{%s}r  #font[22]{(mm)}   ","#varphi");
    sT += string(hT); const char *hTxy = sT.c_str();
    // Innner and outer sections do not yield the exact same Rphi range
    // Let's take the max
    yMx = radii[0][0]*hWidths[0][0];
    xMx = ZHLengths[idet];
  } else if (idet==1) {
    // Since TCanvas' tend to be wider along the horizontal awis, let's have Yr
    // as abscissa and Xr as ordinate.
    snprintf(hT,lT,";#font[32]{%s}r  #font[22]{(mm)}   ","Y");
    sT += string(hT);
    snprintf(hT,lT,";#font[32]{%s}r  #font[22]{(mm)}   ","X");
    sT += string(hT); const char *hTxy = sT.c_str();
    if (isRec) { // RecHits may lie outside detector's frame.
      double modHL = ZHLengths[idet], modHW = hWidths[idet][0];
      xMx = (modHW+modHL)/2; yMx = xMx;
    }
    else {
      yMx = hWidths[idet][0]; xMx = ZHLengths[idet];     
    }
  }
  // Add a 25µ margin to make for imprecisions
  // - Don't know what they come from: l2g? extrapolation?
  // - Turns out to be only needed form SimHits. Yet, in order to have same
  //  binning for SimHits and RecHits...
  xMx += .025; yMx += .025;
}
void recoEvents::SetNSensitiveSurfaces(int nSurfaces)
{
  if      (nSurfaces!=1 && nSurfaces!=5) {
    printf("** SetNSensitiveSurfaces: Invalid arg.(=%d); should be =1 or =5\n",
	   nSurfaces);
  }
  else nSensitiveSurfaces = nSurfaces;
}
void recoEvents::AddRequiredLayerModules(int idet, int index, unsigned long pattern)
{
  // Input required (layer,module)'s
  if (index<0 && index>3) {
    printf("** AddResuiredLayerModules: Invalid <index> arg.(=%d); should be w/in [0,3]\n",
	   index);
  }
  else {
    requireModules = 1;
    requiredLayerModules[idet].patterns[index] = pattern;
    printf("requiredLayerModules[index(i.e. quasi-layer)]:");
    for (int idx = 0; idx<4; idx++) {
      if (requiredLayerModules[idet].patterns[idx])
	printf(" [%d]: pattern = 0x%016lx",
	       idx,requiredLayerModules[idet].patterns[idx]);
    }
    printf("\n");
  }
}
bool recoEvents::isMPGD(int idet)
{
  return 0x1<<idet&MPGDs;
}
bool recoEvents::isBarrel(int idet)
{
  return 0x1<<idet&Barrels;
}

Bool_t recoEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void recoEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t recoEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void recoEvents::DrawphithZR(int iSimRec, // 0: sim, 1: rec, 2: rec 2nd coord of STRIP
			     unsigned int detectorPattern,
			     unsigned int histoPattern, bool decompose,
			     TCanvas *cPrv, int ipad, int col) // Superimpose on pre-existing TCanvas cPrv
{
  // ***** PARSE <iSimRec> ARG.
  if (iSimRec<0 || 2<iSimRec) {
    printf("** DrawphithZR: Invalid <iSimRec> arg.(=%d): neither 0(=sim) nor 1(=rec) nor 2(=2nd coord of STRIP)\n",iSimRec); return;
  }
  if (!reconstruction && iSimRec) {
    printf("** DrawphithZR: Invalid <iSimRec> arg.(=%d): recoEvents is simulationOnly\n",iSimRec); return;
  }
  TDirectory *dSave = gDirectory;
  unsigned int pat = processedDetectors&detectorPattern;
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&pat)) continue;
    Histos *Hs; switch (iSimRec) {
    case 0: Hs = simHs;    break;
    case 1: Hs = recHs[0]; break;
    case 2: if (!(0x1<<idet&stripMode)) {
	printf("** DrawphithZR: requesting <iSimRec> = 2 for idet = 0x%x (\"%s\") which is not among detectors w/ strips(=0x%x)\n",
	       0x1<<idet,detectorNames[idet],stripMode); continue;
      }
      /* */ Hs = recHs[1]; break;
    default: printf("** DrawphithZR: Invalid <iSimRec> (=%d, not in [0,2])\n",
		    iSimRec); return;
    }
    Histos &hs = Hs[idet]; hs.dir->cd();
    TCanvas *cEvents; int pad; if (!cPrv) {
      // ***** TCANVAS
      string cS = string("c")+string(detectorNames[idet]);
      char cTag[] = "0";    // tag to cope w/ several distinct recoEvents object
      int iObj = iObjCreated; if (iObj) { *cTag += iObj%10; cS += string(cTag); }
      if (iSimRec)
	cS += iSimRec==2 ? string("Rec2") : string("Rec");
      const char *cN = cS.c_str();
      cEvents = new TCanvas(cN,cN);
      cEvents->Divide(2,2);
      pad = 1;
      //gStyle->SetOptStat(10);
    }
    else {
      cEvents = cPrv;
      pad = ipad;
    }
    // ***** LOOP ON (selected subset of) HISTOS
    TH2D *h2s[5] = { hs.phi,hs.th,hs.Z,hs.R,hs.eDep};
    unsigned int flags[5] = // 0x1: OptStat = 1110, 0x2: Logx, 0x4: Logy
      /* */        {      0,    0, 0x8, 0x4,    0x3};
    if (idet<2) h2s[3] = hs.Rr;
    int ih, jh; for (ih=jh = 0; ih<5 && jh<4; ih++) {
      if (!(0x1<<ih&histoPattern)) continue;
      cEvents->cd(pad++);
      TH2D *h2 = h2s[ih];
      TH1D *hproj = h2->ProjectionX(); hproj->Draw();
      if      ((flags[ih]&0x4) &&  isBarrel(idet)) {
	hproj->SetMinimum(.5); gPad->SetLogy();
      }
      else if ((flags[ih]&0x8) && !isBarrel(idet)) {
	hproj->SetMinimum(.5); gPad->SetLogy();
      }
      else
	hproj->SetMinimum(0); // Useful for hs.phi
      if (flags[ih]&0x2) gPad->SetLogx();
      TAxis *ay = hproj->GetYaxis();
      ay->SetNdivisions(505); ay->SetLabelFont(62);
      ay->SetLabelSize(.055); ay->SetLabelOffset(.006);
      ay->SetMaxDigits(2);
      int optStat = (flags[ih]&0x1) ? 1110 : 10;
      SetPaveText(hproj,0,optStat);
      if (decompose) {
	char tag[] = "_3";
	if      (idet==0 && (ih==0 || ih==3)) { // CyMBaL phi, R
	  int cols[2] = {bleu,rouge};
	  for (int module = 0; module<8; module++) {
	    snprintf(tag,3,"_%d",module); string hS = h2->GetName(); hS += tag;
	    TH1D *hsum = h2->ProjectionX(hS.c_str(),module+1,module+1);
	    hsum->SetLineColor(cols[module%2]);
	    for (int sector = 1; sector<4; sector++) {
	      int modvle = 8*sector+module;
	      TH1D *helem = h2->ProjectionX("temp",modvle+1,modvle+1);
	      hsum->Add(helem);
	    }
	    hsum->Draw("same");
	    hsum->SetLineColor(cols[module%2]);
	  }
	}
	else if (idet==0) {                     // CyMBaL theta,Z
	  int cols[4] = {violet,vert,orange,marron};
	  for (int module = 0; module<32; module += 8) {
	    snprintf(tag,3,"_%d",module/8); string hS = h2->GetName(); hS += tag;
	    TH1D *h1 = h2->ProjectionX(hS.c_str(),module+1,module+8);
	    h1->Draw("same");
	    h1->SetLineColor(cols[module/8]);
	  }
	}
	else {          // µRWELL, Vertex, Si
	  int col1s[2] = {vert,orange};
	  int col2s[3] = {bleu,vert,rouge};
	  int *cols = idet==1 ? col1s : col2s;
	  int nDivs = h2->GetNbinsY();
	  for (int div = 0; div<nDivs; div++) {
	    snprintf(tag,3,"_%d",div); string hS = h2->GetName(); hS += tag;
	    TH1D *h1 = h2->ProjectionX(hS.c_str(),div+1,div+1);
	    h1->Draw("same");
	    h1->SetLineColor(cols[div%nDivs]);
	  }
	}
      }
      jh++;
    }
  }
  dSave->cd();
}
void recoEvents::DrawModules(int iSimRec, unsigned int detectorPattern, bool decompose)
{
  // ********** DRAW module# FOR EACH DETECTOR IN detectorPattern
  // ***** PARSE ARG.S
  if (iSimRec<0 || 2<iSimRec) {
    printf("** DrawModules: Invalid <iSimRec> arg.(=%d): neither 0(=sim) nor 1(=rec) nor 2(=2nd coord of STRIP)\n",iSimRec); return;
  }
  if (!reconstruction && iSimRec) {
    printf("** DrawphithZR: Invalid <iSimRec> arg.(=%d): recoEvents is simulationOnly\n",iSimRec); return;
  }
  unsigned int pat = processedDetectors&detectorPattern;
  if (!pat) {
    printf("** DrawModules: None of the processed detectors(=0x%x) in arg. <detectorPattern>(=0x%x)\n",
	   processedDetectors,detectorPattern);
    return;
  }
  TDirectory *dSave = gDirectory;
  for (int idet = 0; idet<N_DETs; idet++) {
    // ***** LOOP ON DETECTORS
    if (!(0x1<<idet&pat)) continue;
    Histos *Hs; switch (iSimRec) {
    case 0: Hs = simHs; break;
    case 1: Hs = recHs[0]; break;
    case 2: if (!(0x1<<idet&stripMode)) {
	printf("** DrawModules: requesting <iSimRec> = 2 for idet = 0x%x (\"%s\") which is not among detectors w/ strips(=0x%x)\n",
	       0x1<<idet,detectorNames[idet],stripMode); continue;
      }
      /* */ Hs = recHs[1]; break;
    default: printf("** DrawModules: Invalid <iSimRec> (=%d, not in [0,2])\n",
		    iSimRec); return;
    }
    Histos &hs = Hs[idet]; hs.dir->cd();
    TH2D *h2 = hs.mod;
    TH1D *hproj = h2->ProjectionX(); hproj->Draw();
    hproj->SetMinimum(0);
    // Port TH2D's title to proj. Else title displayed is always that of the first module histo drawn!? As evidenced by instruction below (commented out).
    //printf("<%s> <%s>\n",hproj->GetName(),hproj->GetTitle());
    hproj->SetTitle(h2->GetTitle());
    TAxis *ay = hproj->GetYaxis();
    ay->SetNdivisions(505);
    ay->SetLabelSize(.05); ay->SetLabelOffset(.006);
    ay->SetMaxDigits(2);
    int nDivs = h2->GetNbinsY();
    if (nDivs>1 && decompose) {
      char tag[] = "_3";
      int col1s[2] = {vert,orange};
      int col2s[3] = {bleu,vert,rouge};
      int *cols = idet==1 ? col1s : col2s;
      for (int div = 0; div<nDivs; div++) {
	snprintf(tag,3,"_%d",div); string hS = h2->GetName(); hS += tag;
	TH1D *h1 = h2->ProjectionX(hS.c_str(),div+1,div+1);
	h1->Draw("same");
	h1->SetLineColor(cols[div%nDivs]);
      }  
    }
  }
  dSave->cd();
}
void recoEvents::DrawResiduals(int iRec, // 1: rec, 2: 2nd coord of STRIPS phi/Z U/V
			       unsigned int detectorPattern,
			       unsigned int histoPattern,
			       TCanvas *cPrv, int ipad, int col) // Superimpose on pre-existing TCanvas cPrv
{
  // ********** DRAW (subset of) RESIDUALS FOR EACH DETECTOR IN detectorPattern
  // ***** PARSE ARG.S
  if (iRec<1 || 2<iRec) {
    printf("** DrawResiduals: Invalid <iRec> arg.(=%d): neither 1(=1st coord.) nor 2(=2nd coord.)\n",iRec); return;
  }
  if (!reconstruction) {
    printf("** DrawphithZR: recoEvents is simulationOnly\n"); return;
  }
  int strip = iRec-1;
  unsigned int pat = processedDetectors&detectorPattern;
  if (!pat) {
    printf("** DrawResiduals: None of the processed detectors(=0x%x) in arg. <detectorPattern>(=0x%x)\n",
	   processedDetectors,detectorPattern);
    return;
  }
  // ***** CONTEXT
  //gStyle->SetOptStat(1110);
  string prvFormat = string(gStyle->GetStatFormat());
  TDirectory *dSave = gDirectory;
  // ***** LOOP ON DETECTORS
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&pat)) continue;
    if (strip && !(0x1<<idet&stripMode)) { // 2nd coordinate: skip if not STRIP
      printf("** DrawResiduals: requesting <iRec> = 2 for idet = 0x%x (\"%s\") which is not among detectors w/ strips(=0x%x)\n",
	       0x1<<idet,detectorNames[idet],stripMode); continue;
    }
    TCanvas *cResids; int pad; if (!cPrv) {
      // ***** CANVAS
      string cS = string("c")+string(detectorNames[idet]);
      char cTag[] = "0";    // tag to cope w/ several distinct recoEvents object
      int iObj = iObjCreated; if (iObj) { *cTag += iObj%10; cS += string(cTag); }
      cS += strip ? string("Res2") : string("Res");
      const char *cN = cS.c_str();
      cResids = new TCanvas(cN,cN);
      cResids->Divide(2,2);
      pad = 1;
    }
    else {
      cResids = cPrv;
      pad = ipad;
    }
    gStyle->SetStatFormat("6.3g");
    // ***** LOOP ON RESIDUALS (in subset)
    Resids &rs = resHs[strip][idet]; rs.dir->cd();
    TH1D *r1s[4] = {rs.X, rs.Z, rs.phi, rs.R};
    int printIntegral = -1; // Print Integral intead of Entries, to check for Under/Overflow
    if      (idet==0) { // CyMBaL: precision is on 'Rphir' or 'Z' alternatively,
      // depending upon strip's coord. 'Rr': let's check it's a Dirac.
      r1s[2] = rs.Rphir;
      printIntegral = iRec==1 ? 2 : 1; // Integral on Rphir/Z
    }
    else if (idet==1) { // Outer:  precision is on 'Ur' or 'Vr' alternatively,
      // depending upon strip's coord. 'Rr' = "Rcos(dphi)": check it's a Dirac. 
      r1s[1] = rs.Ur; r1s[2] = rs.Vr;
      printIntegral = iRec==1 ? 1 : 2; // Integral on Ur/Vr
    }
    for (int ih = 0; ih<4; ih++) {
      if (!(0x1<<ih&histoPattern)) continue;
      TH1D *h1 = r1s[ih];
      if (col>=0) h1->SetLineColor(col);
      cResids->cd(pad++);
      bool superImpose = false;
      if (cPrv) {
	TList *l = gPad->GetListOfPrimitives(); TObject *o = l->First();
	while (o) {
	  if (o->IsA()->InheritsFrom(TH1::Class())) {
	    superImpose = true; break;
	  }
	  o = l->After(o);
	}
      }
      int optStat = ih==printIntegral ? 1001100 : 1110;
      if (superImpose) {
	h1->Draw("sames"); SetPaveText(h1,1,optStat); }
      else             {
	h1->Draw();        SetPaveText(h1,0,optStat);
	if      (r1s[ih]==rs.R &&  isBarrel(idet) ||
		 r1s[ih]==rs.Z && !isBarrel(idet)) {
	  h1->SetMinimum(.5); gPad->SetLogy();
	}
      }
      //  Ndivisions policy: 505 -> 510, when Xaxis labels come close to the edge and collide w/ the Yaxis 0 label
      TAxis *ax = h1->GetXaxis();
      double xLabel = ax->GetXmax()/pow(10,int(log10(ax->GetXmax())));
      if (fabs(xLabel-5)<.3)
	ax->SetNdivisions(515);
    }
  }
  gStyle->SetStatFormat(prvFormat.c_str());
  dSave->cd();
}
#include "TPaveText.h"
#include "TPaveStats.h"
void SetPaveText(TH1 *h, int mode, int opt)
{
  TPad *pad = (TPad*)gPad; pad->Update();

  int col = h->GetLineColor();

  string format = string("6.3g");
  int optStat = opt ? opt : 10;
  double dY, X1; if (optStat==10)      { dY = .10; X1 = .70; }
  else           if (optStat==1110)    { dY = .22; X1 = .55; }
  else           if (optStat==1001100) {
    /* */                                dY = .22; X1 = .55;
    format = string("6.0f");
  }

  TPaveStats *st; if ((st = (TPaveStats*)h->GetFunction("stats"))) {
    st->SetTextFont(62);
    st->SetX2NDC(.995); st->SetX1NDC(.60);
    st->SetY2NDC(.995-mode*dY); st->SetY1NDC(.995-(mode+1)*dY);
    st->SetOptStat(optStat);
    st->SetStatFormat(format.c_str());
    st->SetTextColor(col);
    st->Draw();
  }

  TPaveText *tit; if ((tit = (TPaveText*)pad->GetPrimitive("title"))) {
    if (mode==0) { // First histo in pad
      tit->SetTextFont(62);
      tit->SetX1NDC(.25); tit->SetX2NDC(.45);
      tit->SetY1NDC(.92); tit->SetY2NDC(.99);
      tit->Draw();
    }
  }
}
void recoEvents::SetMinima(double min)
{
  // Apply "SetMinimum" on all histos in "gPad"
  TList *l = gPad->GetListOfPrimitives();
  TObject *o = l->First();
  while (o) {
    if (o->IsA()->InheritsFrom("TH1D")) {
      TH1D *h = (TH1D*)o;
      printf(" * SetMinima: %s->SetMinimum(%f);\n",h->GetName(),min);
      h->SetMinimum(min);
    }
    o = l->After(o);
  }
}
#endif // #ifdef recoEvents_cxx
/*
    // HOW TO WRITE TO A PDF FILE
  string pdfName("crater_lake.e-.pdf");
  string pdfOpen(pdfName);  pdfOpen +=  string("(");
  string pdfClose(pdfName); pdfClose += string(")");
  cCyMBaL->Print(pdfOpen.c_str(),"Title:e- CyMBaL");
  //eAGM9_203->Print(pdfName.c_str(),"Title:GM09");
  c2->Print(pdfClose.c_str(),"Title:e- All");
*/
/*
  // All in one (2,2) canvas
  TCanvas *c = new TCanvas("c3"); c->Divide(2,2);
  char hN[] = "hth_9"; int i; TH1D *h;
  gDirectory->cd("/dSim/dCyMBaL");
  c->cd(1); hth_px->Draw();
  for (i=0;i<4;i++) { snprintf(hN,6,"hth_%d",i); h=(TH1D*)gDirectory->Get(hN); h->SetStats(0); h->Draw("same"); }
  gDirectory->cd("/dSim/dOuter");
  c->cd(3); hth_px->Draw();
  for (i=0;i<2;i++) { snprintf(hN,6,"hth_%d",i); h=(TH1D*)gDirectory->Get(hN); h->SetStats(0); h->Draw("same"); }
  gDirectory->cd("/dSim/dSi");
  c->cd(2); hth_px->Draw();
  for (i=0;i<2;i++) { snprintf(hN,6,"hth_%d",i); h=(TH1D*)gDirectory->Get(hN); h->SetStats(0); h->Draw("same"); }
  gDirectory->cd("/dSim/dVertex");
  c->cd(4); hth_px->Draw();
  for (i=0;i<3;i++) { snprintf(hN,6,"hth_%d",i); h=(TH1D*)gDirectory->Get(hN); h->SetStats(0); h->Draw("same"); }  
 */
/*
  // RESIDUALS W/O W/ QUALITY REQUIRED
  TCanvas *c;
  .L recoEvents.so
  TFile *_file0 = TFile::Open("podio_output.root"); TTree *events = (TTree*)gDirectory->Get("events");
  recoEvents ana(events,0xf);
  ana.Loop();
  recoEvents ana2(events,0xf);
  ana2.requireQuality = 1; // REQUIRE quality IN SimTrackerHits
  ana2.Loop();
  // ***** CyMBaL
  ana.DrawResiduals(0x1,0,2);
  c = (TCanvas*)gROOT->FindObject("cCyMBaLRes");
  if (c) ana2.DrawResiduals(0x1,c);
  // ***** Outer
  ana.DrawResiduals(0x2,0,2);
  c = (TCanvas*)gROOT->FindObject("cOuterRes");
  if (c) ana2.DrawResiduals(0x2,c);
  // ***** Vertex
  ana.DrawResiduals(0x4,0,2);
  c = (TCanvas*)gROOT->FindObject("cVertexRes");
  if (c) ana2.DrawResiduals(0x4,c);
  // ***** Si
  ana.DrawResiduals(0x8,0,2);
  c = (TCanvas*)gROOT->FindObject("cSiRes");
  if (c) ana2.DrawResiduals(0x8,c);
  // HOW TO WRITE TO A PDF FILE
  string pdfName("SimRec.pdf");
  string pdfOpen(pdfName);  pdfOpen +=  string("(");
  string pdfClose(pdfName); pdfClose += string(")");
  cCyMBaLRes->Print(pdfOpen.c_str());
  cOuterRes->Print(pdfName.c_str());
  cVertexRes->Print(pdfName.c_str());
  cSiRes->Print(pdfClose.c_str());
  // WRITE TO ROOT FILE
  TDirectory *dSave = gDirectory
  TFile *fout = TFile::Open("SimRes.root","CREATE");
  cCyMBaLRes->Write();
  cOuterRes->Write();
  cVertexRes->Write();
  cSiRes->Write();
  fout->Close();
  dSave->cd();

  // Rr
  gStyle->SetOptStat(10);
  TCanvas *cRr = new TCanvas("cRr");
  cRr->Divide(2,1);
  cRr->cd(1);
  TH1D *h1 = ana.simHs[0].Rr->ProjectionX("Rr");
  h1->Draw(); SetPaveText(h1,0);
  h1->SetLineColor(2);
  TH1D *h1q = ana2.simHs[0].Rr->ProjectionX("Rrq");
  h1q->Draw("sames"); SetPaveText(h1q,1);
  cRr->cd(2);
  ana.simHs[0].XYr[0]->Draw(); SetPaveText(ana.simHs[0].XYr[0],0)
  ana.simHs[0].XYr[0]->Draw("box");
  ana.simHs[0].XYr[0]->SetLineColor(kGray);
  ana.simHs[0].XYr[1]->Draw("boxsames"); SetPaveText(ana.simHs[0].XYr[1],1)
  ana.simHs[0].XYr[2]->Draw("boxsames");
  ana.simHs[0].XYr[3]->Draw("boxsames");
  ana.simHs[0].XYr[3]->SetLineColor(6);

*/
/*
  new TCanvas("cCyMBaLRes"); new TCanvas("cSiRes"); new TCanvas("cOuterRes"); new TCanvas("cVertexRes");
  string pdfName("SimRec.pdf");
  string pdfOpen(pdfName);  pdfOpen +=  string("(");
  string pdfClose(pdfName); pdfClose += string(")");
  cCyMBaLRes->Print(pdfOpen.c_str(),"Title:CyMBaL");
  cOuterRes->Print(pdfName.c_str(),"Title:Outer");
  cVertexRes->Print(pdfName.c_str(),"Title:Vertex");
  cSiRes->Print(pdfClose.c_str(),"Title:Si");
*/
/*
  // PDF w/ EmbededFonts
  string pdfName("Sim.CyMBaL.PR#2177.pdf");
  cCyMBaL->Print(pdfName.c_str(),"EmbedFonts");
  string pdfName("Res.CyMBaL.PR#2177.pdf");
  cCyMBaLRes->Print(pdfName.c_str(),"EmbedFonts");
  string pdfName("Sim.Outer.PR#2177.pdf");
  cOuter->Print(pdfName.c_str(),"EmbedFonts");
  string pdfName("Res.Outer.PR#2177.pdf");
  cOuterRes->Print(pdfName.c_str(),"EmbedFonts");
  string pdfName("AllSim.CyMBaL.PR#2177.pdf");
  cCyMBaL1->Print(pdfName.c_str(),"EmbedFonts");
  string pdfName("AllRes.CyMBaL.PR#2177.pdf");
  cCyMBaL1Res->Print(pdfName.c_str(),"EmbedFonts");
*/
