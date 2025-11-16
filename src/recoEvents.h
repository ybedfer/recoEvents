//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  5 23:53:48 2024 by ROOT version 6.30/02
// from TChain events/
//////////////////////////////////////////////////////////

/*
// USAGE
const char tag[7] = "struv", nEvtsTag[3] = "20"; 
TChain *events = new TChain("events"); char fName[] = "podio.sensor.1234.20,root"; size_t lN = strlen(fName)+1; int seeds[] = {1234,4567,8910,1112}; int nSeeds = sizeof(seeds)/sizeof(int);
for (int i = 0; i<nSeeds; i++) { snprintf(fName,lN,"podio.%s.%s.%4d.root",tag,nEvtsTag,seeds[i]); events->Add(fName); }
.L ../recoEvents/install/librecoEvents.so
recoEvents ana(events,0xf); // Instantiate for all of (0x1:CyMBaL,0x2:Outer,0x4:Vertex,0x8:Si)
ana.select = new TTreeFormula("select", "@MCParticles.size()==1", events);// Add rejection cut
ana.Loop();
ana.DrawphithZR(0,0xf,true); // Draw some simHits, w/ if true, colour highlighting of module type.
ana.DrawphithZR(1,0xf,true); // Draw some recHits (1st coord), w/ if true, colour highlighting of module type.
ana.DrawResiduals(0,0x1); // Draw some residuals for phi (iRec=1) of CyMBaL (=0x1)

ana.recHs[0][0].ZR[0]; h2->SetTitle("CyMBaL Rec"); h2->Draw();  SetPaveText(h2,0)
*/

#ifndef recoEvents_h
#define recoEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
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
#include "edm4hep/TrackerHitData.h"

using namespace std;
using namespace edm4hep;

typedef struct{ TDirectory *dir; TH2D *X, *Y, *Z, *RA, *R, *phi, *phir, *th, *mod, *thphi, *XY, *ZR, *Rr, *XYr[4], *Ur, *Vr; }
  Histos;
typedef struct{ TDirectory *dir; TH1D *X, *Y, *Z, *R, *Rr, *D, *phi, *phir, *Ur, *Vr; }
  Resids;

void SetPaveText(TH1 *h, int mode = 0);

class recoEvents: public TNamed {
public :
   recoEvents(TTree *tree,
	      // Detector pattern: 0x1=CyMBaL, 0x2=Outer, 0x4=Vertex, 0x8=Si
	      // mode: Processed detectors
	      // hasStrips: If not set, detector has pixels.
	      unsigned int mode = 0x1, unsigned int hasStrips = 0x3);
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
   int verbose;

   int getDetHit(int idet, int ih, double &X, double &Y, double &Z,
		 unsigned int &module, unsigned int &div);

   // ***** DETECTOR NAMES
#define N_DETs 4
   const char *detectorNames[N_DETs];
   // STRIP SEGMENTATION
   const char *coordNames[2][N_DETs];
   // Depending upon setup, MPGDs can have more than one sensitive surfaces.
   // Can be 1, 2 or 5.
   int nSensitiveSurfaces;

   // ***** HISTOS
   TTreeFormula *select; // Provides for specifying a rejection cut.
   // Requirements for filling "simHs" and residuals
   int requirePDG;     // Associated MCParticle 
   int requireQuality; // "SimTrackerHit::quality"
   int requireModule; // If >=0, select module "requireModule" (D=-1=All)
   // Cut on #SimHits per detector (and, not yet implemented, per MC particle)
   // >0: #SimHits=="requireNHits", <0: #SimHits>="requireNHits"
   int requireNHits;
   void BookHistos(Histos *Hs, const char *tag);
   unsigned int getStatus(int idet, int ih, int quality);
   void fillHit(int iSimRec, int idet,
		double X, double Y, double Z, unsigned long cellID);
   void fillResids(int idet,
		   const Vector3f &pos, const Vector3d &psim, unsigned long cellID);
   void parseCellID(int idet, unsigned long ID,
		    unsigned int &module, unsigned int &div, unsigned int &strip);
   bool parseStrip(int idet, int simOrRec, unsigned int &strip);
   void printHit(int idet,
		 double X, double Y, double Z, unsigned long cellID);
   static constexpr int violet = kMagenta+2, bleu = kBlue+1, vert = kGreen+2;
   static constexpr int jaune = kOrange+0, orange = kOrange+1, rouge  = kRed+0;
   static constexpr int marron = 50, gris   = 11;
   void DrawphithZR(int iSimRec = 0 /* 0: sim, 1: rec, 2: rec 2nd coord of STRIP */,
		    unsigned int detectorPattern = 0xf, bool decompose = false);
   void DrawModules(int iSimRec = 0, unsigned int detectorPattern = 0x1, bool decompose = false);
   void DrawResiduals(int iRec = 1,   unsigned int detectorPattern = 0x1, TCanvas *cPrv = 0, int col = -1);
   void SetMinima(double min);
   TDirectory *dSim, *dSimDets[N_DETs];
   TDirectory *dRec, *dRecDets[N_DETs];
   Histos simHs[N_DETs], recHs[2][N_DETs];
   Resids resHs[2][N_DETs];
   TH1D *hMult;

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
   vector<MCParticleData> *mcParticles;
   TBranch *MCParticles;
   vector<EventHeaderData> *eventHeader;
   TBranch *EventHeader;

   virtual void     Loop(int nEvents = 0);

 private:
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   ClassDef(recoEvents,1); // Must be the last item before the closing '};'
};

#endif

#ifdef recoEvents_cxx
recoEvents::recoEvents(TTree *tree, unsigned int detectors, unsigned int hasStrips) : fChain(0) 
{
  // Init global settings
  processedDetectors = detectors; verbose = 0; select = 0;
  requirePDG = 0; // Default: do not require any ID
  requireQuality = 0; // Default: do not require "SimTrackerHit::quality"
  requireModule = -1;  // Default: all staves
  requireNHits = 0; // Default: no requirement
  // Check that arg. "hasStrips" fits the pattern of MPGDs
  unsigned int MPGDs = 0x3; // CyMBaL and Outer
  if ((hasStrips&MPGDs)!=hasStrips) {
    printf("** recoEvents: Invalid arg. <hasStrips>(=0x%x): does not fit pattern of MPGDs(=0x%x)! => Aborting...\n",
	   hasStrips,MPGDs);
    return;
  }
  stripMode = hasStrips;
  nSensitiveSurfaces = 1;
  iObjCreated = nObjCreated++;
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
   if (!tree) return;
   fChain = (TChain*)tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1); // Leads to run time crash as soon as one act on mpgdHs

   // ********** BRANCHES
   // ***** SIMHITS
   for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "MPGDBarrel","OuterMPGDBarrel","VertexBarrel","SiBarrel"};
     const char *branchName = branchNames[idet];
     string name(branchName); name += string("Hits");
     fChain->SetBranchAddress(name.c_str(),&hits[idet],&simBranches[idet]);
   }
   // ***** MC Particle
   for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "_MPGDBarrel","_OuterMPGDBarrel","_VertexBarrel","_SiBarrel"};
     const char *branchName = branchNames[idet];
     string name(branchName); name += string("Hits_MCParticle");
     fChain->SetBranchAddress(name.c_str(),&amcs[idet],&amcBranches[idet]);
   }
   // ***** REC BRANCHES
   // Determine whether available for all requested detectors
   unsigned int recBranchPat = 0; for (int idet = 0; idet<N_DETs; idet++) {
     if (!(0x1<<idet&processedDetectors)) continue;
     const char *branchNames[N_DETs] = {
       "MPGDBarrel","OuterMPGDBarrel","SiBarrelVertex","SiBarrelTracker"};
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
	 "MPGDBarrel","OuterMPGDBarrel","SiBarrelVertex","SiBarrelTracker"};
       const char *branchName = branchNames[idet];
       string name(branchName); name += string("RecHits");
       fChain->SetBranchAddress(name.c_str(),&recs[idet],&recBranches[idet]);
     }
     // ***** RECHITS <-> RAWHIT ASSOCIATION
     for (int idet = 0; idet<N_DETs; idet++) {
       const char *branchNames[N_DETs] = {
	 "_MPGDBarrel","_OuterMPGDBarrel","_SiBarrelVertex","_SiBarrel"};
       const char *branchName = branchNames[idet];
       string namr(branchName); namr += string("RawHitAssociations_rawHit");
       string nams(branchName); nams += string("RawHitAssociations_simHit");
       fChain->SetBranchAddress(namr.c_str(),&arhs[idet],&arhBranches[idet]);
       fChain->SetBranchAddress(nams.c_str(),&ashs[idet],&ashBranches[idet]);
     }
   }
   fChain->SetBranchAddress("MCParticles",&mcParticles,&MCParticles);
   fChain->SetBranchAddress("EventHeader",&eventHeader,&EventHeader);

   Notify();

   // ***** DETECTOR NAMES
   const char *detNs[N_DETs] = {"CyMBaL","Outer","Vertex","Si"};
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
// *************** HISTOS
void recoEvents::BookHistos(Histos *Hs, const char* tag)
{
  using namespace ROOT;
  using namespace std;

  // ***** PARSE <tag> ARG.
  // Check internal consistency
  if (strcmp(tag,"s") && strcmp(tag,"r0") && strcmp(tag,"r1")) {
    printf("** BookHistos: Inconsistency: Unvalid <tag> arg.(=\"%s\"). Aborting...\n",tag);
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
    hMult = new TH1D("hMult","@MCparticles.size()",512,0,512);

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
    double dX, dY, RAve, dR, ZMn, ZMx; int nMods, modMn, nDivs, nLayers;
    double dUr = 0;
    if      (idet==0) { // ********** CyMBaL
      dX=dY = 600;
      RAve = 565; dR = 25;
      // ***** Z RANGE: Let extrema of sensitive area fall on bin edge
      // From: epic/compact/tracking/definitions_craterlake.xml
      //<constant name="InnerMPGDBarrel_zmin"            value="105*cm"/> <comment> negative z </comment>
      //<constant name="InnerMPGDBarrel_zmax"            value="143*cm"/> <comment> positive z </comment>
      // From: epic/compact/tracking/mpgd_barrel.xml
      //<constant name="MMOutwardFrameWidth"                    value="5.0*cm"/>
      ZMn = -1000; ZMx = 1380;
      nMods=nDivs = 32; modMn = 0; // 4 sectors * 8 staves, numbering from 0
      nLayers = 1;
    }
    else if (idet==1) { // ********** µRWELL
      dX=dY = 800;
      RAve = 735; dR = 20;
      // <constant name="MPGDOuterBarrelModule_zmin1"     value="164.5*cm"/>
      // <constant name="MPGDOuterBarrelModule_zmin2"     value="174.5*cm"/>
      ZMn = -1645; ZMx = 1745;
      //<constant name="MPGDOuterBarrelModule_length"                  value="0.5*(MPGDOuterBarrelModule_zmin1 + MPGDOuterBarrelModule_zmin2 + MPGDOuterBarrelModule_zoverlap)"/>
      //<constant name="MPGDOuterBarrelModule_zoverlap"                value="0*MPGDOuterBarrelFrame_width"/>
      double zmin1 = 1640.5, zmin2 = 1740.5, modL = (zmin1+zmin2)/2;
      // <constant name="MPGDOuterBarrelModule_width"                   value="36.0*cm"/>
      double modW = 360;
      // UR: 224 for core part + 2*16 bins on the edge to get possible stray hits
      double dRangeUr = (modL+modW)*sqrt(2)/2;
      double deltaUr = dRangeUr/224; dUr = dRangeUr+16*deltaUr;
      nMods = 24; nDivs = 2; modMn = 0; // Numbering from 1
      nLayers = nDivs;
    }
    else if (idet==2) { // ********** VERTEX
      dX=dY = 140;
      // <constant name="VertexBarrel_rmin"               value="3.6*cm"/>
      // <constant name="VertexBarrel_rmax"               value="12.6*cm"/>
      RAve = 80; dR = 50;
      // <constant name="VertexBarrel_length"             value="270.0*mm"/>
      ZMn = -135; ZMx = 135;
      nMods = 128; nDivs = 3; modMn = 1; // 3 layers, numbering from 1
      nLayers = nDivs;
    }
    else if (idet==3) { // ********** Si
      dX=dY = 500;
      // <constant name="SiBarrel1_rmin"                  value="27.0*cm"/>
      // <constant name="SiBarrel2_rmin"                  value="42.0*cm"/>
      RAve = 350; dR = 100;
      // <comment> projective cone at 45 degree </comment>
      // <constant name="SiBarrel_angle"                  value="TrackerPrimaryAngle"/>
      ZMn = -500; ZMx = 500;
      nMods = 68; nDivs = 2; modMn = 1; // Numbering from 1
      nLayers = nDivs;
    }
    // Z: 224 for core part + 2*16 bins on the edge to get possible stray hits
    double dZ = 16*(ZMx-ZMn)/224; ZMn -= dZ; ZMx += dZ;
    double RAMx = 2*RAve, RMn = RAve-dR, RMx = RAve+dR;
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
    snprintf(hN,lN,"%s%s",tag,"RA");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'R');
    hs.RA =  new TH2D(hN,hT,256,  0,RAMx,nDivs,divMn,divMx);
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
      hs.Ur =  new TH2D(hN,hT,1024,-dUr,dUr,nDivs,divMn,divMx);
      snprintf(hN,lN,"%s%s",tag,"Vr");
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'V');
      hs.Vr =  new TH2D(hN,hT,1024,-dUr,dUr,nDivs,divMn,divMx);
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
    TH2D *h2s[] = {hs.X,hs.Y,hs.RA,hs.R,hs.Z,hs.phi,hs.th,hs.mod,
		   hs.Rr,        // CyMBaL/Outer specific
		   hs.Ur,hs.Vr}; // Outer specific
    unsigned int flags[] =
      /* */       { 0xf, 0xf,  0xf, 0xf, 0xf,   0xf,  0xf,   0xf,
		    0x3,
		    0x2, 0x2};
    int nh2s = sizeof(h2s)/sizeof(TH2D*);
    for (int ih = 0; ih<nh2s; ih++) {
      if (!(0x1<<idet&flags[ih])) continue;
      TAxis *ax = h2s[ih]->GetXaxis();
      ax->SetNdivisions(505);
      ax->SetLabelSize(.05); ax->SetLabelOffset(.006);
      ax->SetTitleSize(.06); ax->SetTitleOffset(.8);
      TAxis *ay = h2s[ih]->GetYaxis();
      ay->SetNdivisions(505);
      ay->SetLabelSize(.05); ay->SetLabelOffset(.006);
      ay->SetMaxDigits(2);
    }
    // ***** 2D HISTOS: X vs. Y, R vs. Z, theta vs. phi
    string sT;
    snprintf(hN,lN,"%s%s",tag,"XY");
    snprintf(hT,lT,"%s;#font[32]{%c}  #font[22]{(mm)}   ",dN,'X');
    sT = string(hT);
    snprintf(hT,lT,";#font[32]{%c}  #font[22]{(mm)}   ",'Y');
    sT += string(hT); const char *hTXY = sT.c_str();
    hs.XY = new TH2D(hN,hTXY,256,-dX, dX,256,-dY,dY);
    if (idet==0) { // If CyMBaL
      snprintf(hT,lT,"%s;#font[32]{%c}r  #font[22]{(mm)}   ",dN,'X');
      sT = string(hT);
      snprintf(hT,lT,";#font[32]{%c}r  #font[22]{(mm)}   ",'Y');
      sT += string(hT); const char *hTXY = sT.c_str();
      double Xr0 = dX*cos(pi/8)*.9, dYr = dY*sin(pi/8);
      int cols[4] = {kGray,4,8,2};
      for (int zone = 0; zone<4; zone++) {
	snprintf(hN,lN,"hXYr%d",zone);
	hs.XYr[zone] = new TH2D(hN,hTXY,256,Xr0,dX,256,-dYr,dYr);
	hs.XYr[zone]->SetLineColor(cols[zone]);
      }
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
		   hs.XYr[0],hs.XYr[1],hs.XYr[2],hs.XYr[3]}; // CyMBaL specific
    int nH2s = sizeof(H2s)/sizeof(TH2D*);
    if (idet!=0) nH2s -= 4; // If !CyMBaL, cancel hs.XYr
    for (int ih = 0; ih<nH2s; ih++) {
      TAxis *ax = H2s[ih]->GetXaxis();
      ax->SetNdivisions(505);
      ax->SetLabelSize(.05); ax->SetLabelOffset(.006);
      ax->SetTitleSize(.06); ax->SetTitleOffset(.8);
      ax->SetMaxDigits(2);
      TAxis *ay = H2s[ih]->GetYaxis();
      ay->SetNdivisions(505);
      ay->SetLabelSize(.05); ay->SetLabelOffset(.006);
      ay->SetTitleSize(.06); ay->SetTitleOffset(.8);
      ay->SetMaxDigits(2);
    }
    if (isRec) {
      // ******************** RESIDUALS ********************
      int i01 = tag[1]=='0' ? 0 : 1;
      Resids &rs = resHs[i01][idet]; rs.dir = dSD;
      // ***** RANGES
      double dx = idet<2 ? 512 /* All MPGDs */ : 24 /* Si */; // in µm
      double dr = dx, dd = dx ,dz = dx, du = dx, dv = dx;
      double dphi = .8; // in mrad
      if (0x1<<idet&stripMode) { // ***** STRIPS: UPDATE RANGES DEPENDING ON iStrip
	if (idet==0) { // CyMBaL
	  // From: epic/compact/tracking/mpgd_barrel.xml
	  // <constant name="MMModuleLength"                         value="61.0*cm"/>
	  // 256-16-16 bins for the core part, 2x32 bins outside 
	  double modL = 610, deltaZ = modL/192, dz_phi = modL/2+32*deltaZ;
	  // For phi, X, etc... no way to have limit at bin edge, since we have
	  // 2 distinct sensitive surface radii, viz. 55.6755 and 57.8755
	  double dphi_Z = 2*pi/8/2; dphi_Z *= 1.20;
	  double RsurfMx = 578.755, dx_Z = dphi_Z*RsurfMx;
	  dz_phi *= 1000; dx_Z *= 1000; dphi_Z *= 1000; // mm -> µm, rad -> mrad
	  if (iStrip==1) dz = dz_phi;
	  else { dphi = dphi_Z; dx = dx_Z; }
	  dd = sqrt(4*dx_Z*dx_Z+dz_phi*dz_phi);
	}
	else if (idet==1) { // Outer
	  //<constant name="MPGDOuterBarrelModule_length"                  value="0.5*(MPGDOuterBarrelModule_zmin1 + MPGDOuterBarrelModule_zmin2 + MPGDOuterBarrelModule_zoverlap)"/>
	  //<constant name="MPGDOuterBarrelModule_zoverlap"                value="0*MPGDOuterBarrelFrame_width"/>
	  double zmin1 = 1640.5, zmin2 = 1740.5, modL = (zmin1+zmin2)/2;
	  //<constant name="MPGDOuterBarrelModule_width"                   value="36.0*cm"/>
	  double modW = 360;
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
      if (idet==3) { // Si: Looks like it has a very fine resolution in phi
	dphi = 0.08; 
      }
      double dRr = 1600;
      snprintf(hN,lN,"%c%s",'d',"X");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'X');
      rs.X =   new TH1D(hN,hT,256,-dx,dx);
      snprintf(hN,lN,"%c%s",'d',"Y");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'Y');
      rs.Y =   new TH1D(hN,hT,256,-dx,dx);
      snprintf(hN,lN,"%c%s",'d',"R");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'R');
      rs.R =   new TH1D(hN,hT,256,-dr,dr);
      if      (idet==0) { // If CyMBaL
	snprintf(hN,lN,"%c%s",'d',"Rr");
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'R');
	rs.Rr =   new TH1D(hN,hT,256,-dRr,dRr);
      }
      else if (idet==1) { // If Outer
	snprintf(hN,lN,"%c%s",'d',"Rr");
	snprintf(hT,lT,"%s;d#font[32]{%s}  #font[22]{(#mum)}   ",dN,
		 "Rc#delta#scale[1.2]{#varphi}");
	rs.Rr =   new TH1D(hN,hT,256,-dRr,dRr);
	snprintf(hN,lN,"%c%s",'d',"Ur");
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'U');
	rs.Ur =   new TH1D(hN,hT,256,-du,du);
	snprintf(hN,lN,"%c%s",'d',"Vr");
	snprintf(hT,lT,"%s;d#font[32]{%c}r  #font[22]{(#mum)}   ",dN,'V');
	rs.Vr =   new TH1D(hN,hT,256,-dv,dv);
      }
      snprintf(hN,lN,"%c%s",'d',"Z");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'Z');
      rs.Z =   new TH1D(hN,hT,256,-dz,dz);
      snprintf(hN,lN,"%c%s",'d',"dist");
      snprintf(hT,lT,"%s;d#font[32]{%c}  #font[22]{(#mum)}   ",dN,'D');
      rs.D =   new TH1D(hN,hT,512,-dd,dd);
      snprintf(hN,lN,"%c%s",'d',"phi");
      snprintf(hT,lT,"%s;d#scale[1.2]{#varphi}  #font[22]{(mrad)}   ", dN);
      rs.phi = new TH1D(hN,hT,512,-dphi,dphi); 
      snprintf(hN,lN,"%c%s",'d',"phir");
      snprintf(hT,lT,"%s;d#scale[1.2]{#varphi}r  #font[22]{(mrad)}   ",dN);
      rs.phir = new TH1D(hN,hT,512,-dphi,dphi); 
      TH1D *r1s[] =          {rs.X, rs.Y,rs.R,rs.Z,rs.D,rs.phi,
			      rs.Rr,         // MPGD specific
			      rs.phir,       // CyMBaL specific
			      rs.Ur,rs.Vr};  // Outer specific
      unsigned int flags[] = { 0xf,  0xf, 0xf, 0xf, 0xf,   0xf,
			       0x3,
			       0x1,
			       0x2,  0x2};
      int nr1s = sizeof(r1s)/sizeof(TH1D*);
      for (int ih = 0; ih<nr1s; ih++) {
	if (!(0x1<<idet&flags[ih])) continue;
	TAxis *ax = r1s[ih]->GetXaxis();
	ax->SetNdivisions(510);
	ax->SetLabelSize(.05); ax->SetLabelOffset(.006);
	ax->SetTitleSize(.06); ax->SetTitleOffset(.8);
	ax->SetMaxDigits(3);
	TAxis *ay = r1s[ih]->GetYaxis();
	ay->SetNdivisions(505);
	ay->SetLabelSize(.05); ay->SetLabelOffset(.006);
	ay->SetMaxDigits(3);
      }
    }
  }
  dSave->cd();
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
			     unsigned int detectorPattern, bool decompose)
{
  // ***** PARSE <iSimRec> ARG.
  if (iSimRec<0 || 2<iSimRec) {
    printf("** DrawphithZR: Unvalid <iSimRec> arg.(=%d): neither 0(=sim) nor 1(=rec) nor 2(=2nd coord of STRIP)\n",iSimRec); return;
  }
  if (!reconstruction && iSimRec) {
    printf("** DrawphithZR: Unvalid <iSimRec> arg.(=%d): recoEvents is simulationOnly\n",iSimRec); return;
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
    default: printf("** DrawphithZR: Unvalid <iSimRec> (=%d, not in [0,2])\n",
		    iSimRec); return;
    }
    Histos &hs = Hs[idet]; hs.dir->cd();
    // ***** TCANVAS
    string cS = string("c")+string(detectorNames[idet]);
    char cTag[] = "0";    // tag to cope w/ several distinct recoEvents object
    int iObj = iObjCreated; if (iObj) { *cTag += iObj%10; cS += string(cTag); }
    if (iSimRec)
      cS += iSimRec==2 ? string("Rec2") : string("Rec");
    const char *cN = cS.c_str();
    TCanvas *cEvents = new TCanvas(cN,cN);
    cEvents->Divide(2,2);
    gStyle->SetOptStat(10);
    // ***** LOOP ON (selected subset of) HISTOS
    TH2D *h2s[4] = { hs.phi,hs.th,hs.Z,hs.R };
    if (idet<2) h2s[3] = hs.Rr;
    for (int ih = 0; ih<4; ih++) {
      cEvents->cd(ih+1); TH2D *h2 = h2s[ih];
      TH1D *hproj = h2->ProjectionX(); hproj->Draw();
      hproj->SetMinimum(0); // Useful for hs.phi
      TAxis *ay = hproj->GetYaxis();
      ay->SetNdivisions(505);
      ay->SetLabelSize(.05); ay->SetLabelOffset(.006);
      ay->SetMaxDigits(2);
      SetPaveText(hproj);
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
    }
  }
  dSave->cd();
}
void recoEvents::DrawModules(int iSimRec, unsigned int detectorPattern, bool decompose)
{
  // ********** DRAW module# FOR EACH DETECTOR IN detectorPattern
  // ***** PARSE ARG.S
  if (iSimRec<0 || 2<iSimRec) {
    printf("** DrawModules: Unvalid <iSimRec> arg.(=%d): neither 0(=sim) nor 1(=rec) nor 2(=2nd coord of STRIP)\n",iSimRec); return;
  }
  if (!reconstruction && iSimRec) {
    printf("** DrawphithZR: Unvalid <iSimRec> arg.(=%d): recoEvents is simulationOnly\n",iSimRec); return;
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
    default: printf("** DrawModules: Unvalid <iSimRec> (=%d, not in [0,2])\n",
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
			       TCanvas *cPrv, int col) // Superimpose on pre-existing TCanvas cPrv
{
  // ********** DRAW (subset of) RESIDUALS FOR EACH DETECTOR IN detectorPattern
  // ***** PARSE ARG.S
  if (iRec<1 || 2<iRec) {
    printf("** DrawResiduals: Unvalid <iRec> arg.(=%d): neither 1(=1st coord.) nor 2(=2nd coord.)\n",iRec); return;
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
  gStyle->SetOptStat(1110);
  string prvFormat = string(gStyle->GetStatFormat());
  TDirectory *dSave = gDirectory;
  // ***** LOOP ON DETECTORS
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&pat)) continue;
    if (strip && !(0x1<<idet&stripMode)) { // 2nd coordinate: skip if not STRIP
      printf("** DrawResiduals: requesting <iRec> = 2 for idet = 0x%x (\"%s\") which is not among detectors w/ strips(=0x%x)\n",
	       0x1<<idet,detectorNames[idet],stripMode); continue;
    }
    TCanvas *cResids; if (!cPrv) {
      // ***** CANVAS
      string cS = string("c")+string(detectorNames[idet]);
      char cTag[] = "0";    // tag to cope w/ several distinct recoEvents object
      int iObj = iObjCreated; if (iObj) { *cTag += iObj%10; cS += string(cTag); }
      cS += strip ? string("Res2") : string("Res");
      const char *cN = cS.c_str();
      cResids = new TCanvas(cN,cN);
      cResids->Divide(2,2);
    }
    else
      cResids = cPrv;
    gStyle->SetStatFormat("6.3g");
    // ***** LOOP ON RESIDUALS (in subset)
    Resids &rs = resHs[strip][idet]; rs.dir->cd();
    TH1D *r1s[4] = {rs.X, rs.Z, rs.phi, rs.R};
    if      (idet==0) { // CyMBaL: precision is on 'phir' or 'Z' alternatively,
      // depending upon strip's coord. 'Rr': let's check it's a Dirac.
      r1s[2] = rs.phir;               r1s[3] = rs.Rr;
    }
    else if (idet==1) { // Outer:  precision is on 'Ur' or 'Vr' alternatively,
      // depending upon strip's coord. 'Rr' = "Rcos(dphi)": check it's a Dirac. 
      r1s[1] = rs.Ur; r1s[2] = rs.Vr; r1s[3] = rs.Rr;
    }
    for (int ih = 0; ih<4; ih++) {
      TH1D *h1 = r1s[ih];
      if (col>=0) h1->SetLineColor(col);
      cResids->cd(ih+1);
      if (cPrv) { h1->Draw("sames"); SetPaveText(h1,1); }
      else      { h1->Draw();        SetPaveText(h1,0); }
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
void SetPaveText(TH1 *h, int mode)
{
  TPad *pad = (TPad*)gPad; pad->Update();

  int col = h->GetLineColor();

  int optStat = gStyle->GetOptStat();
  double dY, X1; if (optStat==10)   { dY = .10; X1 = .70; }
  else           if (optStat==1110) { dY = .22; X1 = .55; }

  TPaveStats *st; if ((st = (TPaveStats*)h->GetFunction("stats"))) {
    st->SetX2NDC(.995); st->SetX1NDC(.70);
    st->SetY2NDC(.995-mode*dY); st->SetY1NDC(.995-(mode+1)*dY);
    //st->SetOptStat(1000010);
    st->SetTextColor(col);
    st->Draw();
  }

  TPaveText *tit; if ((tit = (TPaveText*)pad->GetPrimitive("title"))) {
    if (mode==0) { // First histo in pad
      tit->SetX1NDC(.25); tit->SetX2NDC(.40);
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
      printf(" * SetMinima: %s->Setinimum(%f);\n",h->GetName(),min);
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
  ana2.requireQuality = 1; // REQUIRE quality IN SimTrackerHit's
  ana2.Loop();
  // ***** CyMBAL
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
