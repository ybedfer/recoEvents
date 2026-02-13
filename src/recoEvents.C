#define recoEvents_cxx
#include "recoEvents.h"
#include <TH2.h>
#include <TVector3.h>
#include <TRotation.h>

#include <stdio.h>

// Globals
unsigned int recoEvents::nObjCreated = 0;

// Interfaces
unsigned int cExtension(double const* lpos, double const* lmom, // Input subHit
			double rT,                              // Target radius
			int direction, double dZ, double startPhi,
			double endPhi, // Module parameters
			double* lext);

void recoEvents::Loop(int nEvents, int firstEvent)
{
  //   In a ROOT session, you can do:
  //      root> .L recoEvents.so
  //      root> recoEvents ana(events)
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> ana.Loop();       // Loop on all entries
  //

  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
#define FromChain
  /*
    #ifdef FromChain
    fChain->SetBranchStatus("*",0);  // disable all branches
    const char *branchNames[N_DETs] = {
    "MPGDBarrelHits","OuterMPGDBarrelHits","VertexBarrelHits","SiBarrelHits"};
    for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&processedDetectors)) continue;
    fChain->SetBranchStatus(branchNames[idet],1);  // activate branchname
    }
    #endif
  */
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  Long64_t nentries = nEvents==0 ? fChain->GetEntriesFast() : nEvents;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=firstEvent; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
#ifdef FromChain
    nb = fChain->GetEntry(jentry);   nbytes += nb;
#else
    for (int idet = 0; idet<N_DETs; idet++) {
      if (!(0x1<<idet&processedDetectors)) continue;
      branches[idet]->GetEntry(jentry);   nbytes += nb;
    }
#endif

    evtNum = eventHeader->at(0).eventNumber;

    //int treenumber = fChain->GetTreeNumber();
    // if (Cut(ientry) < 0) continue;

    // ***** EVENT SELECTION
    int nMCs = mcParticles->size(); hMult->Fill(nMCs);
    if      (requireNMCs>0)  { if (nMCs!=requireNMCs) continue; }
    else if (requireNMCs<0)  { if (nMCs< requireNMCs) continue; }

    for (int idet = 0; idet<N_DETs; idet++) {
      if (!(0x1<<idet&processedDetectors)) continue;
      initDetEvent();     // Init requirement for current detector
      // *************** LOOP ON SELECTED DETECTORS
      int nHits = hits[idet]->size(); // SimHits
      vector<SimTrackerHitData> coalescedHs; map<int,int> sim2coa;
      for (int ih = 0; ih<nHits; ih++) {
	// ********** LOOP ON sim HITS: COALESCING AND EXTENDING
	// -> "coalescedHs"
	SimTrackerHitData &hit = hits[idet]->at(ih);
	int jh = ih; vector<int> coalesced; int nCoaHs = coalescedHs.size();
	if (0x1<<idet&stripMode) {
	  while (++jh<nHits && samePMO(idet,ih,jh)) {
	    if (extrapolate(idet,ih,jh)) {
	      if (coalesced.empty()) {
		coalesced.push_back(ih); sim2coa[ih] = nCoaHs;
	      }
	      coalesced.push_back(jh); sim2coa[jh] = nCoaHs;
	    }
	    else break;
	  }
	}
	if (coalesced.size()<2) {
	  if (hit.quality==0 &&
	      (0x1<<idet&stripMode)) { // Extend to edge of pathLength
	    SimTrackerHitData hext; extend(idet,ih,hext);
	    coalescedHs.push_back(hext);
	  }
	  else // Do not extend secondary: its position is, a priori, well
	    // represented by its SimTrackerHit position.
	    coalescedHs.push_back(hit);
	  sim2coa[ih] = nCoaHs;
	}
	else {
	  SimTrackerHitData hext; if (!coalesce(idet,coalesced,hext)) return;
	  sim2coa[ih] = nCoaHs; coalescedHs.push_back(hext);
	  ih += coalesced.size()-1;
	}
	updateDetEvent(idet,ih);
      }
      finaliseDetEvent();

      int nCoaHs = coalescedHs.size();
      for (int ih = 0; ih<nCoaHs; ih++) {
	// ********** LOOP ON coalescedHs
	SimTrackerHitData &hit = coalescedHs[ih];
	unsigned int status = getStatus(idet,ih,sim2coa);
	debugHit(idet,ih,nCoaHs,hit,status);
	// ***** MODULE SELECTION
	if (!moduleSelection(hit.cellID)) continue;
	// ***** HIT SELECTION
	if (status!=0x3) continue;
	// ***** FILL sim HISTO
	// - It's rather "coa" rather than "sim"
	const Vector3d &pos = hit.position;
	double X = pos.x, Y = pos.y, Z = pos.z;
	fillHit(0,idet,X,Y,Z,hit.eDep,hit.cellID);
	debugHit(idet,hit);
      }
      if (!reconstruction) continue;

      // ********** ASSOCIATION: BUILD rec -> coa MAP
      int nArhs = arhs[idet]->size(); // Associations raw hit
      int nAshs = ashs[idet]->size(); // Associations sim hit
      int nRecs = recs[idet]->size(); // Rec hits
      if (!(nArhs||nAshs||nHits||nRecs)) continue;
      int nARhs = aRhs[idet]->size(); // Associations Rec Hit -> raw hit
      map<int,vector<int>,less<int>> rec2coas;
      debugAssoc(idet);
      if (nAshs!=nArhs) {
	printf("Evt #%5d Warning: det %d: ash(%d)!=arh(%d)\n",
	       evtNum,idet,nAshs,nArhs);
      }
      else {
	// raw -> Rec
	// raw can only be associated to one Rec, by construction...
	// ...and vice-versa, by convention imprinted in edm4hep::TrackerHit
	map<int,int,less<int>> raw2rec;
	for (int iR = 0; iR<nARhs; iR++) {
	  podio::ObjectID &aRh = aRhs[idet]->at(iR); raw2rec[aRh.index] = iR;
	}
	// Rec <-> coas
	for (int ih = 0; ih<nArhs; ih++) {
	  // ***** LOOP ON raw AND sims
	  // raw: "rIndex"
	  podio::ObjectID &arh = arhs[idet]->at(ih); //RawHitAssociations_rawHit
	  int rIndex = arh.index;
	  // raw -> rec: "RIndex"
	  map<int,int>::const_iterator ir = raw2rec.find(rIndex);
	  if (ir==raw2rec.end())
	    // One raw may be lost here, but that lost one is then associated to
	    // the same sims as the retained one.
	    continue;
	  int RIndex = ir->second;
	  // sim: "sIndex"
	  podio::ObjectID &ash = ashs[idet]->at(ih); //RawHitAssociations_simHit
	  int sIndex = ash.index;
	  // sim -> coa: "cIndex"
	  map<int,int>::const_iterator is = sim2coa.find(sIndex);
	  if (is==sim2coa.end()) continue;
	  int cIndex = is->second;
	  // Rec -> coa
	  map<int,vector<int>>::iterator im = rec2coas.find(RIndex);
	  if (im==rec2coas.end()) {
	    vector<int> coas; coas.push_back(cIndex); rec2coas[RIndex] = coas;
	  } else {
	    // Not twice same SimHit associated to RecHit.
	    vector<int> &coas = im->second;
	    int is, match; for (int is=match = 0; is<(int)coas.size(); is++) {
	      if (coas[is]==cIndex) { match = 1; break; }
	    }
	    if (!match) coas.push_back(cIndex);
	  }
	}
	debugAssoc(idet,raw2rec,sim2coa,rec2coas);
      }

      for (int ir = 0; ir<nRecs; ir++) {
	// ********** LOOP ON Rec HITS
	edm4eic::TrackerHitData &rec = recs[idet]->at(ir);
	const Vector3f &pos = rec.position;
	double X = pos.x, Y = pos.y, Z = pos.z;
	debugRec(idet,ir);
	// ********** RESIDUALS
	// ***** REFERENCE TO coa HIT
	map<int,vector<int>,less<int>>::const_iterator im = rec2coas.find(ir);
	if (im==rec2coas.end()) {
	  printf("#%5d Warning: det %d: rec %d not associated\n",
		 evtNum,idet,ir);
	}
	else {
	  const vector<int> &coas = im->second;
	  int is, selecRec; for (is=selecRec = 0; is<(int)coas.size(); is++) {
	    int cIndex = coas[is];
	    if (cIndex<0 || nHits<=cIndex) {
	      printf("#%5d Warning: det %d: Rec %d -> sim %d\n",
		     evtNum,idet,ir,cIndex);
	    }
	    else {
	      SimTrackerHitData &hit = coalescedHs[cIndex];
	      debugHitRec(idet,is,(int)coas.size(),cIndex,hit,ir);
	      // ***** MODULE SELECTION
	      if (!moduleSelection(rec.cellID)) continue;
	      // ***** HIT SELECTION
	      unsigned int status = getStatus(idet,cIndex,sim2coa);
	      if (status!=0x3) continue;
	      // ***** FILL RESIDUAL
	      selecRec = 1;
	      const Vector3d &psim = hit.position;
	      fillResids(idet,pos,psim,rec.cellID);
	    }
	  }
	  if (selecRec) {  // ***** FILL rec HISTOS
	    fillHit(1,idet,X,Y,Z,rec.edep,rec.cellID);
	    debugRec(idet,ir);
	  }
	}
      }
    }
  }
}
unsigned int recoEvents::getStatus(int idet, int ih)
{
  // Returned status is OR of conformity to quality and PDG requirements.
  // 0x1: PDG OK.
  // 0x2: quality OK.
  unsigned status = 0;
  SimTrackerHitData &hit = hits[idet]->at(ih);
  if (!requireQuality || hit.quality==0)   status |= 0x2;
  podio::ObjectID &amc = amcs[idet]->at(ih); int mcIdx = amc.index;
  MCParticleData &part = mcParticles->at(mcIdx);
  if (!requirePDG || part.PDG==requirePDG) status |= 0x1;
  return status;
}
unsigned int recoEvents::getStatus(int idet, int ih, map<int,int> &sim2coa)
{
  unsigned int status = 0;
  for (auto is = sim2coa.cbegin(); is!=sim2coa.cend(); is++) {
    if (is->second==ih) {
      int jh = is->first; status = getStatus(idet,jh);
      break;
    }
  }
  return status;
}
void recoEvents::fillHit(int simOrRec, int idet,
			 double X, double Y, double Z, double eDep,
			 unsigned long cellID)
{
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  // ***** STRIP: Convert stripID -> strip. Is it valid? 
  if (!parseStrip(idet,simOrRec,strip)) return;
  double R2 = X*X+Y*Y, R = sqrt(R2);
  double phi = atan2(Y,X);
  double rho2 = R2+Z*Z, rho = sqrt(rho2);
  const double pi = TMath::Pi();
  double theta = rho?acos(Z/rho):999*pi;
  Histos *hs; if (simOrRec==0) hs = &(simHs[idet]);
  else                         hs = &(recHs[strip][idet]);
  //  printf("======= %d %d %u %p\n",simOrRec,idet,strip,hs);
  hs->X->Fill(X,div); hs->Y->Fill(Y,div); hs->R->Fill(R,div);
  hs->RA->Fill(R,div);
  hs->Z->Fill(Z,div);
  hs->phi->Fill(phi,div); hs->th->Fill(theta,div);
  hs->mod->Fill(module,div);
  double gain = simOrRec ? gains[idet] : 1;
  hs->eDep->Fill(eDep*1e6/gain,div);
  // 2D histos
  hs->thphi->Fill(phi,theta);
  hs->XY->Fill(X,Y); hs->ZR->Fill(Z,R);

  if      (idet==0) { // Special CyMBaL: fill reduced Radius
    double Xr, Yr, Zr, Rr, phir;
    g2lCyMBaL(X,Y,Z,div,Xr,Yr,Zr,Rr,phir);
    hs->phir->Fill(phir,div);
    hs->Rr->Fill(Rr,div);
    hs->xyr->Fill(Rr*phir,Zr);
  }
  else if (idet==1) { // Special Outer: fill Rcosphi, Ur, Vr
    double Rcdphi, Xr, Yr, Zr, Ur, Vr;
    g2lOuter(X,Y,Z,module,Rcdphi,Xr,Yr,Zr,Ur,Vr);
    hs->Rr->Fill(Rcdphi,div);
    hs->Ur->Fill(Ur,div); hs->Vr->Fill(Vr,div);
    hs->xyr->Fill(Xr,Yr);
  }
}
void recoEvents::fillResids(int idet, const Vector3f &pos, const Vector3d &psim, unsigned long cellID)
{
  int doPrint = 0;
  if (verbose&0x1000<<idet || evtNum==evtToDebug) {
    doPrint = 1; if (verbose&0x10000) doPrint = 2;
  }
  double X =  pos.x,  Y =  pos.y,  Z =  pos.z;
  double Xs = psim.x, Ys = psim.y, Zs = psim.z;
  double dX = 1000*(X-Xs), dY = 1000*(Y-Ys), dZ = 1000*(Z-Zs); // Residuals
  double R2 = X*X+Y*Y,      R = sqrt(R2);
  double Rs2 = Xs*Xs+Ys*Ys, Rs = sqrt(Rs2);
  double dR = 1000*(R-Rs);
  double phi = atan2(Y,X), phis = atan2(Ys,Xs);
  double dphi = 1000*(phi-phis);
  double D = sqrt(R2+Z*Z), Ds = sqrt(Rs2+Zs*Zs), dD = 1000*(D-Ds);
  if (doPrint>1)
    printf(" dX,dY,dZ: %.2f,%.2f,%.2f",dX,dY,dZ);
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  // ***** STRIP: Convert stripID -> strip. Is it valid? 
  if (!parseStrip(idet,1,strip)) return;
  Resids &rs = resHs[strip][idet];
  rs.X->Fill(dX); rs.Y->Fill(dY); rs.Z->Fill(dZ); rs.R->Fill(dR); rs.D->Fill(dD);
  rs.phi->Fill(dphi);
  if      (idet==0) { // CyMBaL specific
    double Xd, Yd, Zd, Rr, phir; g2lCyMBaL(X, Y, Z, div,Xd,Yd,Zd,Rr, phir);
    double Rrs, phirs;           g2lCyMBaL(Xs,Ys,Zs,div,Xd,Yd,Zd,Rrs,phirs);
    if (doPrint>1)
      printf("Rr,Rrs: %.2f,%.2f\n",Rr,Rrs);
    double dRr = 1000*(Rr-Rrs), dphir = 1000*(phir-phirs);
    rs.Rr->Fill(dRr);
    int section = module/8; int staveType = (section==1 || section==2)?0:1;
    double radius = radii[0][staveType];
    rs.Rphir->Fill(radius*dphir);
    if (doPrint) {
      if (strip==1 && fabs(dZ)>480 || strip==0 && fabs(dphir)>1) {
	//if (strip==1 && fabs(dZ)>0) {
	printf("#%5d 0x%08lx,0x%08lx  %7.2f,%7.2f,%8.2f %7.3f mm %7.3fπ\n",
	       evtNum,cellID&0xffffffff,cellID>>32,X,Y,Z,Rr,phir/TMath::Pi());
	char text[] ="dφ 1.810 mrd";
	//           "dZ 1000 µm";
	size_t lt = strlen(text)+1;
	if (strip==1) snprintf(text,lt,"dZ %4.0f µm",dZ);
	else          snprintf(text,lt,"dφ %5.3f mrd",dphir);
	printf("%16s %12s %7.2f,%7.2f,%8.2f R,φ %7.3f mm %7.3fπ  %s\n",
	       "SimHit","",Xs,Ys,Zs,Rrs,phirs/TMath::Pi(),text);
      }
    }
  }
  else if (idet==1) { // Outer specific
    double Rcdphi,  Xr,  Yr,  Zr,  Ur,  Vr;  g2lOuter(X, Y, Z, module,Rcdphi, Xr, Yr, Zr, Ur, Vr);
    double Rcdphis, Xrs, Yrs, Zrs, Urs, Vrs; g2lOuter(Xs,Ys,Zs,module,Rcdphis,Xrs,Yrs,Zrs,Urs,Vrs);
    //if (doPrint)
    //printf("Rr,Rrs: %.2f,%.2f, Ur,Urs: %.2f,%.2f, Vr,Vrs: %.2f,%.2f\n",
    //	     Rcdphi,Rcdphis,Ur,Urs,Vr,Vrs);
    double dRcdphi = 1000*(Rcdphi-Rcdphis);
    rs.Rr->Fill(dRcdphi);
    double dUr = 1000*(Ur-Urs), dVr = 1000*(Vr-Vrs);
    rs.Ur->Fill(dUr); rs.Vr->Fill(dVr);
    //rs.xyr->Fill(Xrs,Yrs,strip?dUr:dVr);
    if (doPrint) {
      if (strip==0 && fabs(dUr)>580 || strip==1 && fabs(dVr)>580) {
      //if (strip==1 && fabs(dU)>0) {
	printf("#%5d 0x%08lx,0x%08lx  %7.2f,%7.2f,%8.2f U,V %8.3f,%8.3f mm\n",
	       evtNum,cellID&0xffffffff,cellID>>32,X,Y,Z,Ur,Vr);
	char text[] ="dU 1000 µm"; size_t lt = strlen(text)+1;
	if (strip==1) snprintf(text,lt,"dV %4.0f µm",dVr);
	else          snprintf(text,lt,"dU %4.0f µm",dUr);
	printf("%16s %12s %7.2f,%7.2f,%8.2f U,V %8.3f %8.3f mm %s\n",
	       "SimHit","",Xrs,Yrs,Zs,Urs,Vrs,text);
      }
    }
  }
  else if (doPrint) printf("\n");
}
void recoEvents::g2lCyMBaL(double X,   double Y,   double Z, unsigned int div,
			   double &Xr, double &Yr, double &Zr,
			   double &Rr, double &phir)
{
  // Transform to centre of curvature of cylindrical tiles
  // and rotate to phi = 0
  int section = div>>3;
  int iphi = div%8, jphi = iphi%2; double phic = iphi*TMath::Pi()/4;
  //<constant name="MMRadial_offset"                        value="1.0*cm"/>
  double offset = 10; // mm
  double rc = (2*jphi-1)*offset/2;    // ...flip sign of offset
  double x1, y1; // Coordinates of the centre of curvature
  x1 = rc * std::cos(phic); y1 = rc * std::sin(phic);
  TVector3 v1(x1,y1,0);
  Xr = X-x1; Yr = Y-y1; Rr = sqrt(Xr*Xr+Yr*Yr);
  TVector3 V(X,Y,0); V -= v1;
  TRotation r;
  r.SetXAxis(TVector3(cos(-phic),sin(-phic),0));
  r.SetYAxis(TVector3(-sin(-phic),cos(-phic),0)); V *= r;
  Xr = V(0); Yr = V(1); phir = atan2(Yr,Xr);
  // Zr
  Zr = Z+sectionDZs[0][section]; if (section>1) Zr *= -1;
}
void recoEvents::g2lCyMBaL(double Px, double Py, double Pz, unsigned int div,
			   double *lmom)
{
  int section = div>>3;
  int iphi = div%8, jphi = iphi%2; double phic = iphi*TMath::Pi()/4;
  double crot = std::cos(phic), srot = std::sin(phic);
  lmom[0] = crot*Px+srot*Py; lmom[1] = -srot*Px+crot*Py;
  lmom[2] = section<2 ? +Pz : -Pz;
}
void recoEvents::l2gCyMBaL(double *lpos, int div, double *gpos)
{
  int section = div>>3;
  int iphi = div%8, jphi = iphi%2; double phic = iphi*TMath::Pi()/4;
  //<constant name="MMRadial_offset"                        value="1.0*cm"/>
  double offset = 10; // mm
  double rc = (2*jphi-1)*offset/2;    // ...flip sign of offset
  double x1, y1; // Coordinates of the centre of curvature
  x1 = rc * std::cos(phic); y1 = rc * std::sin(phic);
  double Xr = lpos[0], Yr = lpos[1], Zr = lpos[2];
  double &X = gpos[0], &Y = gpos[1], &Z = gpos[2];  
  TVector3 V(Xr,Yr,0);
  TRotation r;
  r.SetXAxis(TVector3( cos(phic),sin(phic),0));
  r.SetYAxis(TVector3(-sin(phic),cos(phic),0)); V *= r;
  TVector3 v1(x1,y1,0);
  V += v1;
  X = V(0); Y = V(1);
  // Z
  if (section>1) Zr *= -1;
  Z = Zr-sectionDZs[0][section];
}
void recoEvents::g2lOuter(double X, double Y, double Z, unsigned int module,
			  double &Rcdphi, double &Xr, double &Yr, double &Zr,
			  double &Ur, double &Vr,
			  double *rot)
{
  // Rotate to phi = 0
  int iphi = module>>1; double phic = iphi*TMath::Pi()/6;
  if (rot) *rot = phic;
  //<constant name="MPGDOuterBarrelModule_zmin1"     value="164.5*cm"/>
  //<constant name="MPGDOuterBarrelModule_zmin2"     value="174.5*cm"/>
  TVector3 V(X,Y,Z);
  TRotation r;
  r.SetXAxis(TVector3(cos(-phic),sin(-phic),0));
  r.SetYAxis(TVector3(-sin(-phic),cos(-phic),0));
  V *= r;
  Rcdphi = V(0); Xr = V(1); Yr = V(2);
  int section = module%2;
  if (section==0) { Xr *= -1; Yr *= -1; }
  
  Vr = (Yr+Xr)/sqrt(2); Ur = (Yr-Xr)/sqrt(2);
  // Local Z is perpendicular to (Xr,Yr)=(Ur,Vr) plane
  // => Rcdphi - distance to beam axis.
  Zr = Rcdphi-radii[1][0];
  // Local Y is global Z, w/ a shift depending upon section = module%3
  Yr -= sectionDZs[1][section];
}
void recoEvents::g2lOuter(double Px, double Py, double Pz, unsigned int module,
			  double *lmom)
{
  // Rotate to phi = 0
  int iphi = module>>1; double phic = iphi*TMath::Pi()/6;
  double crot = std::cos(phic), srot = std::sin(phic);
  double Pxr = crot*Px+srot*Py, Pyr = -srot*Px+crot*Py;
  if (module%2==0) { Pz *= -1; Pyr *= -1; }
  lmom[0] = Pyr; lmom[1] = Pz; lmom[2] = Pxr;
}
void recoEvents::parseCellID(int idet, unsigned long ID,
			    unsigned int &module, unsigned int &div,
			    unsigned int &strip)
{
  // ***** MODULE
  // mpgd_barrel.xml:      <id>system:8,layer:4,module:12, ...
  // mpgd_outerbarrel.xml: <id>system:8,layer:4,module:12, ...
  // silicon_barrel.xml:   <id>system:8,layer:4,module:12, ...
  // vertex_barrel.xml:    <id>system:8,layer:4,module:12, ...
  // => All have same module specif.
  // system: CyMBaL = 61, Outer = 64, ... (It's not checked that this is what we get in "ID")
  module = (ID>>12)&0xfff;
  // ***** divISION
  //       iRec = (strip coordinate: 0 = measurement coord, 1 = orthogonal
  if      (idet==0) { // CyMBaL
    //      <id>system:8,layer:4,module:12,sensor:30:2,phi:-16,z:-16</id>
    if (module>15) {
      // Forward sectors are rotated by phi
      int sector = module/8, iphi = module%8; module = 8*sector+(4+iphi)%8;
    }
    div = module;
    strip = nSensitiveSurfaces==5 ? ID>>28&0xf : ID>>30&0x3;
  }
  else if (idet==1) { // µRWELL
    div = module%2;      // "div" = parity
    strip = nSensitiveSurfaces==5 ? ID>>28&0xf : ID>>30&0x3;
  }
  else if (idet==2) {   // Vertex
    int layer = (ID>>8)&0xf; // "div" = log_2(layer)
    int bit; for (bit = 0, div = 3 /* unphysical default */; bit<=2; bit++) {
      if ((0x1<<bit&layer)==layer) { div = bit; break; }
    }
    strip = 0; // Not relevant
  }
  else {
    div = ID%2;
    strip = 0; // Not relevant
  }
}
bool recoEvents::parseStrip(int idet, int simOrRec, unsigned int &strip)
{
  // - Convert input stripID (= 0x1 or 0x2) to "strip" (= 0 or 1)
  // - Check it's valid: Depends upon...
  // ..."stripMode" of current "idet"
  // ...whether it's SimHit or RecHit
  // ..."nSensitiveSurfaces"
  bool hasStrips = stripMode&0x1<<idet;
  if (simOrRec && hasStrips) {
    strip -= 1;
    if (strip<0 || 1<strip) {
      printf("** parseStrip: (stripMode=0x%x) idet=%d stripID=0x%x => strip = %d not in [0,1]\n",
	     stripMode,idet,strip+1,strip);
      strip = 0;
      return false;
    }
  }
  else if (nSensitiveSurfaces==1 || !hasStrips) {
    if (strip!=0) {
      printf("** fillHit: (stripMode=0x%x) idet=%d stripID=0x%x not =0\n",
	     stripMode,idet,strip);
      strip = 0;
      return false;
    }
  }
  return true;
}
bool recoEvents::samePMO(int idet, int ih, int jh) {
  // Check <ih> and <jh> share same P[ARTICLE], M[ODULE], O[RIGIN]...
  // ...TO WHICH WE ADD that SUBVOLUME SHOULD DIFFER. This, to exclude
  // secondaries that otherwise would pass the "extrapolate" step, because on
  // the one hand, they originate from the same primary (and are hence close to
  // one another) and on the other hand, "extrapolate" method is simplistic (
  // at variance to EICrecon's "MPGDTrackerDigi", it does not take into the
  // pathLengths of the hit to be coalesced).
  SimTrackerHitData &hit = hits[idet]->at(ih), &hjt = hits[idet]->at(jh);
  int mcIdx = amcs[idet]->at(ih).index;
  unsigned long cID = hit.cellID, vID = cID&0xffffffff;
  int module = cID>>12&0xfff;
  int quality = hit.quality;
  int mcJdx = amcs[idet]->at(jh).index;
  unsigned long cJD = hjt.cellID, vJD = cJD&0xffffffff;
  int modvle = cJD>>12&0xfff;
  int qualjty = hjt.quality;
  if ((verbose&0x100<<idet) || evtNum==evtToDebug) doDebug |= 0x100;
  if  (doDebug&0x100) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%5d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
    int nHits = hits[idet]->size();
    printf("%2d/%d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx), %2d/%d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx)\n",
	   ih,nHits,mcIdx,module,quality,vID,hit.cellID>>32,
	   jh,nHits,mcJdx,modvle,qualjty,vJD,hjt.cellID>>32);
  }
  return mcJdx==mcIdx && modvle==module && qualjty==quality && vJD!=vID;
}
  //#define DEBUG_EXTRAP
#ifdef DEBUG_EXTRAP
#include <TGraph.h>
#endif
bool recoEvents::extrapolate(int idet, int ih, int jh)
{
  // - Extrapolate hits along their momentum by half-distance.
  // - Return true if compatibility w/in 25 µm (arbitrary value).
  // - Next step is expected to be a coalescing of the compatible hits.
  //  Coalescing which will assign to the resulting (coalesced) hit as a
  //  position, the half-sum of the hits (this, most of the time, see method
  //  "coalesce" for the details), typically sitting on the mid-plane of the
  //  sensitive volume (= union of SUBVOLUMES).
  // => Question: is the above compatibility criterion well grounded? What is
  //  missing is the notion of traversing particle in MPGDTrackerDigi. Which
  //  I illustrate w/ the example of a delta-ray emitted w/in one SUBVOLUME. In
  //  such a case, as far as I understand, DD4hep produces two hits:
  //   + Primary, sitting not on mid-plane,
  //   + Secondary.
  //  The primary has a continuation into the next (or previous) SUBVOLUME,
  //  compatible w/ it in the above sense. Then one has to distinguish two
  //  sub-cases, depending on the ordering along the thickness axis:
  //   i) Secondary < primary   < continuation,
  //  ii) Primary   < secondary < continuation.
  //   Sub-case (i) is the only legitimate primary+continuation coalescing.
  /// - Here, we can identify cases where there is a secondary based on the
  //   sole info contained in the input hit (w/o having to scan the list of all
  //   SimHits): would suffice to check the ''depth'' of the hit along the
  //   thickness axis against expectation = thickness, which will also catch
  //   cases where the particle exits the sensitive volume through the edge (
  //   see "requireTraversing") and compare the hit position to the mid-plane
  //   to determine which sub-case we are in.
  SimTrackerHitData &hit = hits[idet]->at(ih), &hjt = hits[idet]->at(jh);
  const Vector3d &pos = hit.position, &pqs = hjt.position;
  double Mx = pos.x, My = pos.y, Mz = pos.z;
  double Nx = pqs.x, Ny = pqs.y, Nz = pqs.z;
  double dist = sqrt((Mx-Nx)*(Mx-Nx)+(My-Ny)*(My-Ny)+(Mz-Nz)*(Mz-Nz));
  const Vector3f &mom = hit.momentum, &mqm = hjt.momentum;
  double Px = mom.x, Py = mom.y, Pz = mom.z, P = sqrt(Px*Px+Py*Py+Pz*Pz);
  double Qx = mqm.x, Qy = mqm.y, Qz = mqm.z, Q = sqrt(Qx*Qx+Qy*Qy+Qz*Qz);
  double u = dist/2/P, Ex = Mx+u*Px, Ey = My+u*Py, Ez = Mz+u*Pz;
  double v = dist/2/Q, Fx = Nx-v*Qx, Fy = Ny-v*Qy, Fz = Nz-v*Qz;
  double dext = sqrt((Ex-Fx)*(Ex-Fx)+(Ey-Fy)*(Ey-Fy)+(Ez-Fz)*(Ez-Fz));
  // Cut is made depending upon P...
  double pMeV = P*1000, dMx = 2*(.05+1e-5*exp(-(pMeV-25)/2.6));
  // and relaxed when one of the participant is a thin SUBVOLUME
  int subVolID = hit.cellID>>28&0xf, subVolJD = hjt.cellID>>28&0xf;
  if (subVolID<3 || subVolJD<3) dMx *= 1.5;
  bool ok = dext<dMx; // 25 µm
#ifdef DEBUG_EXTRAP
  // Try and determine appropriate parameters for the P dependence of the Cut
  if (dext>.005 && dext<1.00) {
    double pathDepth, depth, pathDfpth, dfpth;
    bool traversing = checkTraversing(idet,hit,pathDepth,depth);
    bool traversjng = checkTraversing(idet,hjt,pathDfpth,dfpth); 
    if (traversing && traversjng) {
      double p = P*1000;
      int ip = P*1000 < 1 ? int(p*5) : 5+int(p/2); ip -= 1;
      static double sps[10], sds[10]; static int first = 1, s1s[10], nOKs[10];
      if (ip<0) ip = 0; if (ip>=10) ip = 9;
      if (first) {
	first = 0; for (int i = 0; i<10; i++) { sps[i]=sds[i] = 0; s1s[i]=nOKs[i] = 0; }
      }
      sps[ip] += p; sds[ip] += dext; s1s[ip]++; if (ok) nOKs[ip]++;
      printf("==== #%5d,%d,%02d,%02d dext = %.3f/%.3f %s P = %6.3f MeV %d %6.3f MeV %.3f %3d  %.2f%%\n",
	     evtNum,idet,ih,jh,dext,dMx,dext<dMx?"OK":"  ",P*1000,ip,sps[ip]/s1s[ip],sds[ip]/s1s[ip],s1s[ip],100.*nOKs[ip]/s1s[ip]);
      if (evtNum==15991) {
	TGraph *g = new TGraph(); g->SetName("gdP");
	for (int ip = 0; ip<10; ip++) {
	  printf(" %d P = %6.3f MeV d %.3f %3d %3d %.2f%%\n",
		 ip,sps[ip]/s1s[ip],sds[ip]/s1s[ip],s1s[ip],nOKs[ip],100.*nOKs[ip]/s1s[ip]);
	  g->SetPoint(ip,sps[ip]/s1s[ip],sds[ip]/s1s[ip]);
	}
	g->Draw("AL");
      }
    }
  }
#endif
  // Debugging
  if ((verbose&0x100<<idet) || evtNum==evtToDebug) doDebug |= 0x100;
  if  (doDebug&0x100) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%5d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
  }
  bool doPrint = doDebug&0x100;
  if (!doPrint && !ok) {
    // Do print debugging message, even if not requested, when extrpolation
    // turns out to fail in case status is 0x3 and despite quality requirement.
    unsigned status = getStatus(idet,ih);
    doPrint = (status==0x3) && requireQuality==2;
  }
  if (doPrint) {
    printf("%d,%d dext %f\n",ih,jh,dext);
    printf("%2d: 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f\n",
	   ih,hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,Ex,Ey,Ez);
    printf("%2d: 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f\n",
	   jh,hjt.cellID&0xffffffff,hjt.cellID>>32,Nx,Ny,Nz,Fx,Fy,Fz);
  }
  if (ok && requireTraversing) {
    if (idet==0) {
      double pathDepth, depth, pathDfpth, dfpth;
      bool traversing = checkTraversing(idet,hit,pathDepth,depth);
      bool traversjng = checkTraversing(idet,hjt,pathDfpth,dfpth); 
      unsigned int continuation = 0;
      if (traversing) continuation |= 0x1;
      if (traversjng) continuation |= 0x2;
      for (int ij = 0; ij<2; ij++) {
	if (continuation&0x1<<ij) continue;
	// Try and rescue cases where traversing to one side.
	double Rr = ij==0 ? depth : dfpth;
	bool is0x3 = subVolID==0x3 || subVolJD==0x3;
	int module = hit.cellID>>12&0xfff;
	int section = module/8; int staveType = (section==1 || section==2)?0:1;
	double midPlane = radii[0][staveType];
	bool isInner = is0x3 ? Rr > midPlane : Rr < midPlane;
	if (isInner) continuation |= 0x4<<ij;
      }
      if (doPrint)
	printf("(Evt%d)0x%04x %d: 0x%x %.3f %.3f %.3f, %d: 0x%x %.3f %.3f %.3f\n",
	       evtNum,continuation,
	       ih,subVolID,depth,hit.pathLength,pathDepth,
	       jh,subVolJD,dfpth,hjt.pathLength,pathDfpth);
      ok = (continuation&0x5) && (continuation&0xa);
    }
    else if (idet==1) {
      double pathDepth, depth, pathDfpth, dfpth;
      bool traversing = checkTraversing(idet,hit,pathDepth,depth);
      bool traversjng = checkTraversing(idet,hjt,pathDfpth,dfpth); 
    }
  }
  return ok;
}
bool recoEvents::coalesce(int idet, vector<int> coalesced, SimTrackerHitData &hext) {
  // Basic idea: simplify the problem by disregarding as much as possible the
  // thin SUBVOLUMES...
  // ...Yet still taking into account their eDep.
  // I)  If SUBVOLUMES 0x3 *and* 0x4, position of calesced hit is half way
  //  between 0x3 and 0x4 extrema.
  // II) If SUBVOLUMES 0x3 *or*  0x4, we fall back on the single SUBVOLUME case
  //  and "extend" the single subHit.
  // III) Else, rare case where only thin SUBVOLUMES, we assign the position
  //  to one of the subHits, preferentially 0x0.
  int ihs[5] = {-1,-1,-1,-1,-1};
  int ic; unsigned int hPat; double eDep; for (ic = 0, eDep = 0, hPat = 0;
					       ic<(int)coalesced.size(); ic++) {
    int ih = coalesced[ic]; SimTrackerHitData &hit = hits[idet]->at(ih);
    int subVolID = hit.cellID>>28&0xf;
    if (subVolID>4) {
      unsigned long cID = hit.cellID;
      printf("Evt #%5d Warning: det %d: SimHit %d(0x%08lx,0x%08lx) has invalid subVolID = %d\n",
	     evtNum,idet,ih,cID&0xffffffff,cID>>32,subVolID);
      continue;
    }
    ihs[subVolID] = ih; hPat |= 0x1<<subVolID;
    eDep += hit.eDep;
  }
  if ((hPat&0x18)==0x18) { // ***** COALESCE SUBVOLUMES 0x3 anD 0x4
    // - Set position = 1/2sum of hits extrapolated by 1/2-pathLength (
    //  disregarding thin SUBVOLUMES, which are still accounted for in "eDep").
    int ih3 = ihs[3], ih4 = ihs[4];
    SimTrackerHitData &h3t = hits[idet]->at(ih3), &h4t = hits[idet]->at(ih4);
    const Vector3d &pos = h3t.position, &pqs = h4t.position;
    double Mx = pos.x, My = pos.y, Mz = pos.z;
    double Nx = pqs.x, Ny = pqs.y, Nz = pqs.z;
    const Vector3f &mom = h3t.momentum, &mqm = h4t.momentum;
    double Px = mom.x, Py = mom.y, Pz = mom.z, P = sqrt(Px*Px+Py*Py+Pz*Pz);
    double Qx = mqm.x, Qy = mqm.y, Qz = mqm.z, Q = sqrt(Qx*Qx+Qy*Qy+Qz*Qz);
    // Extrapolate to extrema
    double Exs[2][2], Eys[2][2], Ezs[2][2]; // [for/backward][3,4]
    int fb; for (fb = 0; fb<2; fb++) {
      int s = fb ? -1 : +1;
      for (int i34 = 0; i34<2; i34++) {
	double &Ex = Exs[fb][i34], &Ey = Eys[fb][i34], &Ez = Ezs[fb][i34];
	if (i34==0) {
	  double u = s*h3t.pathLength/2/P;
	  Ex = Mx+u*Px; Ey = My+u*Py; Ez = Mz+u*Pz;
	} else {
	  double u = s*h4t.pathLength/2/Q;
	  Ex = Nx+u*Qx; Ey = Ny+u*Qy; Ez = Nz+u*Qz;
	}
      }
    }
    // Which pair of extrema do we retain? largest distance
    double d2s[2]; for (fb = 0; fb<2; fb++) {
      int bf = 1-fb;
      double dx = Exs[fb][0]-Exs[bf][1], dy = Eys[fb][0]-Eys[bf][1], dz = Ezs[fb][0]-Ezs[bf][1];
      d2s[fb] = dx*dx+dy*dy+dz*dz;
    }
    fb = d2s[0]>d2s[1] ? 0 : 1; int bf = 1-fb;
    // Position of coalesced = half-sum of extrema
    double Ex = Exs[fb][0], Ey = Eys[fb][0], Ez = Ezs[fb][0];
    double Fx = Exs[bf][1], Fy = Eys[bf][1], Fz = Ezs[bf][1];
    Vector3d pext; pext.x = (Ex+Fx)/2; pext.y = (Ey+Fy)/2; pext.z = (Ez+Fz)/2;
    hext.position = pext; hext.eDep = eDep;
    hext.cellID = h3t.cellID;
    hext.quality = h3t.quality; // Coalesced hits passed samePMO => same quality
    if (doDebug)
      printf("%d(0x%08lx,0x%08lx,%.2f,%.2f,%.2f) + %d(0x%08lx,0x%08lx,%.2f,%.2f,%.2f) -> %.2f,%.2f,%.2f",
	     ih3,h3t.cellID&0xffffffff,h3t.cellID>>32,Mx,My,Mz,
	     ih4,h4t.cellID&0xffffffff,h4t.cellID>>32,Nx,Ny,Nz,
	     pext.x,pext.y,pext.z);
  }
  else if (hPat&0x18) {
    // One of SUBVOLUMES 0x3 or 0x4: extend corresponding subHit (disregarding
    // again thin stuff, which again is nevertheless accounted for in "eDep").
    int ih = (hPat&0x8) ? ihs[3] : ihs[4];
    SimTrackerHitData &hit = hits[idet]->at(ih);
    extend(idet,ih,hext);
    hext.eDep = eDep;
    if (doDebug) {
      unsigned long cID = hit.cellID; const Vector3d &pos = hit.position;
      double Mx = pos.x, My = pos.y, Mz = pos.z;
      printf("%d(0x%08lx,0x%08lx,%.2f,%.2f,%.2f)",
	     ih,cID&0xffffffff,cID>>32,Mx,My,Mz);
    }
  }
  else if (hPat&0x7) {
    int ih; if (hPat&0x1) ih = ihs[0];
    else    if (hPat&0x2) ih = ihs[1];
    else                  ih = ihs[2];
    SimTrackerHitData &hit = hits[idet]->at(ih);
    hext.position = hit.position; hext.eDep = eDep;
    hext.cellID = hit.cellID; hext.quality = hit.quality;
    if (doDebug) {
      unsigned long cID = hit.cellID; const Vector3d &pos = hit.position;
      double Mx = pos.x, My = pos.y, Mz = pos.z;
      printf("%d(0x%08lx,0x%08lx,%.2f,%.2f,%.2f)",
	     ih,cID&0xffffffff,cID>>32,Mx,My,Mz);
    }
  }
  else {
    printf("#%5d: Inconsistency: pattern of hits to be coalesced = 0x%x\n",
	   evtNum,hPat);
    for (int ic = 0; ic<(int)coalesced.size(); ic++) {
      int ih = coalesced[ic]; SimTrackerHitData &hit = hits[idet]->at(ih);
      printf("%d: 0x%08lx,0x%08lx %.2f,%.2f,%.2f\n",
	     ih,hit.cellID&0xffffffff,hit.cellID>>32,
	     hit.position.x,hit.position.y,hit.position.z);
    }
    return false;
  }
  if (doDebug) {
    unsigned long cID = hext.cellID;
    unsigned int module, div, strip; parseCellID(idet,cID,module,div,strip);
    const Vector3d &pos = hext.position;
    double Mx = pos.x, My = pos.y, Mz = pos.z;
    if (idet==0) {
      double Xr, Yr, Zr, Rr, phir;
      g2lCyMBaL(Mx,My,Mz,div,Xr,Yr,Zr,Rr,phir);
      printf(" <-> %.2f,%.2f,%.2f Rr=%.3f\n",Xr,Yr,Zr,Rr);
    } else if (idet==1) {
      double Rcdphi, Xr, Yr, Zr, Ur, Vr;
      g2lOuter(Mx,My,Mz,module,Rcdphi,Xr,Yr,Zr,Ur,Vr);
      printf(" <-> %.2f,%.2f,%.2f Zr=%.3f\n",Xr,Yr,Zr,Zr);
    }
  }
  return true;
}
void recoEvents::extend(int idet, int ih, SimTrackerHitData &hext)
{
  // Extend along P by pathLength/2.
  // - So that hit be sitting close to mid-plane of overall sensitive volume (
  //  this, for a 0x(3|)4 hit; not for a 0x(0|1|2) one, already close enough).
  // - Expected to be called for hit <ih> being...
  //  ...primary,
  //  ...w/o any counterpart.
  //   This corresponds to a particle depositing energy in only one SUBVOLUME.
  // - Extend it, i.e. extrapolate it to the end of its pathLength, to mimic
  //  MPGDTrackerDigi (which in turn tries and emulates what we would get w/ a
  //  single sensitive volume).
  //   It's an approximation: extrapolation should go up to REFERENCE SUBVOLUME.
  //  But it's a good enough approximation.
  // - There are yet cases where extension should not be performed:
  //   + Primary + secondary (delta-ray) w/in SUBVOLUME, see comment in
  //  "extrapolate". => Dependence upon primary vs. secondary ordering. => If
  //  primary hit is sitting on the outer side w.r.t. mid-plane, do not extend.
  //   + Exit through the edge: same thing.
  SimTrackerHitData &hit = hits[idet]->at(ih); unsigned long cID = hit.cellID;
  const Vector3d &pos = hit.position; double Mx = pos.x, My = pos.y, Mz = pos.z;
  const Vector3f &mom = hit.momentum; double Px = mom.x, Py = mom.y, Pz = mom.z;
  double P = sqrt(Px*Px+Py*Py+Pz*Pz), Ex, Ey, Ez;
  if ((verbose&0x100<<idet) || evtNum==evtToDebug) doDebug |= 0x100;
  if  (doDebug&0x100) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%5d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
  }
  bool doExtend = true; int direction = +1;
  int subVolID = cID>>28&0xf; if (subVolID==0x3 || subVolID==0x4) {
    // Determine pathDepth
    if (idet==0) {
      // Cylinder: pathDepth is projected along cylinder's radius
      double pathDepth, Rr; doExtend = checkTraversing(idet,hit,pathDepth,Rr);
      if ((doDebug&0x100) && doExtend)
	printf(" 0x%08lx,0x%08lx Traversing? %d: pathDepth %.3f Rr %.3f\n",
	       cID&0xffffffff,cID>>32,doExtend,pathDepth,Rr);
      if (!doExtend) {
	// Try and rescue cases where traversing to one side.
	int module = cID>>12&0xfff;
	int section = module/8; int staveType = (section==1 || section==2)?0:1;
	double midPlane = radii[0][staveType];
	bool isInner = subVolID==0x3 ? Rr > midPlane : Rr < midPlane;
	doExtend = isInner;
	if (isInner && subVolID==0x4) {
	  // Is it special case of reEntrance w/ hit position outside SUBVOLUME?
	  // As of 2016/01, typically, MPGDTrackerDigi grants this case no
	  // continuation (because of its erratic pathLength, which fails to
	  // get the hit to reach the ( lower) wall precisely enough for the
	  // crossing point to be validated in "cTraversing". Here we don't
	  // have the position+/-pathLength/2==wall kind of check. Therefore we
	  // get into a situation where we will have to extendHit starting from
	  // outside. Which situation triggers an inconsistency warning. In
	  // order to avoid the warning, let's cancel "extendHit" if hit
	  // position lies outside, keeping in mind that the outcome will be the
	  // same: no entension, reproducing MPGDTrackerDigi's stance.
	  double rMin = midPlane+volumeThicknesses[0]/2-radiatorThicknesses[0];
	  if (Rr<rMin) doExtend = false;
	}
	if (doDebug&0x100)
	  printf(" 0x%08lx,0x%08lx Traversing? %d: pathDepth %.3f Rr %.3f/%.3f\n",
		 cID&0xffffffff,cID>>32,doExtend,pathDepth,Rr,midPlane);
      }
      direction = pathDepth>0 ? +1 : -1;
    }
    else if (idet==1) {
      // Cylinder: pathDepth is projected along cylinder's radius
      double pathDepth, Rr; doExtend = checkTraversing(idet,hit,pathDepth,Rr);
      if ((doDebug&0x100)) // && doExtend
	printf(" 0x%08lx,0x%08lx Traversing? %d: pathDepth %.3f Zr %.3f %.3f\n",
	       cID&0xffffffff,cID>>32,doExtend,pathDepth,Rr,hit.pathLength);
    }
    if (doExtend) {
      double gext[3]; doExtend = extendHit(idet,ih,direction,gext);
      if (doExtend) { Ex = gext[0]; Ey = gext[1]; Ez = gext[2]; }
    }
  }
  if (!doExtend) {
    Ex = Mx; Ey = My; Ez = Mz;
  }
  Vector3d pext; pext.x = Ex; pext.y = Ey; pext.z = Ez;
  hext.position = pext; hext.eDep = hit.eDep;
  hext.cellID = hit.cellID; hext.quality = hit.quality;
  if (doDebug&0x100) {
    printf("%d extend\n",ih);
    if (idet==0) {
      unsigned int module, div, strip; parseCellID(idet,hit.cellID,module,div,strip);
      double Xr, Yr, Zr, Rr, phir;
      g2lCyMBaL(Mx,My,Mz,div,Xr,Yr,Zr,Rr,phir); double MR = Rr;
      g2lCyMBaL(Ex,Ey,Ez,div,Xr,Yr,Zr,Rr,phir); double ER = Rr;
      printf(" 0x%08lx 0x%08lx %.2f,%.2f,%.2f (%.3f) -> %.2f,%.2f,%.2f (%.3f) mm\n",
	     hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,MR,Ex,Ey,Ez,ER);
    }
    else
      printf(" 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f mm\n",
	     hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,Ex,Ey,Ez);
      
  }
}
bool recoEvents::checkTraversing(int idet, SimTrackerHitData &hit,
				 double &pathDepth, double &depth)
{
  bool traversing = true;
  const Vector3d &pos = hit.position;
  double Mx = pos.x, My = pos.y, Mz = pos.z;
  const Vector3f &mom = hit.momentum;
  double Px = mom.x, Py = mom.y, Pz = mom.z, P = sqrt(Px*Px+Py*Py+Pz*Pz);
  unsigned int module, div, strip; parseCellID(idet,hit.cellID,module,div,strip);
  int subVolID = hit.cellID>>28&0xf;
  if (idet==0  && (subVolID==0x3 || subVolID==0x4)) {
    double Xr, Yr, Zr, Rr, phir;
    g2lCyMBaL(Mx,My,Mz,div,Xr,Yr,Zr,Rr,phir);
    double lmom[3]; g2lCyMBaL(Px,Py,Pz,div,lmom);
    double Pxr = lmom[0], Pyr = lmom[1];
    double cTheta = (Pxr*Xr+Pyr*Yr)/P/sqrt(Xr*Xr+Yr*Yr);
    pathDepth = hit.pathLength*cTheta; depth = Rr;
    double thickness = radiatorThicknesses[idet];
    traversing = fabs(fabs(pathDepth)-thickness)<.075; // Somewhat large 75 µm tolerance
  }
  else if (idet==1) {
    double Rcdphi, Xr, Yr, Zr, Ur, Vr;
    g2lOuter(Mx,My,Mz,module,Rcdphi,Xr,Yr,Zr,Ur,Vr);
    double lmom[3]; g2lOuter(Px,Py,Pz,module,lmom);
    double Pxr = lmom[2];
    double cTheta = Pxr*abs(Zr)/P/sqrt(Zr*Zr);
    pathDepth = hit.pathLength*cTheta; depth = Zr;
    double thickness = radiatorThicknesses[idet];
    traversing = fabs(fabs(pathDepth)-thickness)<.075; // Somewhat large 75 µm tolerance
  }
  return traversing;
}
// ***** EVENT CONTROL, DEBUGGING
void recoEvents::initDetEvent()
{
  prvPat = 0;
  modPats = 0; // Module-pattern of all modules w/ hit
  // Modules where LONE PRIMARIES and only that
  modulesLP = 0; moduleLP = 0;
  modulesLP1 = 0; // In addition, single hit
  // debugging
  doDebug = 0;
}
void recoEvents::updateDetEvent(int idet, int ih)
{
  // LONE PRIMARY
  // - Is primary hit w/o secondary offsprings.
  // - I.e. primary not followed by secondary.
  // - Identified by unbroken sequence of primary hits (in a broken sequence,
  //  some of the hits may be LONE PRIMARIES, but we're interested here in
  //  singling out modules where all hits are lone primaries).
  //   This, PROVIDED PRIMARY IS REQUIRED.
  SimTrackerHitData &hit = hits[idet]->at(ih);
  int module = (hit.cellID>>12)&0xfff;
  unsigned long modPat = ((unsigned long)0x1)<<module;
  if (modPat!=prvPat) { // New module: update module-patterns w/ previous one
    modulesLP |= moduleLP; modulesLP1 |= moduleLP;
  }
  if (prvPat&modPats)   // Already encountered module...
    modulesLP1 &= ~prvPat; // ...Remove from module-pattern of SINGLE LONE PRIMARY
  // LONE PRIMARY?
  unsigned int status = getStatus(idet,ih);
  if (status&0x2) {
    if (modPat!=prvPat) // New module has primary...
      moduleLP = modPat;  // ...=> Set "moduleLP"
  }
  else                // Interfering secondary
    moduleLP = 0;        // ...=> Reset "moduleLP"
  modPats |= prvPat;
  prvPat = modPat;
}
void recoEvents::finaliseDetEvent()
{ // Last iteration: update module-patterns w/ last module
  modulesLP |= moduleLP; modulesLP1 |= moduleLP;
  if (prvPat&modPats) modulesLP1 &= ~prvPat;
}
bool recoEvents::moduleSelection(unsigned long cellID)
{
  int module = cellID>>12&0xfff;
  unsigned long modPat = ((unsigned long)0x1)<<module;
  if (requireModules && !(modPat&requireModules)) return false;
  bool ok = true;
  if (requireQuality) {
    ok = modPat&modulesLP;
    if (requireQuality>1) ok = modPat&modulesLP1;
  }
  return ok;
}
void recoEvents::debugHit(int idet, int ih, int nHs, SimTrackerHitData &hit, unsigned int status)
{
  if ((verbose&0x1<<idet) || evtNum==evtToDebug) doDebug |= 0x1;
  if (doDebug&0x1) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%5d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
    const Vector3d &pos = hit.position;
    printf("hit %2d/%-2d: 0x%08lx,0x%08lx X,Y,Z %7.2f,%7.2f,%8.2f status 0x%x\n",
	   ih,nHs,hit.cellID&0xffffffff,hit.cellID>>32,
	   pos.x,pos.y,pos.z,status);
  }
}
void recoEvents::debugHit(int idet, SimTrackerHitData &hit)
{
  if (doDebug&0x1) {
    const Vector3d &pos = hit.position;
    double X = pos.x, Y = pos.y, Z = pos.z;
    printHit(idet,X,Y,Z,hit.cellID);
  }
}
void recoEvents::printHit(int idet,
			  double X, double Y, double Z, unsigned long cellID)
{
  // Print: global, then local
  double R2 = X*X+Y*Y, R = sqrt(R2); double phi = atan2(Y,X);
  const double pi = TMath::Pi();
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  unsigned int cell = cellID>>32;
  if (verbose&0x10000) {
    printf(" 0x%08lx,0x%08x  %7.2f,%7.2f,%8.2f %7.2f cm %6.3fπ\n",
	   cellID&0xfffffff,cell,X,Y,Z,R,phi/pi);
    if      (idet==1) printf(" %5.1f 0x%02x 0x%x",phi/pi*12,module,module>>1);
    else if (idet==0) printf(" %5.1f 0x%02x 0x%x",phi/pi*8, module,module>>1);
  }
  else if (idet==2) {
    unsigned int layer = module&0xf, tile = module>>4;
    printf(" %d,%2d,0x%08x  %7.2f,%7.2f,%8.2f %7.2f cm %6.3fπ\n",
	   layer,tile,cell,X,Y,Z,R,phi/pi);
  }
  else if (idet==0) {
    double Xr, Yr, Zr, Rr, phir; g2lCyMBaL(X,Y,Z,div,Xr,Yr,Zr,Rr,phir);
    printf("%10s 0x%08lx,0x%08x X,Y,Z %7.2f,%7.2f,%8.2f Rr %7.3f mm phir %6.3fπ\n",
	   "",cellID&0xffffffff,cell,X,Y,Z,Rr,phir/pi);
  }
  else if (idet==1) {
    double Rcdphi, Xr, Yr, Ur, Vr, Zr; g2lOuter(X,Y,Z,module,Rcdphi,Xr,Yr,Zr,Ur,Vr);
    printf("%10s 0x%08lx,0x%08x X,Y,Z %7.2f,%7.2f,%8.2f Rr %7.3f U,V %7.2f,%7.2f mm\n",
	   "",cellID&0xffffffff,cell,X,Y,Z,Rcdphi,Ur,Vr);
  }
  else {
    printf("** recoEvents::getDetHit: Invalid det# = %d\n",idet);
    exit(1);
  }
}
void recoEvents::debugRec(int idet, int ir)
{
  if (!(verbose&0x1<<idet)) return;
  int nRecs = recs[idet]->size();
  edm4eic::TrackerHitData &rec = recs[idet]->at(ir);
  const Vector3f &pos = rec.position;
  printf("rec %d/%d: 0x%08lx,0x%08lx %6.1f,%6.1f,%6.1f\n",
	 ir,nRecs,rec.cellID&0xffffffff,rec.cellID>>32,pos.x,pos.y,pos.z);
}
void recoEvents::debugHitRec(int idet, int is, int ncoas, int cIndex, SimTrackerHitData &hit, int ir)
{
  if ((verbose&0x10<<idet) || evtNum==evtToDebug) {
    int nRecs = recs[idet]->size();
    edm4eic::TrackerHitData &rec = recs[idet]->at(ir);
    printf("SimHit %d/%d:%d(0x%08lx,0x%08lx) RecHit %d/%d 0x%08lx,0x%08lx\n",
	   is,ncoas,cIndex,hit.cellID&0xffffffff,hit.cellID>>32,
	   ir,nRecs,rec.cellID&0xffffffff,rec.cellID>>32);
  }
}
void recoEvents::debugAssoc(int idet)
{
  if (verbose&0x10<<idet) doDebug |= 0x10;
  if (doDebug&0x10) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%5d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
    int nArhs = arhs[idet]->size(); // Associations raw hit
    int nAshs = ashs[idet]->size(); // Associations sim hit
    int nARhs = aRhs[idet]->size(); // Associations Rec Hit -> raw hit
    printf("arhs %d\n",nArhs);
    for (int ih = 0; ih<nArhs; ih++) {
      podio::ObjectID &arh = arhs[idet]->at(ih);
      int rIndex = arh.index; unsigned int cID = arh.collectionID;
      printf("%d: %d 0x%x\n",ih,rIndex,cID);
    }
    printf("ashs %d\n",nAshs);
    for (int ih = 0; ih<nAshs; ih++) {
      podio::ObjectID &ash = ashs[idet]->at(ih);
      int sIndex = ash.index; unsigned int cID = ash.collectionID;
      printf("%d: %d 0x%x\n",ih,sIndex,cID);
    }
    printf("aRhs %d\n",nARhs);
    for (int iR = 0; iR<nARhs; iR++) {
      podio::ObjectID &aRh = aRhs[idet]->at(iR);
      printf("%d: %d\n",iR,aRh.index);
    }
  }
}
void recoEvents::debugAssoc(int idet, map<int,int>raw2rec, map<int,int> sim2coa, map<int,vector<int>> rec2coas)
{
  if (!(verbose&0x10<<idet)) return;
  map<int,int>::const_iterator ir;
  for (ir = raw2rec.cbegin(); ir != raw2rec.cend(); ir++) {
    printf("raw2rec %d -> %d\n",ir->first,ir->second);
  }
  for (ir = sim2coa.cbegin(); ir != sim2coa.cend(); ir++) {
    printf("sim2coa %d -> %d\n",ir->first,ir->second);
  }
  map<int,vector<int>>::const_iterator im;
  for (im = rec2coas.cbegin(); im != rec2coas.cend(); im++) {
    printf("rec2coas: %d ->",im->first);
    const vector<int> &coas = im->second;
    for (int is = 0; is<(int)coas.size(); is++) printf(" %d",coas[is]);
    printf("\n");
  }
}
// ***** EXTENSION
bool recoEvents::extendHit(int idet, int ih, int direction, double *gext)
{
  // - Extend only SUBVOLUMES 0x3 and 0x4 (other ones are presumed to be close
  //  enough to mid-plane.
  // - Mid-way between:
  //   - Extrapolation by pathLengh/2 opposite to counterpart,
  //   - Extreme of counterpart.
  SimTrackerHitData &hit = hits[idet]->at(ih); unsigned long cID = hit.cellID;
  const Vector3d &pos = hit.position; double Mx = pos.x, My = pos.y, Mz = pos.z;
  const Vector3f &mom = hit.momentum; double Px = mom.x, Py = mom.y, Pz = mom.z;
  double P = sqrt(Px*Px+Py*Py+Pz*Pz);
  int subVolID = cID>>28&0xf;
  // Determine ini extremum = extrememum facing away from conterpart SUBVOLUME.
  double gini[3]; // Global ini position
  double &Ax = gini[0], &Ay = gini[1], &Az = gini[2], u;
  if (subVolID==3)       u = -direction*hit.pathLength/2/P;
  else /* subVolID==4 */ u =  direction*hit.pathLength/2/P;
  Ax = Mx+u*Px; Ay = My+u*Py; Az = Mz+u*Pz;
  // Determine end extremum = extremum into counterpart SUBVOLUME.
  double gend[3]; // Global end position
  // Transform to local where limitations to extension are conveniently defined.
  // If CyMBaL, limitations in (phir,Zr); if Outer, limitations in (Xr,Yr).
  double lpos[3], &Xr = lpos[0], &Yr = lpos[1], &Zr = lpos[2];
  if (idet==0) {
    unsigned int module, div, strip; parseCellID(idet,cID,module,div,strip);
    double Rr, phir; g2lCyMBaL(Mx,My,Mz,div,Xr,Yr,Zr,Rr,phir);
    double lmom[3]; g2lCyMBaL(Px,Py,Pz,div,lmom);
    // Limitations
    int section = div>>3, staveType = section==1 || section==2 ? 0 : 1;
    double endPhi = hWidths[0][staveType], startPhi = -endPhi;
    double dZ = ZHLengths[0];
    // Target = extreme edge of counterpart SUBVOLUME.
    int pm = subVolID==3 ? +1 : -1;
    double rT = radii[0][staveType]+pm*volumeThicknesses[0]/2;
    // Local lend
    double lend[3];
    unsigned int status =
      cExtension(lpos,lmom,rT,pm*direction,dZ,startPhi,endPhi,lend);
    if (!(status&0x1)) {
      printf("#%5d: extendHit(%d,%d(0x%08lx,0x%08lx)) %.4f-(%d)->%.4f inconsistency: 0x%x\n",
	     evtNum,idet,ih,cID&0xffffffff,cID>>32,Rr,pm*direction,rT,status);
      printf("%.2f,%.2f(%.4f),%.2f mm %.4f,%.4f,%.4f GeV\n",
	     Xr,Yr,Rr,Zr,lmom[0],lmom[1],lmom[2]);
      return false;
    }
    else {
      double gend[3]; l2gCyMBaL(lend,div,gend);
      for (int i = 0; i<3; i++) gext[i] = (gini[i]+gend[i])/2;
    }
  }
  else {
    int pm = subVolID==3 ? +1 : -1; double u = pm*direction*hit.pathLength/2/P;
    gext[0] = Mx+u*Px; gext[1] = My+u*Py; gext[2] = Mz+u*Pz;
  }
  return true;
}
unsigned int cExtension(double const* lpos, double const* lmom, // Input subHit
			double rT,                              // Target radius
			int direction, double dZ, double startPhi,
			double endPhi, // Module parameters
			double* lext)
{
  unsigned int status = 0;
  double Mx = lpos[0], My = lpos[1], Mz = lpos[2];
  double Px = lmom[0], Py = lmom[1], Pz = lmom[2], norm = sqrt(Px * Px + Py * Py + Pz * Pz);
  double M2 = Mx * Mx + My * My, rIni = sqrt(M2), rLow, rUp;
  if (rIni < rT) {
    rLow = rIni;
    rUp  = rT;
  } else {
    rLow = rT;
    rUp  = rIni;
  }
  // Intersection w/ the edge in phi
  double tF = 0;
  for (double phi : {startPhi, endPhi}) {
    // M+t*P = 0 + t'*U. t = (My*Ux-Mx*Uy)/(Px*Uy-Py*Ux);
    double Ux = cos(phi), Uy = sin(phi);
    double D = Px * Uy - Py * Ux;
    if (D) { // If P not // to U
      double t = (My * Ux - Mx * Uy) / D;
      if (t * direction < 0)
        continue;
      double Ex = Mx + t * Px, Ey = My + t * Py, rE = sqrt(Ex * Ex + Ey * Ey), Ez = Mz + t * Pz;
      if (rLow < rE && rE < rUp && fabs(Ez) < dZ) {
        status |= 0x1;
        tF = t;
      }
    }
  }
  // Intersection w/ the edge in Z
  double zLow = -dZ, zUp = +dZ;
  for (double Z : {zLow, zUp}) {
    // Mz+t*Pz = Z
    if (Pz) {
      double t = (Z - Mz) / Pz;
      if (t * direction < 0)
        continue;
      double Ex = Mx + t * Px, Ey = My + t * Py, rE = sqrt(Ex * Ex + Ey * Ey);
      double phi = atan2(Ey, Ex);
      if (rLow < rE && rE < rUp && startPhi < phi && phi < endPhi) {
        if (t < 0) {
          if (!status || (status && t > tF)) {
            status |= 0x1;
            tF = t;
          }
        } else if (t > 0) {
          if (!status || (status && t < tF)) {
            status |= 0x1;
            tF = t;
          }
        }
      }
    }
  }
  // Else intersection w/ target radius
  if (!status) {
    double a = Px * Px + Py * Py, b = Px * Mx + Py * My, c = M2 - rT * rT;
    if (!a) {           // P is // to Z (while it did no intersect the edge in Z)
      status |= 0x1000; // Inconsistency
    } else if (!c) {    // Hit is on target (while we've moved away from it)
      status |= 0x2000; // Inconsistency
    } else {
      double det = b * b - a * c;
      if (det >= 0) {
        double sqdet = sqrt(det);
        for (int is = 0; is < 2; is++) {
          int s    = 1 - 2 * is;
          double t = (-b + s * sqdet) / a;
          if (t * direction < 0)
            continue;
          double Ix = Mx + t * Px, Iy = My + t * Py, Iz = Mz + t * Pz, phi = atan2(Iy, Ix);
          if (fabs(Iz) > dZ || phi < startPhi || endPhi < phi)
            continue;
          if (!(status & 0x1) ||
              // Two intersects: let's retain the earliest one.
              ((status & 0x1) && fabs(t) < fabs(tF))) {
            tF = t;
            status |= 0x1;
          }
        }
      }
    }
  }
  if (status&0x1) {
    lext[0] = Mx + tF * Px;
    lext[1] = My + tF * Py;
    lext[2] = Mz + tF * Pz;
  }
  return status;
}
void cDrawHit(double *pars, double *lpos, double *lmom, double path, double *r2s);
void bDrawHit(double *pars, double *lpos, double *lmom, double path, double *r2s);
void AddHit(double *lpos, double *lmom, double path);
void recoEvents::DrawSimHit(int jentry, unsigned int detectorPattern, int ih,
			    bool addToPreExisting)
{
  // ***** PARSE ARGs
  if (!fChain->GetEntry(jentry)) {
    printf("** DrawSimHit: Invalid <jentry> arg. = %d\n",jentry);
    return;
  }
  unsigned int pat = processedDetectors&detectorPattern;
  if (!pat) {
    printf("** DrawSimHit: None of the processed detectors(=0x%x) in arg. <detectorPattern>(=0x%x)\n",
	   processedDetectors,detectorPattern);
    return;
  }
  // ***** LOOP ON DETECTORS
  for (int idet = 0; idet<N_DETs; idet++) {
    if (!(0x1<<idet&pat)) continue;
    int nHits = hits[idet]->size(); // #SimHits
    if (ih>=nHits) {
      printf("** DrawSimHit: Invalid <ih>(=%d) arg., >= #SimHits[%d](=%d)\n",
	     ih,idet,nHits);
      return;
    }
    SimTrackerHitData &hit = hits[idet]->at(ih); unsigned long cID = hit.cellID;
    // Globals
    const Vector3d &pos = hit.position; double Mx = pos.x, My = pos.y, Mz = pos.z;
    const Vector3f &mom = hit.momentum; double Px = mom.x, Py = mom.y, Pz = mom.z;
    // Locals
    double Xr, Yr, Zr, lpos[3];
    if (idet==0) {
      // Globals -> Locals
      unsigned int module, div, strip; parseCellID(idet,cID,module,div,strip);
      double Rr, phir, rot; g2lCyMBaL(Mx,My,Mz,div,Xr,Yr,Zr,Rr,phir);
      lpos[0] = Xr/10; lpos[1] = Yr/10; lpos[2] = Zr/10;
      double lmom[3]; g2lCyMBaL(Px,Py,Pz,div,lmom);
      // Pathlength
      double path = hit.pathLength/10;
      // Parameters
      int section = div>>3, staveType = section==1 || section==2 ? 0 : 1;
      double r0 = radii[0][staveType], rMin, rMax;
      rMin = r0-volumeThicknesses[0]/2;
      rMax = r0-volumeThicknesses[0]/2+radiatorThicknesses[0];
      rMin /= 10; rMax /= 10;
      double dZ = ZHLengths[0]/10;
      double endPhi = hWidths[0][staveType], startPhi = -endPhi;
      double pars[5] = {rMin,rMax,startPhi,endPhi,dZ};
      rMin = r0+volumeThicknesses[0]/2-radiatorThicknesses[0];
      rMax = r0+volumeThicknesses[0]/2;
      rMin /= 10; rMax /= 10;
      double r2s[2] = {rMin,rMax};
      
      printf("SimHit 0x%08lx,0x%08lx\n",cID&0xffffffff,cID>>32);
      printf("double lpos[3] = {%.6f,%.6f,%.6f};\n",
	     lpos[0]/10,lpos[1]/10,lpos[2]/10);
      printf("double lmom[3] = {%.6f,%.6f,%.6f};\n",lmom[0],lmom[1],lmom[2]);
      printf("double path = %.6f;\n",hit.pathLength/10);
      printf("double pars[5] = {%.5f,%.5f,%.5f,%.5f,%.5f};\n",
	rMin,rMax,startPhi,endPhi,dZ);
      printf("double r2s[2] = {%.5f,%.5f};\n",r2s[0],r2s[1]);

      // DrawHit
      if (addToPreExisting) AddHit(       lpos,lmom,path);
      else                  cDrawHit(pars,lpos,lmom,path,r2s);
    }
    else if (idet==1) {
      // Globals -> Locals
      unsigned int module, div, strip; parseCellID(idet,cID,module,div,strip);
      double Rcdphi, Ur, Vr; g2lOuter(Mx,My,Mz,module,Rcdphi,Xr,Yr,Zr,Ur,Vr);
      lpos[0] = Zr/10; lpos[1] = Xr/10; lpos[2] = Yr/10;
      double lmom[3]; g2lOuter(Px,Py,Pz,module,lmom);
      // Pathlength
      double path = hit.pathLength/10;
      // Parameters
      double r0 = 0, rMin, rMax;
      rMin = r0-volumeThicknesses[1]/2;
      rMax = r0-volumeThicknesses[1]/2+radiatorThicknesses[1];
      rMin /= 10; rMax /= 10;
      double dZ = ZHLengths[0]/10;
      double endX = hWidths[1][0]/10, startX = -endX;
      double pars[5] = {rMin,rMax,startX,endX,dZ};
      rMin = r0+volumeThicknesses[1]/2-radiatorThicknesses[1];
      rMax = r0+volumeThicknesses[1]/2;
      rMin /= 10; rMax /= 10;
      double r2s[2] = {rMin,rMax};

      printf("gpos %.6f,%.6f,%.6f\n",Mx/10,My/10,Mz/10);      
      printf("gmom %.6f,%.6f,%.6f\n",Px,Py,Pz);

      printf("double lpos[3] = {%.6f,%.6f,%.6f};\n",
	     lpos[0]/10,lpos[1]/10,lpos[2]/10);
      printf("double lmom[3] = {%.6f,%.6f,%.6f};\n",lmom[0],lmom[1],lmom[2]);
      printf("double path = %.6f;\n",hit.pathLength/10);
      printf("double pars[5] = {%.5f,%.5f,%.5f,%.5f,%.5f};\n",
	rMin,rMax,startX,endX,dZ);
      printf("double r2s[2] = {%.5f,%.5f};\n",r2s[0],r2s[1]);

      // DrawHit
      if (addToPreExisting) AddHit(      lpos,lmom,path);
      else                  bDrawHit(pars,lpos,lmom,path,r2s);
    }
  }
}
#include "TEllipse.h"
#include "TLine.h"
void cDrawHit(double *pars, double *lpos, double *lmom, double path, double *r2s)
{
  double rMin = pars[0], rMax = pars[1], dZ = pars[4]; 
  double startPhi = pars[2], endPhi = pars[3];
  double phMn = startPhi*180/TMath::Pi(), phMx = endPhi*180/TMath::Pi();

  TCanvas *c = new TCanvas("cEvt");

  // TH2
  double xLow = rMin*cos(endPhi), xUp = rMax, dx = (xUp-xLow)/10;
  xLow -= dx; xUp += dx;
  double yUp  = rMax*sin(endPhi), yLow = -yUp, dy = (yUp-yLow)/10; 
  yLow -= dy; yUp += dy;
  TH2D *h2 = new TH2D("hEvt",";X (cm)  ;Y (cm)  ",2000,xLow,xUp,2000,yLow,yUp);
  h2->Draw();
  h2->SetStats(0); h2->SetFillColorAlpha(0,0); h2->SetLineColor(0);
  TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
  ax->SetNdivisions(505); ay->SetNdivisions(505);
  ay->SetTitleOffset(1.3);
  int nBinsX = h2->GetNbinsX(), nBinsY = h2->GetNbinsY();

  // Draw cylinders
  const double ringWidth = 2;
  if (r2s) {
    int orange = kOrange+1, vert = kGreen+2;
    TEllipse *eMin2 = new TEllipse(0,0,r2s[0],r2s[0],phMn,phMx);
    //eMin2->SetLineColor(orange); eMin2->SetLineWidth(ringWidth);
    eMin2->SetLineColor(2); eMin2->SetLineWidth(ringWidth);
    eMin2->SetFillColorAlpha(0,0);
    TEllipse *eMax2 = new TEllipse(0,0,r2s[1],r2s[1],phMn,phMx);
    //eMax2->SetLineColor(vert);   eMax2->SetLineWidth(ringWidth);
    eMax2->SetLineColor(2);   eMax2->SetLineWidth(ringWidth);
    eMax2->SetFillColorAlpha(0,0);
    eMax2->Draw();
    eMin2->Draw();
  }
  TEllipse *eMin = new TEllipse(0,0,rMin,rMin,phMn,phMx);
  //eMin->SetLineColor(2); eMin->SetLineWidth(ringWidth);
  eMin->SetLineColor(4); eMin->SetLineWidth(ringWidth);
  eMin->SetFillColorAlpha(0,0);
  TEllipse *eMax = new TEllipse(0,0,rMax,rMax,phMn,phMx);
  eMax->SetLineColor(4); eMax->SetLineWidth(ringWidth);
  eMax->SetFillColorAlpha(0,0);
  eMax->Draw();
  eMin->Draw();
  
  // Draw <lpos>
  double w = .0125;
  int  orange = kOrange+1;
  TLine *tl = new TLine; tl->SetLineColor(orange);
  double Mx = lpos[0], My = lpos[1];
  tl->DrawLine(Mx-w,My,Mx+w,My); tl->DrawLine(Mx,My-w,Mx,My+w);

  // Draw <lpos>+/-path/2
  // Extrema
  TLine *tp = new TLine; tp->SetLineColor(1);
  double xMn = xUp, xMx = xLow, yMn = yUp, yMx = yLow;
  double Px = lmom[0], Py = lmom[1], Pz = lmom[2];
  double norm = sqrt(Px*Px+Py*Py+Pz*Pz), at = path/2/norm;
#define DEBUG_DRAWHIT
#ifdef DEBUG_DRAWHIT
  printf("x,yLow/Up %.2f,%.2f %.2f,%.2f path/norm => at %.2f,%.2f,%f\n",
	 xLow,xUp,yLow,yUp,path,norm,at);
#endif
  for (int s = -1; s<=+1; s += 2) {
    double Nx = Mx+s*at*Px, Ny = My+s*at*Py;
    tp->DrawLine(Mx,My,Nx,Ny);
    tp->DrawLine(Nx-w,Ny,Nx+w,Ny); tp->DrawLine(Nx,Ny-w,Nx,Ny+w);
    if (Nx-w<xMn) xMn = Nx-w; if (Ny-w<yMn) yMn = Ny-w;
    if (Nx+w>xMx) xMx = Nx+w; if (Ny+w>yMx) yMx = Ny+w;
#ifdef DEBUG_DRAWHIT
    printf("%d: Mx,y %.2f,%.2f/%.2fx,yMn/Mx %.2f,%.2f %.2f,%.2f\n",
	   s,Mx,My,w,xMn,xMx,yMn,yMx);
#endif
  }

  // Re-draw h2
  // (This givess interactive access to its TAxis objects, otherwise hidden by
  // the TEllipse.)
  //  h2->Draw("same");

  // Range
  //dx = (xMx-xMn)/4; xMn -= dx; xMx += dx;
  //dy = (yMx-yMn)/4; yMn -= dy; yMx += dy;
  dx = (xMx-xMn)/1; xMn -= dx; xMx += dx;
  dy = (yMx-yMn)/2; yMn -= dy; yMx += dy;
  // If pathLength short, range may be too narrow: expand it, that it's >> 0.3cm.
  if (xMx-xMn<.3) {
    double xMean = (xMn+xMx)/2; xMn = xMean-.4; xMx = xMean+.4;
  }
  if (yMx-yMn<.3) {
    double yMean = (yMn+yMx)/2; yMn = yMean-.4; yMx = yMean+.4;
  }
  int bxMn = int(1+nBinsX*(xMn-xLow)/(xUp-xLow)+.5), bxMx = int(1+nBinsX*(xMx-xLow)/(xUp-xLow)+.5);
  int byMn = int(1+nBinsY*(yMn-yLow)/(yUp-yLow)+.5), byMx = int(1+nBinsY*(yMx-yLow)/(yUp-yLow)+.5);
  printf("%.2f,%.2f  %.2f,%.2f  [%d,%d] [%d,%d]\n",
	 xMn,xMx,yMn,yMx,bxMn,bxMx,byMn,byMx);
  ax->SetRange(bxMn,bxMx);
  ay->SetRange(byMn,byMx);
  //ax->Draw(); ay->Draw();
  
  gPad->Modified(); gPad->Update();
}
void bDrawHit(double *pars, double *lpos, double *lmom, double path, double *r2s)
{
  double rMin = pars[0], rMax = pars[1], dZ = pars[4];
  double startX = pars[2], endX = pars[3];
  double XMn = startX, XMx = endX;
  double rMi2 = r2s[0],  rMa2 = r2s[1];

  TCanvas *c = new TCanvas("bEvt");

  // TH2
  double xLow = rMin, xUp = rMa2, dx = (xUp-xLow)/10;
  xLow -= dx; xUp += dx;
  double yUp  = XMx, yLow = XMn,  dy = (yUp-yLow)/10; 
  yLow -= dy; yUp += dy;
  TH2D *h2 = new TH2D("hEvt",";Z (cm)  ;X (cm)  ",2000,xLow,xUp,2000,yLow,yUp);
  h2->Draw();
  h2->SetStats(0); h2->SetFillColorAlpha(0,0); h2->SetLineColor(0);
  TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
  ax->SetNdivisions(505); ay->SetNdivisions(505);
  ay->SetTitleOffset(1.3);
  int nBinsX = h2->GetNbinsX(), nBinsY = h2->GetNbinsY();

  // Draw boxes
  const double ringWidth = 2;
  TLine *box1 = new TLine;
  box1->SetLineColor(2); box1->SetLineWidth(ringWidth);
  box1->DrawLine(rMin,XMn,rMin,XMx); box1->DrawLine(rMax,XMn,rMax,XMx);
  TLine *box2 = new TLine;
  box2->SetLineColor(4); box2->SetLineWidth(ringWidth);
  box2->DrawLine(rMi2,XMn,rMi2,XMx); box2->DrawLine(rMa2,XMn,rMa2,XMx);

  // Draw <lpos>
  double w = .05;
  double Mx = lpos[0], My = lpos[1];
  int orange = kOrange+1, vert = kGreen+2, col = Mx<rMax ? vert : orange;
  TLine *tl = new TLine; tl->SetLineColor(orange);
  tl->DrawLine(Mx-w,My,Mx+w,My); tl->DrawLine(Mx,My-w,Mx,My+w);

  // Draw <lpos>+/-path/2
  // Extrema
  TLine *tp = new TLine; tp->SetLineColor(orange);
  double xMn = xUp, xMx = xLow, yMn = yUp, yMx = yLow;
  double Px = lmom[0], Py = lmom[1], Pz = lmom[2];
  double norm = sqrt(Px*Px+Py*Py+Pz*Pz), at = path/2/norm;
  for (int s = -1; s<=+1; s += 2) {
    double Nx = Mx+s*at*Px, Ny = My+s*at*Py;
    tp->DrawLine(Mx,My,Nx,Ny);
    tp->DrawLine(Nx-w,Ny,Nx+w,Ny); tp->DrawLine(Nx,Ny-w,Nx,Ny+w);
    if (Nx-w<xMn) xMn = Nx-w; if (Ny-w<yMn) yMn = Ny-w;
    if (Nx+w>xMx) xMx = Nx+w; if (Ny+w>yMx) yMx = Ny+w;
#ifdef DEBUG_DRAWHIT
    printf("%d: Mx,y %.2f,%.2f/%.2fx,yMn/Mx %.2f,%.2f %.2f,%.2f\n",
	   s,Mx,My,w,xMn,xMx,yMn,yMx);
#endif
  }

  // Re-draw h2
  // (This givess interactive access to its TAxis objects, otherwise hidden by
  // the TEllipse.)
  //  h2->Draw("same");

  // Range
  //dx = (xMx-xMn)/4; xMn -= dx; xMx += dx;
  //dy = (yMx-yMn)/4; yMn -= dy; yMx += dy;
  dx = (xMx-xMn)/1; xMn -= dx; xMx += dx;
  dy = (yMx-yMn)/2; yMn -= dy; yMx += dy;
  // If pathLength short, range may be too narrow: expand it, that it's >> 0.3cm.
  if (xMx-xMn<.3) {
    double xMean = (xMn+xMx)/2; xMn = xMean-.4; xMx = xMean+.4;
  }
  if (yMx-yMn<.3) {
    double yMean = (yMn+yMx)/2; yMn = yMean-.4; yMx = yMean+.4;
  }
  int bxMn = int(1+nBinsX*(xMn-xLow)/(xUp-xLow)+.5), bxMx = int(1+nBinsX*(xMx-xLow)/(xUp-xLow)+.5);
  int byMn = int(1+nBinsY*(yMn-yLow)/(yUp-yLow)+.5), byMx = int(1+nBinsY*(yMx-yLow)/(yUp-yLow)+.5);
  printf("%.2f,%.2f  %.2f,%.2f  [%d,%d] [%d,%d]\n",
	 xMn,xMx,yMn,yMx,bxMn,bxMx,byMn,byMx);
  ax->SetRange(bxMn,bxMx);
  ay->SetRange(byMn,byMx);
  //ax->Draw(); ay->Draw();
  
  gPad->Modified(); gPad->Update();
}
void AddHit(double *lpos, double *lmom, double path)
{
  // Draw <lpos>
  //double w = .0125;
  double w = .025;
  int  vert = kGreen+2;
  TLine *tl = new TLine; tl->SetLineColor(vert);
  double Mx = lpos[0], My = lpos[1];
  tl->DrawLine(Mx-w,My,Mx+w,My); tl->DrawLine(Mx,My-w,Mx,My+w);

  // Draw <lpos>+/-path/2
  TLine *tp = new TLine; tp->SetLineColor(1);
  double Px = lmom[0], Py = lmom[1], Pz = lmom[2];
  double norm = sqrt(Px*Px+Py*Py+Pz*Pz), at = path/2/norm;
  for (int s = -1; s<=+1; s += 2) {
    double Nx = Mx+s*at*Px, Ny = My+s*at*Py;
    tp->DrawLine(Mx,My,Nx,Ny);
    tp->DrawLine(Nx-w,Ny,Nx+w,Ny); tp->DrawLine(Nx,Ny-w,Nx,Ny+w);
  }
}
