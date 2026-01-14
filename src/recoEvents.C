#define recoEvents_cxx
#include "recoEvents.h"
#include <TH2.h>
#include <TVector3.h>
#include <TRotation.h>

#include <stdio.h>

// Globals
unsigned int recoEvents::nObjCreated = 0;

// Interfaces
void getReducedCyMBaL(double X, double Y, unsigned int div,
		      double &Xr, double &Yr, double &Rr, double &phir,
		      int *zone = 0, // Whether in peak, L/R edge, else
		      double *rot = 0);   // Rotation angle
void getReducedOuter(double X, double Y, double Z, unsigned int div,
		     double &Rcphi, double &Xr, double &Yr,
		     double &Ur, double &Vr);

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
	fillHit(0,idet,X,Y,Z,hit.cellID);
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
	printf("Evt #%d Warning: %3d det %d: ash(%d)!=arh(%d)\n",
	       evtNum,(int)jentry,idet,nAshs,nArhs);
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
	  printf("Warning: Entry %3d det %d: rec %d not associated\n",
		 (int)jentry,idet,ir);
	}
	else {
	  const vector<int> &coas = im->second;
	  int is, selecRec; for (is=selecRec = 0; is<(int)coas.size(); is++) {
	    int cIndex = coas[is];
	    if (cIndex<0 || nHits<=cIndex) {
	      printf("#%d Warning: det %d: Rec %d -> sim %d\n",
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
	  if (selecRec) {
	    fillHit(1,idet,X,Y,Z,rec.cellID);     // ***** FILL rec HISTOS
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
			 double X, double Y, double Z, unsigned long cellID)
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
  hs->thphi->Fill(phi,theta);
  hs->XY->Fill(X,Y); hs->ZR->Fill(Z,R);

  if      (idet==0) { // Special CyMBaL: fill reduced Radius
    double Xr, Yr, Rr, phir; int zone; getReducedCyMBaL(X,Y,div,Xr,Yr,Rr,phir,&zone);
    hs->phir->Fill(phir,div);
    hs->Rr->Fill(Rr,div);
    hs->XYr[zone]->Fill(Xr,Yr);
  }
  else if (idet==1) { // Special Outer: fill Rcosphi, Ur, Vr
    double Rcdphi, Xr, Yr, Ur, Vr; int zone; getReducedOuter(X,Y,Z,module,Rcdphi,Xr,Yr,Ur,Vr);
    hs->Rr->Fill(Rcdphi,div);
    hs->Ur->Fill(Ur,div); hs->Vr->Fill(Vr,div);
  }
}
void recoEvents::fillResids(int idet, const Vector3f &pos, const Vector3d &psim, unsigned long cellID)
{
  int doPrint = 0;
  if (verbose&0x1000<<idet) {
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
    double Xd, Yd, Rr, phir;   getReducedCyMBaL(X,Y,div,Xd,Yd,Rr,phir);
    double Rrs, phirs;         getReducedCyMBaL(Xs,Ys,div,Xd,Yd,Rrs,phirs);
    if (doPrint>1)
      printf("Rr,Rrs: %.2f,%.2f\n",Rr,Rrs);
    double dRr = 1000*(Rr-Rrs), dphir = 1000*(phir-phirs);
    rs.Rr->Fill(dRr); rs.phir->Fill(dphir);
    if (doPrint) {
      if (strip==1 && fabs(dZ)>580) {
      //if (strip==1 && fabs(dZ)>0) {
	printf("#%5d 0x%08lx,0x%08lx  %7.2f,%7.2f,%8.2f %7.3f mm %7.3fπ\n",
	       evtNum,cellID&0xffffffff,cellID>>32,X,Y,Z,Rr,phir/TMath::Pi());
	printf("%16s %12s %7.2f,%7.2f,%8.2f %7.3f mm %7.3fπ  dZ %4.0f µm\n",
	       "SimHit","",Xs,Ys,Zs,Rrs,phirs/TMath::Pi(),dZ);
      }
    }
  }
  else if (idet==1) { // Outer specific
    double Rcdphi,  Xr,  Yr,  Ur,  Vr;  getReducedOuter(X, Y, Z, module,Rcdphi, Xr, Yr, Ur, Vr);
    double Rcdphis, Xrs, Yrs, Urs, Vrs; getReducedOuter(Xs,Ys,Zs,module,Rcdphis,Xrs,Yrs,Urs,Vrs);
    if (doPrint)
      printf("Rr,Rrs: %.2f,%.2f, Ur,Urs: %.2f,%.2f, Vr,Vrs: %.2f,%.2f\n",
	     Rcdphi,Rcdphis,Ur,Urs,Vr,Vrs);
    double dRcdphi = 1000*(Rcdphi-Rcdphis);
    rs.Rr->Fill(dRcdphi);
    double dUr = 1000*(Ur-Urs), dVr = 1000*(Vr-Vrs);
    rs.Ur->Fill(dUr); rs.Vr->Fill(dVr);
  }
  else if (doPrint) printf("\n");
}
void getReducedCyMBaL(double X, double Y, unsigned int div,
		      double &Xr, double &Yr, double &Rr, double &phir,
		      int *zone, // Whether in peak, L/R edge, else
		      double *rot)   // Rotation angle
{
  // Transform to centre of curvature of cylindrical tiles
  // and rotate to phi = 0
  unsigned int iz = div>>3;
  int iphi = div%8, jphi = iphi%2; double phic = iphi*TMath::Pi()/4;
  if (rot) *rot = phic;
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
  if (zone) {
    // *****zone: Whether in peak, L/R edge, else
    // bin low up 
    // Inner: peak = 344 556.748 556.797
    //        low  = 313 555.234 555.283
    //        up   = 374 558.213 558.262
    // Outer: peak = 794 578.721 578.770
    //        low  = 763 577.207 577.256
    //        up   = 825 580.234 580.283
    const double cls[2] = {556.748,578.721}, cus[2] = {556.797,578.770};
    const double lls[2] = {555.234,577.207}, lus[2] = {555.283,577.256};
    const double uls[2] = {558.213,580.234}, uus[2] = {558.262,580.283};
    int io = (1<=iz && iz<=2)?0:1;
    if      (cls[io]<Rr && Rr<cus[io]) *zone = 0;
    else if (lls[io]<Rr && Rr<lus[io]) *zone = 1;
    else if (uls[io]<Rr && Rr<uus[io]) *zone = 2;
    else                               *zone = 3;
  }
}
void getReducedOuter(double X, double Y, double Z, unsigned int module,
		     double &Rcdphi, double &Xr, double &Yr,
		     double &Ur, double &Vr)
{
  // Rotate to phi = 0
  int iphi = module>>1; double phic = iphi*TMath::Pi()/6;
  //<constant name="MPGDOuterBarrelModule_zmin1"     value="164.5*cm"/>
  //<constant name="MPGDOuterBarrelModule_zmin2"     value="174.5*cm"/>
  double zmin1 = 1640.5, zmin2 = 1740.5;
  double modL = (zmin1+zmin2)/2, z0 = (zmin2-zmin1)/2;
  double dZ = (module&0x1)?z0+modL/2:z0-modL/2; Z -= dZ;
  TVector3 V(X,Y,Z);
  TRotation r;
  r.SetXAxis(TVector3(cos(-phic),sin(-phic),0));
  r.SetYAxis(TVector3(-sin(-phic),cos(-phic),0));
  V *= r;
  Rcdphi = V(0); Xr = V(1); Yr = V(2);
  Vr = (Yr+Xr)/sqrt(2); Ur = (Yr-Xr)/sqrt(2);
  
  /*
  printf(" 0x%02x %2d %6.1f,%6.1f,%4.1f %4.1f %6.1f,%6.1f,%4.1f\n",
	 module,iphi,X,Y,atan2(Y,X)/TMath::Pi()*12,phic/TMath::Pi()*12,
	 V(0),V(1),atan2(V(1),V(0)/TMath::Pi()*12));
  */
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
      printf("Evt #%d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
    int nHits = hits[idet]->size();
    printf("%2d/%d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx), %2d/%d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx)\n",
	   ih,nHits,mcIdx,module,quality,vID,hit.cellID>>32,
	   jh,nHits,mcJdx,modvle,qualjty,vJD,hjt.cellID>>32);
  }
  return mcJdx==mcIdx && modvle==module && qualjty==quality && vJD!=vID;
}
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
  const Vector3f &mom = hit.momentum, mqm = hjt.momentum;
  double Px = mom.x, Py = mom.y, Pz = mom.z, P = sqrt(Px*Px+Py*Py+Pz*Pz);
  double Qx = mqm.x, Qy = mqm.y, Qz = mqm.z, Q = sqrt(Qx*Qx+Qy*Qy+Qz*Qz);
  double u = dist/2/P, Ex = Mx+u*Px, Ey = My+u*Py, Ez = Mz+u*Pz;
  double v = dist/2/Q, Fx = Nx-v*Qx, Fy = Ny-v*Qy, Fz = Nz-v*Qz;
  double dext = sqrt((Ex-Fx)*(Ex-Fx)+(Ey-Fy)*(Ey-Fy)+(Ez-Fz)*(Ez-Fz));
  bool ok = dext<.025;
  // Debugging
  if ((verbose&0x100<<idet) || evtNum==evtToDebug) doDebug |= 0x100;
  if  (doDebug&0x100) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
  }
  bool doPrint = doDebug&0x100;
  if (!doPrint && !ok) {
    // Do print debugging message, even if not requested, when extrpolation
    // turns out to fail in case status is 0x3.
    unsigned status = getStatus(idet,ih); doPrint = status==0x3;
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
      int subVolID = hit.cellID>>28&0xf, subVolJD = hjt.cellID>>28&0xf;
      unsigned int continuation = 0;
      if (traversing) continuation |= 0x1;
      if (traversjng) continuation |= 0x2;
      for (int ij = 0; ij<2; ij++) {
	if (continuation&0x1<<ij) continue;
	// Try and rescue cases where traversing to one side.
	double Rr = ij==0 ? depth : dfpth;
	bool is0x3 = subVolID==0x3 || subVolJD==0x3;
	int module = hit.cellID>>12&0xfff;
	int sector = module/8; int inOut = (sector==1 || sector==2) ? 0 : 1;
	double midPlane = radii[0][inOut];
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
  }
  return ok;
}
bool recoEvents::coalesce(int idet, vector<int> coalesced, SimTrackerHitData &hext) {
  int ic, ih3, ih4, ih0; double eDep; for (ic = 0, ih3=ih4=ih0 = -1, eDep = 0;
					   ic<(int)coalesced.size(); ic++) {
    int ih = coalesced[ic]; SimTrackerHitData &hit = hits[idet]->at(ih);
    int subVolID = hit.cellID>>28&0xf;
    if      (subVolID==0x3) ih3 = ih;
    else if (subVolID==0x4) ih4 = ih;
    else if (subVolID==0x0) ih0 = ih;
    else if (ih0<0)         ih0 = ih;
    eDep += hit.eDep;
  }
  if (ih3>=0 && ih4>=0) {
    SimTrackerHitData &hit = hits[idet]->at(ih3), &hjt = hits[idet]->at(ih4);
    const Vector3d &pos = hit.position, &pqs = hjt.position;
    double Mx = pos.x, My = pos.y, Mz = pos.z;
    double Nx = pqs.x, Ny = pqs.y, Nz = pqs.z;
    Vector3d pext; pext.x = (Mx+Nx)/2; pext.y = (My+Ny)/2; pext.z = (Mz+Nz)/2;
    hext.position = pext; hext.eDep = eDep;
    hext.cellID = hit.cellID; hext.quality = hit.quality;
    return true;
  }
  else if (ih0>=0) {
    SimTrackerHitData &hit = hits[idet]->at(ih0);
    hext.position = hit.position; hext.eDep = eDep;
    hext.cellID = hit.cellID; hext.quality = hit.quality;
    return true;
  }
  else {
    printf("#%d: Inconsistency\n",evtNum);
    for (int ic = 0; ic<(int)coalesced.size(); ic++) {
      int ih = coalesced[ic]; SimTrackerHitData &hit = hits[idet]->at(ih);
      printf("%d: 0x%08lx,0x%08lx %.2f,%.2f,%.2f\n",
	     ih,hit.cellID&0xffffffff,hit.cellID>>32,
	     hit.position.x,hit.position.y,hit.position.z);
    }
    return false;
  }
  return false;
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
      printf("Evt #%d det %d\n",evtNum,idet); doDebug |= 0x20000;
    }
  }
  bool doExtend = true; int direction = +1;
  int subVolID = hit.cellID>>28&0xf; if (subVolID==0x3 || subVolID==0x4) {
    // Determine pathDepth, i.e. path projected along cylinder's radius
    if (idet==0) {
      double pathDepth, Rr; doExtend = checkTraversing(idet,hit,pathDepth,Rr);
      if ((doDebug&0x100) && doExtend)
	printf(" 0x%08lx,0x%08lx Traversing? %d: pathDepth %.3f Rr %.3f\n",
	       cID&0xffffffff,cID>>32,doExtend,pathDepth,Rr);
      if (!doExtend) {
	// Try and rescue cases where traversing to one side.
	int module = hit.cellID>>12&0xfff;
	int sector = module/8; int inOut = (sector==1 || sector==2) ? 0 : 1;
	double midPlane = radii[0][inOut];
	bool isInner = subVolID==0x3 ? Rr > midPlane : Rr < midPlane;
	doExtend = isInner;
	if (doDebug&0x100)
	  printf(" 0x%08lx,0x%08lx Traversing? %d: pathDepth %.3f Rr %.3f/%.3f\n",
		 cID&0xffffffff,cID>>32,doExtend,pathDepth,Rr,midPlane);
      }
      direction = pathDepth>0 ? +1 : -1;
    }
    if (doExtend) {
      if (subVolID==3) {
	double u =  direction*hit.pathLength/2/P;
	Ex = Mx+u*Px; Ey = My+u*Py; Ez = Mz+u*Pz;
      }
      else if (subVolID==4) {
	double u = -direction*hit.pathLength/2/P;
	Ex = Mx+u*Px; Ey = My+u*Py; Ez = Mz+u*Pz;
      }
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
      double Xr, Yr, Rr, phir;
      getReducedCyMBaL(Mx,My,div,Xr,Yr,Rr,phir); double MR = Rr;
      getReducedCyMBaL(Ex,Ey,div,Xr,Yr,Rr,phir); double ER = Rr;
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
    double Xr, Yr, Rr, phir, rot;
    getReducedCyMBaL(Mx,My,div,Xr,Yr,Rr,phir,0,&rot);
    double crot = std::cos(rot), srot = std::sin(rot);
    double Pxr = crot*Px+srot*Py, Pyr = -srot*Px+crot*Py;
    double cTheta = (Pxr*Xr+Pyr*Yr)/P/sqrt(Xr*Xr+Yr*Yr);
    pathDepth = hit.pathLength*cTheta; depth = Rr;
    double thickness = thicknesses[idet];
    traversing = fabs(fabs(pathDepth)-thickness)<.200; // Somewhat large 200 µm tolerance
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
      printf("Evt #%d det %d\n",evtNum,idet); doDebug |= 0x20000;
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
    double Xr, Yr, Rr, phir; getReducedCyMBaL(X,Y,div,Xr,Yr,Rr,phir);
    printf("%10s 0x%08lx,0x%08x X,Y,Z %7.2f,%7.2f,%8.2f Rr %7.3f mm phir %6.3fπ\n",
	   "",cellID&0xffffffff,cell,X,Y,Z,Rr,phir/pi);
  }
  else if (idet==1) {
    double Rcdphi, Xr, Yr, Ur, Vr; getReducedOuter(X,Y,Z,module,Rcdphi,Xr,Yr,Ur,Vr);
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
    printf("#%d %d/%d:%d(0x%08lx,0x%08lx) %d/%d 0x%08lx,0x%08lx\n",evtNum,
	   is,ncoas,cIndex,hit.cellID&0xffffffff,hit.cellID>>32,
	   ir,nRecs,rec.cellID&0xffffffff,rec.cellID>>32);
  }
}
void recoEvents::debugAssoc(int idet)
{
  if (verbose&0x10<<idet) doDebug |= 0x10;
  if (doDebug&0x10) {
    if (!(doDebug&0x20000)) {
      printf("Evt #%d det %d\n",evtNum,idet); doDebug |= 0x20000;
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
