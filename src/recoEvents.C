#define recoEvents_cxx
#include "recoEvents.h"
#include <TH2.h>
#include <TVector3.h>
#include <TRotation.h>

// Globals
unsigned int recoEvents::nObjCreated = 0;

// Interfaces
void getReducedCyMBaL(double X, double Y, unsigned int div,
		      double &Xr, double &Yr, double &Rr, double &phir,
		      int *zone = 0); // Whether in peak, L/R edge, else
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
    if (verbose&0xfff) printf("Event #%d\n",evtNum);

    //int treenumber = fChain->GetTreeNumber();
    // if (Cut(ientry) < 0) continue;

    // ***** #MCs
    int nMCs = mcParticles->size(); hMult->Fill(nMCs);
    if      (requireNMCs>0)  { if (nMCs!=requireNMCs) continue; }
    else if (requireNMCs<0)  { if (nMCs< requireNMCs) continue; }

    // Init requirement
    unsigned int modID = 0, prvID = 0;
    int nGoodCoaHs = 0, nGoodCHs = 0; // #coalesced-hits fulfilling requirements
    int nHsPerModule = 0, nHsPM = 0;  // #coalesced-hits per module
    // Modules where all hits fulfill quality requirement
    unsigned int modulesOK = 0, moduleOK = 0;
    for (int idet = 0; idet<N_DETs; idet++) {
      if (!(0x1<<idet&processedDetectors)) continue;
      // *************** LOOP ON SELECTED DETECTORS
      int nHits = hits[idet]->size(); // Sim hits
      vector<SimTrackerHitData> coalescedHs; map<int,int> sim2coa;
      for (int ih = 0; ih<nHits; ih++) {
	// ********** LOOP ON sim HITS
	SimTrackerHitData &hit = hits[idet]->at(ih);
	int jh = ih; vector<int> coalesced; int nCoaHs = coalescedHs.size();
	while (++jh<nHits && samePMO(idet,ih,jh)) {
	  if (extrapolate(idet,ih,jh)) {
	    if (coalesced.empty()) {
	      coalesced.push_back(ih); sim2coa[ih] = nCoaHs;
	    }
	    coalesced.push_back(jh); sim2coa[jh] = nCoaHs;
	  }
	  else break;
	}
	if (coalesced.size()<2) {
	  if (hit.quality==0) { // Extend to edge of pathLength
	    SimTrackerHitData hext; extend(idet,ih,hext);
	    coalescedHs.push_back(hext);
	  }
	  else // Do not extend secondary: its position is, in general, well
	    // represented by its SimTrackerHit position.
	    coalescedHs.push_back(hit);
	  sim2coa[ih] = nCoaHs;
	}
	else {
	  SimTrackerHitData hext; if (!coalesce(idet,coalesced,hext)) return;
	  sim2coa[ih] = nCoaHs; coalescedHs.push_back(hext);
	  ih += coalesced.size()-1;
	}
	// LONE PRIMARY?
	unsigned int status = getStatus(idet,ih,hit.quality);
	modID = hit.cellID&0xfffffff;
	if (modID==prvID)
	  nHsPM++;
	else {
	  nHsPM = 1; modulesOK |= moduleOK;
	}
	if (nHsPM>nHsPerModule) nHsPerModule = nHsPM;
	if (status==0x3) {
	  if (modID==prvID)
	    nGoodCHs++;
	  else {
	    nGoodCHs = 1; moduleOK = modID;
	  }
	  if (nGoodCHs>nGoodCoaHs) nGoodCoaHs = nGoodCHs;
	}
	else if (!(status&0x2))
	  moduleOK = 0;
	prvID = modID;
      }
      if (moduleOK) modulesOK |= moduleOK;

      bool doDebug = nMCs==1 && (nGoodCoaHs>2 || evtNum==94 || evtNum==2402) && idet==0;
      if (doDebug) printf("Evt %d\n",evtNum);
      int nCoaHs = coalescedHs.size();
      // ***** EVENT SELECTION
      if      (requireNHits>0) { if (nHsPerModule!=requireNHits) continue; }
      else if (requireNHits<0) { if (nHsPerModule< requireNHits) continue; }
      for (int ih = 0; ih<nCoaHs; ih++) {
	SimTrackerHitData &hit = coalescedHs[ih];
	// ***** HIT SELECTION
	unsigned int status = getStatus(idet,ih,sim2coa);
	if (verbose&0x1<<idet) debugHit(idet,ih,hit,status);
	if (status!=0x3) continue;
	if (requireQuality>1) {
	  modID = hit.cellID&0xfffffff; if ((modID&modulesOK)!=modID) continue;
	}
	// Fill hit
	const Vector3d &pos = hit.position;
	double X = pos.x, Y = pos.y, Z = pos.z;
	fillHit(0,idet,X,Y,Z,hit.cellID);       // ***** FILL sim HISTOS
	if (verbose&0x1000) printHit(idet,X,Y,Z,hit.cellID);
	if (doDebug) printHit(idet,X,Y,Z,hit.cellID);
      }
      if (verbose&0x1000) printf("\n");
      if (doDebug) printf("\n");
      if (!reconstruction) continue;
      int nArhs = arhs[idet]->size(); // Associations raw hit
      int nAshs = ashs[idet]->size(); // Associations sim hit
      int nARhs = aRhs[idet]->size(); // Associations Rec Hit -> raw hit
      int nRecs = recs[idet]->size(); // Rec hits
      if (!(nArhs||nAshs||nHits||nRecs)) continue;
      map<int,vector<int>,less<int>> rec2sims;
      if (verbose&0x10<<idet) debugAssoc(idet);
      // ********** BUILD rec <-> sim MAP
      if (nAshs!=nArhs) {
	printf("Warning: %3d det %d: ash(%d)!=arh(%d)\n",
	       (int)jentry,idet,nAshs,nArhs);
      }
      else {
	// raw -> Rec
	// raw can only be associated to one Rec, by construction...
	// ...and vice-versa, by convention imprinted in edm4hep::TrackerHit
	map<int,int,less<int>> raw2rec;
	for (int iR = 0; iR<nARhs; iR++) {
	  podio::ObjectID &aRh = aRhs[idet]->at(iR); raw2rec[aRh.index] = iR;
	}
	// raw <-> sims
	for (int ih = 0; ih<nArhs; ih++) {
	  podio::ObjectID &arh = arhs[idet]->at(ih); //RawHitAssociations_rawHit
	  int rIndex = arh.index;
	  map<int,int>::const_iterator ir = raw2rec.find(rIndex);
	  if (ir==raw2rec.end())
	    // One raw may be lost here, but that lost one is then associated to
	    // the same sims as the retained one.
	    continue;
	  int RIndex = ir->second;
	  podio::ObjectID &ash = ashs[idet]->at(ih); //RawHitAssociations_simHit
	  int sIndex = ash.index;
	  map<int,int>::const_iterator is = sim2coa.find(sIndex);
	  if (is==sim2coa.end()) continue;
	  int cIndex = is->second;
	  map<int,vector<int>>::iterator im = rec2sims.find(RIndex);
	  if (im==rec2sims.end()) {
	    vector<int> sims; sims.push_back(cIndex); rec2sims[RIndex] = sims;
	  } else {
	    // Not twice same SimHit associated to RecHit.
	    vector<int> &sims = im->second;
	    int is, match; for (int is=match = 0; is<(int)sims.size(); is++) {
	      if (sims[ih]==cIndex) { match = 1; break; }
	    }
	    if (!match) sims.push_back(cIndex);
	  }
	}
	if (verbose&0x10<<idet) debugAssoc(raw2rec,sim2coa,rec2sims);
      }
      for (int ih = 0; ih<nRecs; ih++) {
	// ********** LOOP ON rec HITS
	edm4eic::TrackerHitData &rec = recs[idet]->at(ih);
	const Vector3f &pos = rec.position;
	double X = pos.x, Y = pos.y, Z = pos.z;    
	if (verbose&0x1<<idet) debugRec(idet,ih,rec);
	if (verbose&0x1000) printHit(idet,X,Y,Z,rec.cellID);
	// ********** RESIDUALS
	map<int,vector<int>,less<int>>::const_iterator im = rec2sims.find(ih);
	if (im==rec2sims.end()) {
	  printf("Warning: Entry %3d det %d: rec %d not associated\n",
		 (int)jentry,idet,ih);
	}
	else {
	  const vector<int> &sims = im->second;
	  int is, selecRec; for (is=selecRec = 0; is<(int)sims.size(); is++) {
	    int sIndex = sims[is]; // ***** REFERENCE TO sim HIT
	    if (evtNum==19 || evtNum==29 || evtNum==2402) {
	      SimTrackerHitData &hit = coalescedHs[sIndex];
	      printf("#%d %d/%d:%d(0x%08lx,0x%08lx) %d/%d 0x%08lx,0x%08lx\n",evtNum,
		     is,(int)sims.size(),sIndex,hit.cellID&0xffffffff,hit.cellID>>32,
		     ih,nRecs,rec.cellID&0xffffffff,rec.cellID>>32);
	    }
	    if (sIndex<0 || nHits<=sIndex) {
	      printf("Warning: %3d det %d: raw %d <-> sim %d\n",
		     (int)jentry,idet,ih,sIndex);
	    }
	    else {
	      SimTrackerHitData &hit = coalescedHs[sIndex];
	      // ***** HIT SELECTION
	      unsigned int status = getStatus(idet,sIndex,sim2coa);
	      if (status!=0x3) continue;
	      if (requireQuality>1) {
		modID = rec.cellID&0xfffffff; if ((modID&modulesOK)!=modID) continue;
	      }
	      selecRec = 1;
	      const Vector3d &psim = hit.position;
	      fillResids(idet,pos,psim,rec.cellID); // ***** FILL RESIDUAL
	    }
	  }
	  if (selecRec)
	    fillHit(1,idet,X,Y,Z,rec.cellID);     // ***** FILL rec HISTOS
	}
      }
      if (verbose&0x1000) printf("\n");
    }
  }
}
unsigned int recoEvents::getStatus(int idet, int ih, int quality)
{
  // Returned status is OR of conformity to quality and PDG requirements.
  // 0x1: PDG OK.
  // 0x2: quality OK.
  unsigned status = 0;
  podio::ObjectID &amc = amcs[idet]->at(ih); int mcIdx = amc.index;
  MCParticleData &part = mcParticles->at(mcIdx);
  if (!requirePDG || part.PDG==requirePDG) status |= 0x1;
  if (!requireQuality || quality==0)       status |= 0x2;
  return status;
}
unsigned int recoEvents::getStatus(int idet, int ih, map<int,int> &sim2coa)
{
  unsigned int status = 0;
  for (auto is = sim2coa.cbegin(); is!=sim2coa.cend(); is++) {
    if (is->second==ih) {
      int jh = is->first; SimTrackerHitData &sim = hits[idet]->at(jh);
      status = getStatus(idet,jh,sim.quality);
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
  // ***** MODULE SELECTION
  if (requireModule>=0 && module!=requireModule) return;
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
  double X =  pos.x,  Y =  pos.y,  Z =  pos.z;
  double Xs = psim.x, Ys = psim.y, Zs = psim.z;
  double dX = 1000*(X-Xs), dY = 1000*(Y-Ys), dZ = 1000*(Z-Zs); // Residuals
  double R2 = X*X+Y*Y,      R = sqrt(R2);
  double Rs2 = Xs*Xs+Ys*Ys, Rs = sqrt(Rs2);
  double dR = 1000*(R-Rs);
  double phi = atan2(Y,X), phis = atan2(Ys,Xs);
  double dphi = 1000*(phi-phis);
  double D = sqrt(R2+Z*Z), Ds = sqrt(Rs2+Zs*Zs), dD = 1000*(D-Ds);
  if (verbose&0x2000) printf(" dX,dY,dZ: %.2f,%.2f,%.2f",dX,dY,dZ);
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  // ***** STRIP: Convert stripID -> strip. Is it valid? 
  if (!parseStrip(idet,1,strip)) return;
  Resids &rs = resHs[strip][idet];
  rs.X->Fill(dX); rs.Y->Fill(dY); rs.Z->Fill(dZ); rs.R->Fill(dR); rs.D->Fill(dD);
  rs.phi->Fill(dphi);
  if      (idet==0) { // CyMBaL specific
    double Xd, Yd, Rr, phir;   getReducedCyMBaL(X,Y,div,Xd,Yd,Rr,phir);
    double Rrs, phirs;         getReducedCyMBaL(Xs,Ys,div,Xd,Yd,Rrs,phirs);
    if (verbose&0x2000) printf("Rr,Rrs: %.2f,%.2f\n",Rr,Rrs);
    double dRr = 1000*(Rr-Rrs), dphir = 1000*(phir-phirs);
    rs.Rr->Fill(dRr); rs.phir->Fill(dphir);
    if (strip==1 && dZ<-600) {
      printf("#%5d 0x%08lx,0x%08lx  %7.2f,%7.2f,%8.2f %7.3f mm\n",
	     evtNum,cellID&0xffffffff,cellID>>32,X,Y,Z,Rr);
      printf("#   SimHit                    %7.2f,%7.2f,%8.2f %7.3f mm\n",
	     Xs,Ys,Zs,Rrs);
    }
  }
  else if (idet==1) { // Outer specific
    double Rcdphi,  Xr,  Yr,  Ur,  Vr;  getReducedOuter(X, Y, Z, module,Rcdphi, Xr, Yr, Ur, Vr);
    double Rcdphis, Xrs, Yrs, Urs, Vrs; getReducedOuter(Xs,Ys,Zs,module,Rcdphis,Xrs,Yrs,Urs,Vrs);
    if (verbose&0x2000)
      printf("Rr,Rrs: %.2f,%.2f, Ur,Urs: %.2f,%.2f, Vr,Vrs: %.2f,%.2f\n",
	     Rcdphi,Rcdphis,Ur,Urs,Vr,Vrs);
    double dRcdphi = 1000*(Rcdphi-Rcdphis);
    rs.Rr->Fill(dRcdphi);
    double dUr = 1000*(Ur-Urs), dVr = 1000*(Vr-Vrs);
    rs.Ur->Fill(dUr); rs.Vr->Fill(dVr);
  }
  else if (verbose&0x2000) printf("\n");
}
void getReducedCyMBaL(double X, double Y, unsigned int div,
		      double &Xr, double &Yr, double &Rr, double &phir,
		      int *zone) // Whether in peak, L/R edge, else
{
  // Transform to centre of curvature of cylindrical tiles
  // and rotate to phi = 0
  unsigned int iz = div>>3;
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
void recoEvents::printHit(int idet,
			  double X, double Y, double Z, unsigned long cellID)
{
  double R2 = X*X+Y*Y, R = sqrt(R2); double phi = atan2(Y,X);
  const double pi = TMath::Pi();
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  int cell = cellID>>32;
  if (verbose&0x4000) {
    printf(" 0x%08lx,0x%08x  %7.2f,%7.2f,%8.2f %7.2f cm %6.3fπ",
	   cellID&0xfffffff,cell,X,Y,Z,R,phi/pi);
    if      (idet==1) printf(" %5.1f 0x%02x 0x%x",phi/pi*12,module,module>>1);
    else if (idet==0) printf(" %5.1f 0x%02x 0x%x",phi/pi*8, module,module>>1);
  }
  else if (idet==2) {
    unsigned int layer = module&0xf, tile = module>>4;
    printf(" %d,%2d,0x%08x  %7.2f,%7.2f,%8.2f %7.2f cm %6.3fπ",
	   layer,tile,cell,X,Y,Z,R,phi/pi);
  }
  else if (idet==0) {
    double Xr, Yr, Rr, phir; int zone; getReducedCyMBaL(X,Y,div,Xr,Yr,Rr,phir,&zone);
    //printf(" %2d,0x%08x  %7.2f,%7.2f,%8.2f cm %6.3fπ  %7.2f,%7.2f %7.2f cm %6.3fπ",
    //	   module,cell,X,Y,Z,phi/pi,Xr,Yr,Rr,phir/pi);
    printf(" 0x%8lx,0x%08x  %7.2f,%7.2f,%8.2f %7.2f cm\n",
	   cellID&0xffffffff,cell,X,Y,Z,Rr);
  }
  else if (idet==1) {
    double Rcdphi, Xr, Yr, Ur, Vr; int zone; getReducedOuter(X,Y,Z,module,Rcdphi,Xr,Yr,Ur,Vr);
    printf(" %2d,0x%08x  %7.2f,%7.2f,%8.2f %7.2f,%7.2f %7.2f,%7.2f",
	   module,cell,X,Y,Z,Xr,Yr,Ur,Vr);
  }
  else {
    printf("** recoEvents::getDetHit: Invalid det# = %d\n",idet);
    exit(1);
  }
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
  if (verbose&0x100<<idet)
    printf("%d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx), %d: PMO=%2d,%2d,%d (0x%08lx,0x%08lx)\n",
	   ih,mcIdx,module,quality,vID,hit.cellID>>32,
	   jh,mcJdx,modvle,qualjty,vJD,hjt.cellID>>32);
  return mcJdx==mcIdx && modvle==module && qualjty==quality && vJD!=vID;
}
bool recoEvents::extrapolate(int idet, int ih, int jh)
{
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
  bool doPrint = verbose&0x100<<idet;
  if (!doPrint && !ok) {
    unsigned status = getStatus(idet,ih,hit.quality); doPrint = status==0x3;
  }
  if (doPrint) {
    printf("#%4d %d,%d dext %f\n",evtNum,ih,jh,dext);
    printf("%2d: 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f\n",
	   ih,hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,Ex,Ey,Ez);
    printf("%2d: 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f\n",
	   jh,hjt.cellID&0xffffffff,hjt.cellID>>32,Nx,Ny,Nz,Fx,Fy,Fz);
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
    eDep += hit.EDep;
  }
  if (ih3>=0 && ih4>=0) {
    SimTrackerHitData &hit = hits[idet]->at(ih3), &hjt = hits[idet]->at(ih4);
    const Vector3d &pos = hit.position, &pqs = hjt.position;
    double Mx = pos.x, My = pos.y, Mz = pos.z;
    double Nx = pqs.x, Ny = pqs.y, Nz = pqs.z;
    Vector3d pext; pext.x = (Mx+Nx)/2; pext.y = (My+Ny)/2; pext.z = (Mz+Nz)/2;
    hext.position = pext; hext.EDep = eDep;
    hext.cellID = hit.cellID; hext.quality = hit.quality;
    return true;
  }
  else if (ih0>=0) {
    SimTrackerHitData &hit = hits[idet]->at(ih0);
    hext.position = hit.position; hext.EDep = eDep;
    hext.cellID = hit.cellID; hext.quality = hit.quality;
    return true;
  }
  else {
    printf("#%d: Inconsistency\n",evtNum);
    for (int ic = 0; ic<(int)coalesced.size(); ic++) {
      int ih = coalesced[ic]; SimTrackerHitData &hit = hits[idet]->at(ih);
      printf("%d: 0x%08lx 0x%lx %.2f,%.2f,%.2f\n",
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
  // - Expected to be called for hit <ih> being...
  //  ...primary,
  //  ...w/o any counterpart.
  //   This corresponds to a particle depositing energy in only one SUBVOLUME.
  // - Extend it, i.e. extrapolate it to the end of its pathLength, to mimic
  //  MPGDTrackerDigi (which in turn tries and emulates what we would get w/ a
  //  single sensitive volume).
  //   It's an approximation: extrapolation should go up to REFERENCE SUBVOLUME.
  //  But it's the best we can do, w/o knowledge of the detector geometry.
  SimTrackerHitData &hit = hits[idet]->at(ih); unsigned long cID = hit.cellID;
  const Vector3d &pos = hit.position; double Mx = pos.x, My = pos.y, Mz = pos.z;
  const Vector3f &mom = hit.momentum; double Px = mom.x, Py = mom.y, Pz = mom.z;
  double P = sqrt(Px*Px+Py*Py+Pz*Pz), Ex, Ey, Ez;
  int subVolID = hit.cellID>>28&0xf; if (subVolID==3) {
    double u =  hit.pathLength/2/P;  Ex = Mx+u*Px; Ey = My+u*Py; Ez = Mz+u*Pz;
  }
  else if (subVolID==4) {
    double u = -hit.pathLength/2/P; Ex = Mx+u*Px; Ey = My+u*Py; Ez = Mz+u*Pz;
  }
  else {
    Ex = Mx; Ey = My; Ez = Mz;
  }
  Vector3d pext; pext.x = Ex; pext.y = Ey; pext.z = Ez;
  hext.position = pext; hext.EDep = hit.EDep;
  hext.cellID = hit.cellID; hext.quality = hit.quality;
  if (verbose&0x100<<idet) {
    printf("#%5d: %d extend\n",evtNum,ih);
    if (idet==0) {
      unsigned int module, div, strip; parseCellID(idet,hit.cellID,module,div,strip);
      double Xr, Yr, Rr, phir;
      getReducedCyMBaL(Mx,My,div,Xr,Yr,Rr,phir); double Mr = Rr;
      getReducedCyMBaL(Ex,Ey,div,Xr,Yr,Rr,phir); double Er = Rr;
      printf(" 0x%08lx 0x%08lx %.2f,%.2f,%.2f (%.3f) -> %.2f,%.2f,%.2f (%.3f) mm\n",
	     hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,Mr,Ex,Ey,Ez,Er);
    }
    else
      printf(" 0x%08lx 0x%08lx %.2f,%.2f,%.2f -> %.2f,%.2f,%.2f mm\n",
	     hit.cellID&0xffffffff,hit.cellID>>32,Mx,My,Mz,Ex,Ey,Ez);
      
  }
}
void recoEvents::debugHit(int idet, int ih, SimTrackerHitData &hit, unsigned int status)
{
  int nHits = hits[idet]->size();
  const Vector3d &pos = hit.position;
  printf("hit %d/%d: %6.1f,%6.1f,%6.1f 0x%lx 0x%x\n",
	 ih,nHits,pos.x,pos.y,pos.z,hit.cellID>>32,status);
}
void recoEvents::debugRec(int idet, int ih, edm4eic::TrackerHitData &rec)
{
  int nRecs = recs[idet]->size();
  const Vector3f &pos = rec.position;
  printf("rec %d/%d: %6.1f,%6.1f,%6.1f 0x%lx\n",
	 ih,nRecs,pos.x,pos.y,pos.z,rec.cellID>>32);
}
void recoEvents::debugAssoc(int idet)
{
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
void recoEvents::debugAssoc(map<int,int>raw2rec, map<int,int> sim2coa, map<int,vector<int>> rec2sims)
{				  
  map<int,int>::const_iterator ir;
  for (ir = raw2rec.cbegin(); ir != raw2rec.cend(); ir++) {
    printf("raw2rec %d -> %d\n",ir->first,ir->second);
  }
  for (ir = sim2coa.cbegin(); ir != sim2coa.cend(); ir++) {
    printf("sim2coa %d -> %d\n",ir->first,ir->second);
  }
  map<int,vector<int>>::const_iterator im;
  for (im = rec2sims.cbegin(); im != rec2sims.cend(); im++) {
    printf("rec2sims: %d ->",im->first);
    const vector<int> &sims = im->second;
    for (int is = 0; is<(int)sims.size(); is++) printf(" %d",sims[is]);
    printf("\n");
  }
}
