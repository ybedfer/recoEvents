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

void recoEvents::Loop(int nEvents)
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
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
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

    //int treenumber = fChain->GetTreeNumber();
    // if (Cut(ientry) < 0) continue;
      
    hMult->Fill(mcParticles->size());

    for (int idet = 0; idet<N_DETs; idet++) {
      if (!(0x1<<idet&processedDetectors)) continue;
      // *************** LOOP ON SELECTED DETECTORS
      int nHits = hits[idet]->size(); // Sim hits
      // ***** SELECTION
      if      (requireNHits>0) { if (nHits!=requireNHits) continue; }
      else if (requireNHits<0) { if (nHits< requireNHits) continue; }
      for (int ih = 0; ih<nHits; ih++) {
	// ********** LOOP ON sim HITS
	SimTrackerHitData &hit = hits[idet]->at(ih);
	const Vector3d &pos = hit.position;
	double X = pos.x, Y = pos.y, Z = pos.z;
	// ***** STATUS
	unsigned int status = getStatus(idet,ih,hit.quality); if (status!=0x3) {
	  if (verbose&0x1) printf("hit %d: %6.1f,%6.1f,%6.1f 0x%lx 0x%x\n",
				  ih,X,Y,Z,hit.cellID>>32,status);
	  continue;
	}
	fillHit(0,idet,X,Y,Z,hit.cellID);       // ***** FILL sim HISTOS
	if (verbose&0x1) printf("hit %d: %6.1f,%6.1f,%6.1f 0x%lx\n",
				ih,X,Y,Z,hit.cellID>>32);
	if (verbose&0x6) printHit(idet,X,Y,Z,hit.cellID);
      }
      if (verbose&0x6) printf("\n");
      if (!reconstruction) continue;
      int nArhs = arhs[idet]->size(); // Associations raw hit
      int nAshs = ashs[idet]->size(); // Associations sim hit
      int nRecs = recs[idet]->size(); // Rec hits
      if (!(nArhs||nAshs||nHits||nRecs)) continue;
      if (verbose&0x1) printf("%3d det %d: arh,ash,hit,rec %d,%d,%d,%d\n",
			      (int)jentry,idet,nArhs,nAshs,nHits,nRecs);
      map<int,int,less<int>> raw2sim;
      if (verbose&0x1) {
	for (int ih = 0; ih<nArhs; ih++) {
	  podio::ObjectID &arh = arhs[idet]->at(ih);
	  int rIndex = arh.index; unsigned int cID = arh.collectionID;
	  printf("arh %d: %d 0x%x\n",ih,rIndex,cID);
	}
	for (int ih = 0; ih<nAshs; ih++) {
	  podio::ObjectID &ash = ashs[idet]->at(ih);
	  int sIndex = ash.index; unsigned int cID = ash.collectionID;
	  printf("ash %d: %d 0x%x\n",ih,sIndex,cID);
	}
      }
      // ********** BUILD raw <-> sim MAP
      if (nAshs!=nArhs) {
	printf("Warning: %3d det %d: ash(%d)!=arh(%d)\n",
	       (int)jentry,idet,nAshs,nArhs);
      }
      else {
	for (int ih = 0; ih<nArhs; ih++) {
	  podio::ObjectID &arh = arhs[idet]->at(ih);
	  int rIndex = arh.index;
	  podio::ObjectID &ash = ashs[idet]->at(ih);
	  int sIndex = ash.index;
	  raw2sim[rIndex] = sIndex;
	}
      }
      for (int ih = 0; ih<nRecs; ih++) {
	// ********** LOOP ON rec HITS
	edm4eic::TrackerHitData &rec = recs[idet]->at(ih);
	const Vector3f &pos = rec.position;
	double X = pos.x, Y = pos.y, Z = pos.z;    
	fillHit(1,idet,X,Y,Z,rec.cellID);       // ***** FILL rec HISTOS
	if (verbose&0x1) printf("rec %d: %6.1f,%6.1f,%6.1f 0x%lx\n",
				ih,X,Y,Z,rec.cellID>>32);
	if (verbose&0x6) printHit(idet,X,Y,Z,rec.cellID);
	// ********** RESIDUALS
	map<int,int,less<int>>::const_iterator im = raw2sim.find(ih);
	if (im==raw2sim.end()) {
	  printf("Warning: %3d det %d: raw %d not associated\n",
		 (int)jentry,idet,ih);
	}
	else {
	  int sIndex = im->second; // ***** REFERENCE TO sim HIT
	  if (sIndex<0 || nHits<=sIndex) {
	    printf("Warning: %3d det %d: raw %d <-> sim %d\n",
		   (int)jentry,idet,ih,sIndex);
	  }
	  else {
	    SimTrackerHitData &hit = hits[idet]->at(sIndex);
	    if (!requireQuality || hit.quality==0) {
	      const Vector3d &psim = hit.position;
	      fillResids(idet,pos,psim,rec.cellID); // ***** FILL RESIDUAL
	    }
	  }
	}
      }
      if (verbose&0x6) printf("\n");
    }
  }
}
unsigned int recoEvents::getStatus(int idet, int ih, int quality)
{
  // Returned status is OR of conformity to quality and PDG requirements.
  // 0x3 means OK.
  unsigned status = 0;
  podio::ObjectID &amc = amcs[idet]->at(ih); int mcIdx = amc.index;
  MCParticleData &part = mcParticles->at(mcIdx);
  if (!requirePDG || part.PDG==requirePDG) status |= 0x1;
  if (!requireQuality || quality==0)       status |= 0x2;
  return status;
}
void recoEvents::fillHit(int simOrRec, int idet,
			 double X, double Y, double Z, unsigned long cellID)
{
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  // ***** STRIP: Convert stripID -> strip. Is it valid? 
  if (!parseStrip(idet,simOrRec,strip)) return;
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
  if (verbose&0x2) printf(" dX,dY,dZ: %.2f,%.2f,%.2f",dX,dY,dZ);
  unsigned int module, div, strip; parseCellID(idet,cellID,module,div,strip);
  // ***** STRIP: Convert stripID -> strip. Is it valid? 
  if (!parseStrip(idet,1,strip)) return;
  Resids &rs = resHs[strip][idet];
  rs.X->Fill(dX); rs.Y->Fill(dY); rs.Z->Fill(dZ); rs.R->Fill(dR); rs.D->Fill(dD);
  rs.phi->Fill(dphi);
  if      (idet==0) { // CyMBaL specific
    double Xd, Yd, Rr, phir;   getReducedCyMBaL(X,Y,div,Xd,Yd,Rr,phir);
    double Rrs, phirs;         getReducedCyMBaL(Xs,Ys,div,Xd,Yd,Rrs,phirs);
    if (verbose&0x2) printf("Rr,Rrs: %.2f,%.2f\n",Rr,Rrs);
    double dRr = 1000*(Rr-Rrs), dphir = 1000*(phir-phirs);
    rs.Rr->Fill(dRr); rs.phir->Fill(dphir);
  }
  else if (idet==1) { // Outer specific
    double Rcdphi,  Xr,  Yr,  Ur,  Vr;  getReducedOuter(X, Y, Z, module,Rcdphi, Xr, Yr, Ur, Vr);
    double Rcdphis, Xrs, Yrs, Urs, Vrs; getReducedOuter(Xs,Ys,Zs,module,Rcdphis,Xrs,Yrs,Urs,Vrs);
    if (verbose&0x2)
      printf("Rr,Rrs: %.2f,%.2f, Ur,Urs: %.2f,%.2f, Vr,Vrs: %.2f,%.2f\n",
	     Rcdphi,Rcdphis,Ur,Urs,Vr,Vrs);
    double dRcdphi = 1000*(Rcdphi-Rcdphis);
    rs.Rr->Fill(dRcdphi);
    double dUr = 1000*(Ur-Urs), dVr = 1000*(Vr-Vrs);
    rs.Ur->Fill(dUr); rs.Vr->Fill(dVr);
  }
  else if (verbose&0x2) printf("\n");
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
  // system: CyMBal = 61, Outer = 64, ... (It's not checked that this is what we get in "ID")
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
  // ...whether it's simHit or recHit
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
  if (verbose&0x4) {
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
    printf(" %2d,0x%08x  %7.2f,%7.2f,%8.2f cm %6.3fπ  %7.2f,%7.2f %7.2f cm %6.3fπ",
	   module,cell,X,Y,Z,phi/pi,Xr,Yr,Rr,phir/pi);
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
