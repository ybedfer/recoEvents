## recoEvents: Project the `events` TTree of the `podio_output` of `eicrecon`

### Installation
```
mkdir build; mkdir install; cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_BUILD_TYPE=Debug ..
make
make install
```

### Usage, from the ROOT command line:
```
.L ../recoEvents/install/librecoEvents.so
TFile *_file0 = TFile::Open("podio_output.root");
TTree *events = (TTree*)gDirectory->Get("events");
// Instantiate a recoEvents object w/ as arg.s:
// - pointer to events TTree
// - bit pattern: 0x1:CyMBaL, 0x2:Outer, 0x4:Vertex, 0x8:Si.
recoEvents ana(events,0xf);

// Loop on events in tree
ana.Loop();

// Draw Histos
ana.DrawResiduals(0x1);   // TCanvas of various residuals for 0x1:CyMBaL
ana.DrawphithZR(1,0x1);   // TCanvas of Rec hits of 0x1:CyMBaL...
ana.DrawphithZR(1,0x1,1); // ...same w/ some decorations
ana.DrawphithZR(1,0x2,1); // ...same for Sim hits of 0x2:outer
new TCanvas("c2D");
ana.recHs[0].XY->Draw();  // Draw Y vs. X for Rec hits for [0]:CyMBal
ana.simHs[2].ZR->Draw();  // Draw R vs. Z for Sim hits for [2]:Vertex

// Instantiate a second recoEvents object
TFile *_file0 = TFile::Open("podio_output.2.root");
TTree *event2 = (TTree*)gDirectory->Get("events");
recoEvents ana2(event2,0x1);

// Histos can also be accessed and listed from the ROOT file system.
```


