from DDSim.DD4hepSimulation import DD4hepSimulation 
from g4units import mm, GeV, MeV, keV, deg
SIM = DD4hepSimulation()
SIM.gun.thetaMin = 2*deg
SIM.gun.thetaMax = 178*deg
SIM.gun.distribution="uniform"
SIM.gun.momentumMin = 10*GeV
SIM.gun.momentumMax = 10*GeV
SIM.gun.particle = "mu-"
SIM.gun.multiplicity = 1
SIM.random.enableEventSeed = True
SIM.outputFile = "epic_craterlake_tracking_only.edm4hep.root"
SIM.filter.filters['edep0'] = dict(name="EnergyDepositMinimumCut/Cut0", parameter={"Cut": 0.0})
SIM.filter.mapDetFilter['InnerMPGDBarrel'] = "edep0"
SIM.filter.mapDetFilter['MPGDOuterBarrel'] = "edep0"
## SIM.filter.filters['edep5ev'] = dict(name="EnergyDepositMinimumCut/5eV", parameter={"Cut": 0.005*keV})
## SIM.filter.mapDetFilter['InnerMPGDBarrel'] = "edep5ev"
## SIM.filter.mapDetFilter['MPGDOuterBarrel'] = "edep5ev"
## SIM.filter.filters['edep1dev'] = dict(name="EnergyDepositMinimumCut/1deV", parameter={"Cut": 0.010*keV})
## SIM.filter.mapDetFilter['InnerMPGDBarrel'] = "edep1dev"
## SIM.filter.mapDetFilter['MPGDOuterBarrel'] = "edep1dev"
