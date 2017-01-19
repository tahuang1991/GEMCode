import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomPtGunProducer",
  PGunParameters = cms.PSet(
    MaxPt = cms.double(100.0),
    MinPt = cms.double(1.0),
    PartID = cms.vint32(-13),
    MinPhi = cms.double(-3.14159265359),
    MaxPhi = cms.double(3.14159265359),
    MinEta = cms.double(-3.0),
    MaxEta = cms.double(3.0)
  ),
  Verbosity = cms.untracked.int32(0),
  psethack = cms.string('single mu pt 40 forward'),
  AddAntiParticle = cms.bool(True),
  firstRun = cms.untracked.uint32(1)
)
