import FWCore.ParameterSet.Config as cms

SimTrackMatching = cms.PSet(
    # common
    useCSCChamberTypes = cms.untracked.vint32(0,1,2,3,4,5,6,7,8,9,10),
    cscStations = cms.vstring('ALL','ME11','ME1a','ME1b',
                              'ME12','ME13','ME21','ME22',
                              'ME31','ME32','ME41','ME42'),
    ntupleTrackChamberDelta = cms.bool(True),
    ntupleTrackEff = cms.bool(True),
    overrideminNHitsChamber = cms.bool(False),
    minNHitsChamber = cms.untracked.int32(4),
    verbose = cms.bool(False),
    ## per collection params
    simTrack = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag('g4SimHits'),
        minPt = cms.double(1.5),
        maxPt = cms.double(999.),
        minEta = cms.double(1.45),
        maxEta = cms.double(4.0),
        onlyMuon = cms.bool(True),
        requireVertex = cms.bool(True),
        requireGenPart = cms.bool(True),
    ),
    ## GEM
    gemSimHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag('g4SimHits','MuonGEMHits'),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    gemStripDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonGEMDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    gemPadDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonGEMCSCPadDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
     ),
    gemCoPadDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonGEMCSCPadDigis", "Coincidence"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    gemRecHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("gemRecHits"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    ## ME0
    me0SimHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag('g4SimHits','MuonME0Hits'),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    me0DigiPreReco = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonME0Digis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    me0RecHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("me0RecHits"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    me0Segment = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("me0Segments"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    me0Muon = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("me0MuonConverter"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    ## RPC
    rpcSimHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag('g4SimHits','MuonRPCHits'),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    rpcStripDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonRPCDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    rpcRecHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("rpcRecHits"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    ## CSC
    cscSimHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag('g4SimHits','MuonCSCHits'),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
        minNHitsChamber = cms.int32(4),
    ),
    cscStripDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonCSCDigis", "MuonCSCComparatorDigi"),
        minBX = cms.int32(3),
        maxBX = cms.int32(9),
        matchDeltaStrip = cms.int32(1),
        minNHitsChamber = cms.int32(4),
    ),
    cscWireDigi = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simMuonCSCDigis", "MuonCSCWireDigi"),
        minBX = cms.int32(3),
        maxBX = cms.int32(8),
        matchDeltaWG = cms.int32(1),
        minNHitsChamber = cms.int32(4),
    ),
    cscCLCT = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCscTriggerPrimitiveDigis"),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
    ),
    cscALCT = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCscTriggerPrimitiveDigis"),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
    ),
    cscLCT = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCscTriggerPrimitiveDigis"),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
        addGhosts = cms.bool(True),
        matchAlctGem = cms.bool(False),
        hsFromSimHitMean = cms.bool(True),
    ),
    cscMPLCT = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCscTriggerPrimitiveDigis"),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
        addGhosts = cms.bool(True),
    ),
    cscRecHit = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("csc2DRecHits"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    cscSegment = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("cscSegments"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    ## tracks
    tfTrack = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCsctfTrackDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    tfCand = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simCsctfDigis", "CSC"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    gmtRegCand = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simGmtDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    gmtCand = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("simGmtDigis"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    l1Extra = cms.PSet(
        verbose = cms.int32(0),
        input = cms.InputTag("l1extraParticles"),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
)

## additional utilities
def noRPCCollections(analyzer):
    analyzer.simTrackMatching.rpcSimHit.input = ""
    analyzer.simTrackMatching.rpcStripDigi.input = ""
    analyzer.simTrackMatching.rpcRecHit.input = ""
    return analyzer

def noGEMCollections(analyzer):
    analyzer.simTrackMatching.gemSimHit.input = ""
    analyzer.simTrackMatching.gemStripDigi.input = ""
    analyzer.simTrackMatching.gemPadDigi.input = ""
    analyzer.simTrackMatching.gemCoPadDigi.input = ""
    analyzer.simTrackMatching.gemRecHit.input = ""
    return analyzer

def noCSCCollections(analyzer):
    analyzer.simTrackMatching.cscSimHit.input = ""
    analyzer.simTrackMatching.cscStripDigi.input = ""
    analyzer.simTrackMatching.cscWireDigi.input = ""
    analyzer.simTrackMatching.cscCLCT.input = ""
    analyzer.simTrackMatching.cscALCT.input = ""
    analyzer.simTrackMatching.cscLCT.input = ""
    analyzer.simTrackMatching.cscMPLCT.input = ""
    analyzer.simTrackMatching.cscRecHit.input = ""
    analyzer.simTrackMatching.cscSegment.input = ""
    return analyzer

def noME0Collections(analyzer):
    analyzer.simTrackMatching.me0SimHit.input = ""
    analyzer.simTrackMatching.me0DigiPreReco.input = ""
    analyzer.simTrackMatching.me0RecHit.input = ""
    analyzer.simTrackMatching.me0Segment.input = ""
    analyzer.simTrackMatching.me0Muon.input = ""
    return analyzer

def noTrackCollections(analyzer):
    analyzer.simTrackMatching.tfTrack.input = ""
    analyzer.simTrackMatching.tfCand.input = ""
    analyzer.simTrackMatching.gmtRegCand.input = ""
    analyzer.simTrackMatching.gmtCand.input = ""
    analyzer.simTrackMatching.l1Extra.input = ""
    return analyzer

def noSimHitCollections(analyzer):    
    analyzer.simTrackMatching.rpcSimHit.input = ""
    analyzer.simTrackMatching.gemSimHit.input = ""
    analyzer.simTrackMatching.cscSimHit.input = ""
    analyzer.simTrackMatching.me0SimHit.input = ""
    return analyzer
    
def noDigiCollections(analyzer):    
    analyzer.simTrackMatching.me0DigiPreReco.input = ""
    analyzer.simTrackMatching.rpcStripDigi.input = ""
    analyzer.simTrackMatching.gemStripDigi.input = ""
    analyzer.simTrackMatching.gemPadDigi.input = ""
    analyzer.simTrackMatching.gemCoPadDigi.input = ""
    analyzer.simTrackMatching.cscStripDigi.input = ""
    analyzer.simTrackMatching.cscWireDigi.input = ""
    return analyzer

def noL1Collections(analyzer):    
    analyzer.simTrackMatching.cscCLCT.input = ""
    analyzer.simTrackMatching.cscALCT.input = ""
    analyzer.simTrackMatching.cscLCT.input = ""
    analyzer.simTrackMatching.cscMPLCT.input = ""
    return analyzer

def noRecoCollections(analyzer):    
    analyzer.simTrackMatching.me0RecHit.input = ""
    analyzer.simTrackMatching.me0Segment.input = ""
    analyzer.simTrackMatching.me0Muon.input = ""
    analyzer.simTrackMatching.gemRecHit.input = ""
    analyzer.simTrackMatching.rpcRecHit.input = ""
    analyzer.simTrackMatching.cscRecHit.input = ""
    analyzer.simTrackMatching.cscSegment.input = ""
    return analyzer

def onlySimHitCollections(analyzer):
    analyzer = noDigiCollections(analyzer)
    analyzer = noL1Collections(analyzer)
    analyzer = noRecoCollections(analyzer)    
    analyzer = noTrackCollections(analyzer)
    return analyzer

def upToDigiCollections(analyzer):
    analyzer = noL1Collections(analyzer)
    analyzer = noRecoCollections(analyzer)
    analyzer = noTrackCollections(analyzer)
    return analyzer

def upToL1Collections(analyzer):
    analyzer = noRecoCollections(analyzer)
    analyzer = noTrackCollections(analyzer)
    return analyzer

def upToRecoCollections(analyzer):
    analyzer = noTrackCollections(analyzer)
    return analyzer
