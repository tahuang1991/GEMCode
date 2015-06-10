import FWCore.ParameterSet.Config as cms
from L1Trigger.CSCTrackFinder.csctfTrackDigis_cfi import*

SimTrackMatching = cms.PSet(
    # common
    useCSCChamberTypes = cms.untracked.vint32(0),
    ## endcap stations
    cscStations = cms.vstring('ALL','ME11','ME1a','ME1b','ME12','ME13',
                              'ME21','ME22','ME31','ME32','ME41','ME42'),
    gemStations = cms.vstring('GE11','GE21'),
    me0Stations = cms.vstring('ME0'),
    rpcStations = cms.vstring('RE12','RE13','RE22','RE23','RE31',
                              'RE32','RE33','RE41','RE42','RE43'),
    ntupleTrackChamberDelta = cms.bool(True),
    ntupleTrackEff = cms.bool(True),
    overrideminNHitsChamber = cms.bool(False),
    minNHitsChamber = cms.untracked.int32(4),
    verbose = cms.bool(True),
    ## per collection params
    simTrack = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits')),
        minPt = cms.double(1.5),
        maxPt = cms.double(999.),
        minEta = cms.double(0.0),
        maxEta = cms.double(2.4),
        onlyMuon = cms.bool(True),
        requireVertex = cms.bool(True),
        requireGenPart = cms.bool(True),
    ),
    ## GEM
    gemSimHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits','MuonGEMHits')),
        run = cms.bool(True),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    gemStripDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonGEMDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    gemPadDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonGEMPadDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
     ),
    gemCoPadDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCscTriggerPrimitiveDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    gemRecHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("gemRecHits")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    ## ME0
    me0SimHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits','MuonME0Hits')),
        run = cms.bool(True),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    me0DigiPreReco = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonME0Digis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    me0RecHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("me0RecHits")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    me0Segment = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("me0Segments")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    me0Muon = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("me0MuonConverter")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    ## RPC
    rpcSimHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits','MuonRPCHits')),
        run = cms.bool(True),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
    ),
    rpcStripDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonRPCDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    rpcRecHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("rpcRecHits")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaStrip = cms.int32(1),
    ),
    ## CSC
    cscSimHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits','MuonCSCHits')),
        run = cms.bool(True),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
        minNHitsChamber = cms.int32(4),
    ),
    cscStripDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonCSCDigis", "MuonCSCComparatorDigi")),
        run = cms.bool(True),
        minBX = cms.int32(3),
        maxBX = cms.int32(9),
        matchDeltaStrip = cms.int32(2),
        minNHitsChamber = cms.int32(4),
    ),
    cscWireDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonCSCDigis", "MuonCSCWireDigi")),
        run = cms.bool(True),
        minBX = cms.int32(3),
        maxBX = cms.int32(8),
        matchDeltaWG = cms.int32(2),
        minNHitsChamber = cms.int32(4),
    ),
    cscCLCT = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCscTriggerPrimitiveDigis"),
                                       cms.InputTag("hltMuonCSCDigis","MuonCSCCLCTDigi")),
        run = cms.bool(True),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
    ),
    cscALCT = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCscTriggerPrimitiveDigis"),
                                       cms.InputTag("hltMuonCSCDigis","MuonCSCALCTDigi")),
        run = cms.bool(True),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
    ),
    cscLCT = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCscTriggerPrimitiveDigis"),
                                       cms.InputTag("hltMuonCSCDigis","MuonCSCCorrelatedLCTDigi")),
        run = cms.bool(True),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
        addGhosts = cms.bool(False),
        matchAlctGem = cms.bool(True),
        matchAlctRpc = cms.bool(True),
        matchClctGem = cms.bool(False),
        matchClctRpc = cms.bool(False),
        hsFromSimHitMean = cms.bool(True),
    ),
    cscMPLCT = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCscTriggerPrimitiveDigis")),
        run = cms.bool(True),
        minBX = cms.int32(5),
        maxBX = cms.int32(7),
        minNHitsChamber = cms.int32(4),
        addGhosts = cms.bool(True),
    ),
    cscRecHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("csc2DRecHits")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    cscSegment = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("cscSegments")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    ## DT
    dtSimHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag('g4SimHits','MuonDTHits')),
        run = cms.bool(True),
        simMuOnly = cms.bool(True),
        discardEleHits = cms.bool(True),
        minNHitsChamber = cms.int32(4),
    ),
    dtDigi = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonDTDigis")),
        run = cms.bool(False),
        ## not sure which BX is the central one
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        matchDeltaWire = cms.int32(1),
        minNHitsChamber = cms.int32(4),
    ),
    dtStubs = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simDtTriggerPrimitiveDigis")),
        run = cms.bool(False),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        minNHitsChamber = cms.int32(4),
    ),
    dtRecHit = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("dt2DRecHits")),
        run = cms.bool(False),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
    ),
    ## TrackFinder tracks
    cscTfTrack = cms.PSet(
        verbose = cms.int32(0),
        run = cms.bool(True),
        validInputTags = cms.VInputTag(cms.InputTag("simCsctfTrackDigis")),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.5),
    ),
    dtTfTrack = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simDttfDigis", "DTTF")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),        
    ),
    rpcPAC = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simMuonRPCDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),        
    ),
    ## TrackFinder candidates
    cscTfCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simCsctfDigis", "CSC")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
    dtTfCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simDttfDigis", "DT")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
    rpcfTfCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simRpcTriggerDigis","RPCf")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
    rpcbTfCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simRpcTriggerDigis","RPCb")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),    
    sectorProcessor = csctfTrackDigis.SectorProcessor,
#       SRLUT = cms.PSet(
#			Binary = cms.untracked.bool(False),
#			ReadLUTs = cms.untracked.bool(False),
#			LUTPath = cms.untracked.string('./'),
#			UseMiniLUTs = cms.untracked.bool(True),
#		),
#        PTLUT = cms.PSet(
#    			LowQualityFlag = cms.untracked.uint32(4),
#			ReadPtLUT = cms.bool(False),
#			PtMethod = cms.untracked.uint32(32),
#	       ),
#	CoreLatency = cms.uint32(7),
#	gangedME1a = cms.untracked.bool(True),
#	MinBX = cms.int32(3),
#	MaxBX = cms.int32(9),
#	initializeFromPSet = cms.bool(True),
#    ),
    ## GMT and L1Extra
    gmtRegCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simGmtDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
    gmtCand = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(cms.InputTag("simGmtDigis")),
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
    l1Extra = cms.PSet(
        verbose = cms.int32(0),
        validInputTags = cms.VInputTag(), #cms.InputTag("hltL1extraParticles")
        run = cms.bool(True),
        minBX = cms.int32(-1),
        maxBX = cms.int32(1),
        deltaR = cms.double(0.05),
    ),
)

## additional utilities
def noRPCCollections(analyzer):
    analyzer.simTrackMatching.rpcSimHit.validInputTags = ""
    analyzer.simTrackMatching.rpcStripDigi.validInputTags = ""
    analyzer.simTrackMatching.rpcRecHit.validInputTags = ""
    return analyzer

def noGEMCollections(analyzer):
    analyzer.simTrackMatching.gemSimHit.validInputTags = ""
    analyzer.simTrackMatching.gemStripDigi.validInputTags = ""
    analyzer.simTrackMatching.gemPadDigi.validInputTags = ""
    analyzer.simTrackMatching.gemCoPadDigi.validInputTags = ""
    analyzer.simTrackMatching.gemRecHit.validInputTags = ""
    return analyzer

def noCSCCollections(analyzer):
    analyzer.simTrackMatching.cscSimHit.validInputTags = ""
    analyzer.simTrackMatching.cscStripDigi.validInputTags = ""
    analyzer.simTrackMatching.cscWireDigi.validInputTags = ""
    analyzer.simTrackMatching.cscCLCT.validInputTags = ""
    analyzer.simTrackMatching.cscALCT.validInputTags = ""
    analyzer.simTrackMatching.cscLCT.validInputTags = ""
    analyzer.simTrackMatching.cscMPLCT.validInputTags = ""
    analyzer.simTrackMatching.cscRecHit.validInputTags = ""
    analyzer.simTrackMatching.cscSegment.validInputTags = ""
    return analyzer

def noME0Collections(analyzer):
    analyzer.simTrackMatching.me0SimHit.validInputTags = ""
    analyzer.simTrackMatching.me0DigiPreReco.validInputTags = ""
    analyzer.simTrackMatching.me0RecHit.validInputTags = ""
    analyzer.simTrackMatching.me0Segment.validInputTags = ""
    analyzer.simTrackMatching.me0Muon.validInputTags = ""
    return analyzer

def noTrackCollections(analyzer):
    analyzer.simTrackMatching.cscTfTrack.validInputTags = ""
    analyzer.simTrackMatching.cscTfCand.validInputTags = ""
    analyzer.simTrackMatching.gmtRegCand.validInputTags = ""
    analyzer.simTrackMatching.gmtCand.validInputTags = ""
    analyzer.simTrackMatching.l1Extra.validInputTags = ""
    return analyzer

def noSimHitCollections(analyzer):    
    analyzer.simTrackMatching.rpcSimHit.validInputTags = ""
    analyzer.simTrackMatching.gemSimHit.validInputTags = ""
    analyzer.simTrackMatching.cscSimHit.validInputTags = ""
    analyzer.simTrackMatching.me0SimHit.validInputTags = ""
    return analyzer
    
def noDigiCollections(analyzer):    
    analyzer.simTrackMatching.me0DigiPreReco.validInputTags = ""
    analyzer.simTrackMatching.rpcStripDigi.validInputTags = ""
    analyzer.simTrackMatching.gemStripDigi.validInputTags = ""
    analyzer.simTrackMatching.gemPadDigi.validInputTags = ""
    analyzer.simTrackMatching.gemCoPadDigi.validInputTags = ""
    analyzer.simTrackMatching.cscStripDigi.validInputTags = ""
    analyzer.simTrackMatching.cscWireDigi.validInputTags = ""
    return analyzer

def noL1Collections(analyzer):    
    analyzer.simTrackMatching.cscCLCT.validInputTags = ""
    analyzer.simTrackMatching.cscALCT.validInputTags = ""
    analyzer.simTrackMatching.cscLCT.validInputTags = ""
    analyzer.simTrackMatching.cscMPLCT.validInputTags = ""
    return analyzer

def noRecoCollections(analyzer):    
    analyzer.simTrackMatching.me0RecHit.validInputTags = ""
    analyzer.simTrackMatching.me0Segment.validInputTags = ""
    analyzer.simTrackMatching.me0Muon.validInputTags = ""
    analyzer.simTrackMatching.gemRecHit.validInputTags = ""
    analyzer.simTrackMatching.rpcRecHit.validInputTags = ""
    analyzer.simTrackMatching.cscRecHit.validInputTags = ""
    analyzer.simTrackMatching.cscSegment.validInputTags = ""
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
