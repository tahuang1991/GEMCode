import FWCore.ParameterSet.Config as cms

def customizeL1TrackME0ForMuonTDR(process):
  process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
  process.load('Configuration.StandardSequences.Reconstruction_cff')

  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TTTrack*_Level1TTTracks_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_L1TkMuons_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*Muon*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_me0*_*_*')

  # don't keep the track clusters and stubs
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTCluster*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_TTStub*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simEcal*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simHcal*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_mix_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiStripDigis_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiPixelDigis_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simCalo*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simCalo*_TrackerHits*_*')

  process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
  from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *

  ## customization provided by Sebastien Viret

  ## This crash is related to a geometry point. The branch you merged was
  ## based on a tracker geometry with 11 rings in the TIB1, whereas I suppose
  ## that in D17 you're using the last one, with 12 rings. Therefore the stub
  ## producer joboption is bad.

  # Overwrite stub builder config
  #
  # Baseline Threshold
  # 2Gev thresh 99% in the barrel flat and tilted except in TIB1 where 3GeV 99% and 95% in two last tilted rings
  # 2GeV thresh 99% in the disks except in ring1 1 and 2 of disk 1 and 3 where 95% 3GeV
  #
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(True)
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.zMatching2S = cms.bool(True)
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.NTiltedRings = cms.vdouble( 0., 12., 12., 12., 0., 0., 0.) #Number of tilted rings per side in barrel layers (for tilted geom only)
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.BarrelCut = cms.vdouble( 0, 2., 2., 3, 4, 5, 6.5)
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.TiltedBarrelCutSet = cms.VPSet(
    cms.PSet( TiltedCut = cms.vdouble( 0 ) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 2, 2, 2, 2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1., 1.) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 2.5, 2.5, 2.5, 2.5, 2, 2, 2.5, 2.5, 2, 2, 2, 2) ),
    cms.PSet( TiltedCut = cms.vdouble( 0, 3.5, 3.5, 3, 3, 3, 3, 2.5, 2.5, 2.5, 2, 2, 2) ),
    )
  process.TTStubAlgorithm_official_Phase2TrackerDigi_.EndcapCutSet = cms.VPSet(
    cms.PSet( EndcapCut = cms.vdouble( 0 ) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2.5, 3, 3, 3.5, 4, 2.5, 3, 3.5, 4.5, 5.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1, 1, 1.5, 2, 2, 2, 2.5, 3, 3, 3, 2, 3, 4, 5, 5.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 2, 2, 2.5, 2.5, 2.5, 3.5, 2.5, 5, 5.5, 6) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 2, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
    cms.PSet( EndcapCut = cms.vdouble( 0, 1.5, 1.5, 1.5, 1.5, 2, 2, 2, 3, 3, 6, 6, 6.5) ),
    )

  process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)

  # L1 tracking
  process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
  process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
  process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

  ## ME0 segments
  process.simMuonME0ReDigis192 = process.simMuonME0ReDigis.clone(
    numberOfStrips = cms.uint32(192)
  )
  process.RandomNumberGeneratorService.simMuonME0ReDigis192 = cms.PSet(
    initialSeed = cms.untracked.uint32(2234567),
    engineName = cms.untracked.string('HepJamesRandom')
  )
  process.me0RecHits192 = process.me0RecHits.clone(
    me0DigiLabel = cms.InputTag("simMuonME0ReDigis192")
  )
  process.me0Segments192 = process.me0Segments.clone(
    me0RecHitLabel = cms.InputTag("me0RecHits192")
  )
  process.me0DigiRecoSequence = cms.Sequence(
    process.simMuonME0ReDigis192 *
    process.me0RecHits192 *
    process.me0Segments192
  )

  process.me0Segments192.algo_psets[1].algo_pset.maxPhiAdditional = cms.double(1.2*0.35/192)
  process.me0Segments192.algo_psets[1].algo_pset.maxPhiSeeds = cms.double(1.2*0.35/192)

  process.L1simulation_step = cms.Path(process.SimL1TMuon*process.me0DigiRecoSequence)
  process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)

  process.schedule = cms.Schedule(process.digitisation_step,
                                  process.L1simulation_step,
                                  process.L1TrackTrigger_step,
                                  process.TTClusterStub,
                                  process.TTTracksWithTruth,
                                  process.endjob_step,
                                  process.FEVTDEBUGHLToutput_step)
  return process
