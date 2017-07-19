import FWCore.ParameterSet.Config as cms

def customizeL1TrackME0ForMuonTDR(process):
  process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
  process.load('Configuration.StandardSequences.Reconstruction_cff')

  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TTTrack*_Level1TTTracks_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TTCluster*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_TTStub*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_L1TkMuons_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*Muon*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('keep *_me0*_*_*')

  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simEcal*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simHcal*_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_mix_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiStripDigis_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simSiPixelDigis_*_*')
  process.FEVTDEBUGHLToutput.outputCommands.append('drop *_simCalo*_*_*')

  process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
  from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
  process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)

  # L1 tracking
  process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
  process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
  process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

  process.load("L1Trigger.L1TTrackMatch.L1TkMuonProducer_cfi")
  process.pL1TkMuon = cms.Path( process.L1TkMuons )

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

  process.L1simulation_step = cms.Path(process.SimL1Emulator*process.me0DigiRecoSequence)
  process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)

  process.schedule = cms.Schedule(process.digitisation_step,
                                  process.L1simulation_step,
                                  process.L1TrackTrigger_step,
                                  process.TTClusterStub,
                                  process.TTTracksWithTruth,
                                  process.endjob_step,
                                  process.FEVTDEBUGHLToutput_step)
  return process
