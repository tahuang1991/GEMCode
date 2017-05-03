#ifndef GEMCode_GEMValidation_TFTrack_h
#define GEMCode_GEMValidation_TFTrack_h

// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>

// user include files
#include "GEMCode/GEMValidation/interface/GenericDigi.h"
#include "GEMCode/GEMValidation/interface/BaseMatcher.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
/* #include "L1Trigger/CSCTrackFinder/interface/CSCTFSectorProcessor.h" */
/* #include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h" */
/* #include "L1Trigger/CSCTrackFinder/interface/CSCTrackFinderDataTypes.h" */
/* #include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h" */
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>
#include "L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

class TFTrack
{
 public:
  /// constructor
  /* TFTrack(const csc::L1Track* t, const CSCCorrelatedLCTDigiCollection*); */

  TFTrack(const l1t::EMTFTrack *t);
  /// copy constructor
  TFTrack(const TFTrack&);
  /// destructor
  ~TFTrack();

  //l1t::tftype = enum{bmtf, omtf_neg, omtf_pos, emtf_neg, emtf_pos};
  enum {CSCTF_Track, EMTF_Track, OMTF_Track};
  /* void init(edm::ESHandle< L1MuTriggerScales > &muScales, */
	/*     edm::ESHandle< L1MuTriggerPtScale > &muPtScale); */

  void setDR(double dr);

  /// L1 track
  /* const csc::L1Track* getL1Track() const {return l1track_;} */
  /// collection of trigger digis
  std::vector<const CSCCorrelatedLCTDigi*> getTriggerDigis() const {return triggerDigis_;}
  /// collection of MPC LCTs
  std::vector<CSCDetId> getTriggerDigisIds() const {return triggerIds_;}
  std::vector<std::pair<float, float>> getTriggerEtaPhis() {return triggerEtaPhis_;}
  //std::vector<csctf::TrackStub> getTriggerStubs() const {return triggerStubs_;}
  std::vector<matching::Digi*> getTriggerMPLCTs() const {return mplcts_;}
  std::vector<CSCDetId> getChamberIds() const {return ids_;}

  unsigned int digiInME(int st, int ring) const;

  void addTriggerDigi(const CSCCorrelatedLCTDigi*);
  void addTriggerDigiId(const CSCDetId&);
  void addTriggerEtaPhi(const std::pair<float,float>&);
  /* void addTriggerStub(const csctf::TrackStub&); */

  int tracktype() const { return trackType_; }
  unsigned int nStubs() const {return l1track_->NumHits();}
  /// how many stubs?
  /* unsigned int nStubs(bool mb1, bool me1, bool me2, bool me3, bool me4) const; */
  /* unsigned int nStubs() const {return nstubs;} */
  /* /// how many stubs in CSC? */
  /* unsigned int nStubsCSCOk(bool me1, bool me2, bool me3, bool me4) const; */
  /* /// has stub in muon barrel/endcap */
  /* bool hasStubStation(int wheel) const; */
  /* /// has stub in muon barrel? */
  /* bool hasStubBarrel() const; */
  /* /// has stub in muon endcap? */
  bool hasStubEndcap(int station) const;
  /* /// matches CSC stubs? */
  /* /\* bool hasStubCSCOk(int st) const; *\/ */
  /* /// has stubs that pass match? */
  /* bool passStubsMatch(double eta, int minLowHStubs, int minMidHStubs, int minHighHStubs) const */;

  /// print some information
  void print();



  /// bending angles
  unsigned dPhi12() const { return dPhi12_;}
  unsigned dPhi23() const { return dPhi23_;}

  unsigned ptPacked() const {return pt_packed_;}
  unsigned etaPacked() const {return eta_packed_;}
  unsigned phiPacked() const {return phi_packed_;}
  unsigned qPacked() const {return q_packed_;}
  unsigned int chargesign() const {return chargesign_;}
  int charge() const {return charge_;}
  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}
  double phi_local() const {return phi_local_;}
  double dr() const {return dr_;}
  int quality() const {return quality_;}
  int bx() const {return bx_;}
  std::vector<bool> deltaOk();
  bool debug() const {return debug_;}
  bool passDPhicutTFTrack(int st, float pt) const;

 private:
  int trackType_;
  const l1t::EMTFTrack* l1track_;
  std::vector<const CSCCorrelatedLCTDigi*> triggerDigis_;
  const l1t::EMTFHitCollection * trackHits_; //similar to triggerDigis_ + triggerDigis_
  std::vector<CSCDetId> triggerIds_;
  std::vector<std::pair<float, float>> triggerEtaPhis_;
  /* std::vector<csctf::TrackStub> triggerStubs_; */
  std::vector<matching::Digi*> mplcts_;
  std::vector<CSCDetId> ids_; // chamber ids
  unsigned phi_packed_;
  unsigned eta_packed_;
  unsigned pt_packed_;
  unsigned q_packed_;
  unsigned dPhi12_;
  unsigned dPhi23_;
  unsigned int chargesign_;
  int charge_;
  double phi_;
  double phi_local_;
  double eta_;
  double pt_;
  double dr_;
  int quality_;
  int bx_;
  unsigned int nstubs;
  std::vector<bool> deltaOk_;
  bool debug_;
};

#endif
