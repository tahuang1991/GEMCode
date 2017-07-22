#ifndef GEMCode_GEMValidation_L1TrackTriggerVeto_cc
#define GEMCode_GEMValidation_L1TrackTriggerVeto_cc

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
/* #include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h" */
/* #include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h" */

class L1TrackTriggerVeto
{
 public:
  L1TrackTriggerVeto(const edm::ParameterSet& ps,
                     const edm::EventSetup& es,
                     const edm::Event& iEvent,
                     float eta, float phi);

  // TT Track veto
  void setEtaPhiReference(double eta, double phi) {etaReference_ = eta; phiReference_ = phi;}
  void calculateTTIsolation(const std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >&);
  bool isLooseVeto() const {return isLooseVeto_;}
  bool isMediumVeto() const {return isMediumVeto_;}
  bool isTightVeto() const {return isTightVeto_;}

 private:
  const edm::ParameterSet& ps_;
  const edm::Event& ev_;
  const edm::EventSetup& es_;
  float etaReference_;
  float phiReference_;

  std::vector<edm::InputTag> trackInput_;
  int verbose_;
  bool run_;

  GlobalPoint extrapolateGP(const TTTrack< Ref_Phase2TrackerDigi_ > &tk, int station=2);
  TrajectoryStateOnSurface propagateToZ(const GlobalPoint &, const GlobalVector &, double, double) const;
  TrajectoryStateOnSurface propagateToR(const GlobalPoint &, const GlobalVector &, double, double) const;

  // veto
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > tttracks_;
  bool isLooseVeto_;
  bool isMediumVeto_;
  bool isTightVeto_;


  // propagators
  edm::ESHandle<MagneticField> magfield_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<Propagator> propagatorOpposite_;
  edm::ESHandle<Propagator> propagatorAny_;
};

#endif
