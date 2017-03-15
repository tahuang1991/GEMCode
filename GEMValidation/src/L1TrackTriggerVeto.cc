#include "GEMCode/GEMValidation/interface/L1TrackTriggerVeto.h"
#include "GEMCode/GEMValidation/interface/Helpers.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

L1TrackTriggerVeto::L1TrackTriggerVeto(const edm::ParameterSet& ps,
                                       const edm::EventSetup& es,
                                       const edm::Event& ev,
                                       float eta, float phi)
  : ps_(ps), ev_(ev), es_(es), etaReference_(eta), phiReference_(phi)
{
  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);

  isLooseVeto_ = 0;
  isMediumVeto_ = 0;
  isTightVeto_ = 0;

  auto l1track = ps_.getParameter<edm::ParameterSet>("l1track");
  trackInput_ = l1track.getParameter<std::vector<edm::InputTag>>("validInputTags");
  verbose_ = l1track.getParameter<int>("verbose");
  run_ = l1track.getParameter<bool>("run");

  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  if (gemvalidation::getByLabel(trackInput_, TTTrackHandle, ev_) and run_) {
    calculateTTIsolation(*TTTrackHandle.product());
  }
}


void L1TrackTriggerVeto::calculateTTIsolation(const std::vector< TTTrack< Ref_PixelDigi_ > >& TTTracks)
{
  for (unsigned int j=0; j<TTTracks.size(); ++j) {
    auto l1Tk = TTTracks[j];
    const double l1Tk_pt = l1Tk.getMomentum().perp();

    double l1Tk_eta_prop = -99;
    double l1Tk_phi_prop = -99;
    GlobalPoint ex_point(extrapolateGP(l1Tk));
    if (!(ex_point == GlobalPoint())) {
      l1Tk_eta_prop = ex_point.eta();
      l1Tk_phi_prop = normalizedPhi(ex_point.phi());
      const double dR_l1Mu_l1Tk_prop = reco::deltaR(l1Tk_eta_prop, l1Tk_phi_prop,
                                                    etaReference_, phiReference_);

      if (dR_l1Mu_l1Tk_prop <= 0.12 and l1Tk_pt >= 4) isLooseVeto_ = 1;
      if (dR_l1Mu_l1Tk_prop <= 0.12 and l1Tk_pt >= 3) isMediumVeto_ = 1;
      if (dR_l1Mu_l1Tk_prop <= 0.12 and l1Tk_pt >= 2) isTightVeto_ = 1;
    }
  }
}


GlobalPoint
L1TrackTriggerVeto::extrapolateGP(const TTTrack< Ref_PixelDigi_ > &tk, int station)
{
  TrajectoryStateOnSurface tsos;
  GlobalPoint inner_point(tk.getPOCA());
  GlobalVector inner_vec (tk.getMomentum());
  double charge(tk.getRInv()>0? 1: -1);
  double R, Zmin, Zmax;
  if (station == 1){
    R = 440.; Zmax = 600.; Zmin = -600.;
  }
  else if (station == 2){
    R = 523.; Zmax = 828.; Zmin = -828.;
  }
  else {
    R = 0.; Zmax = 0.; Zmin = 0.;
  }

  if (std::abs(tk.getMomentum().eta())<1.2) tsos = propagateToR(inner_point, inner_vec, charge, R);
  else if (tk.getMomentum().eta()>1.2)      tsos = propagateToZ(inner_point, inner_vec, charge, Zmax);
  else if (tk.getMomentum().eta()<-1.2)     tsos = propagateToZ(inner_point, inner_vec, charge, Zmin);
  else                                      tsos = TrajectoryStateOnSurface();

  if (tsos.isValid()) return tsos.globalPosition();
  else                return GlobalPoint();
}


TrajectoryStateOnSurface
L1TrackTriggerVeto::propagateToZ(const GlobalPoint &inner_point, const GlobalVector &inner_vec, double charge, double z) const
{
  Plane::PositionType pos(0.f, 0.f, z);
  Plane::RotationType rot;
  Plane::PlanePointer my_plane(Plane::build(pos, rot));

  FreeTrajectoryState state_start(inner_point, inner_vec, charge, &*magfield_);

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_plane));
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_plane);
  return tsos;
}


TrajectoryStateOnSurface
L1TrackTriggerVeto::propagateToR(const GlobalPoint &inner_point, const GlobalVector &inner_vec, double charge, double R) const
{
  Cylinder::CylinderPointer my_cyl(Cylinder::build(Surface::PositionType(0,0,0), Surface::RotationType(), R));

  FreeTrajectoryState state_start(inner_point, inner_vec, charge, &*magfield_);

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_cyl));
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_cyl);
  return tsos;
}
