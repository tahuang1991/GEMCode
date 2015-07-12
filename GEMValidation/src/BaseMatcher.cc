#include "GEMCode/GEMValidation/interface/BaseMatcher.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"


BaseMatcher::BaseMatcher(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es)
: trk_(t), vtx_(v), conf_(ps), ev_(ev), es_(es), verbose_(0)
{
  // list of CSC chamber type numbers to use
  std::vector<int> csc_types = conf().getParameter<std::vector<int> >("useCSCChamberTypes");
  for (int i=0; i <= CSC_ME42; ++i) useCSCChamberTypes_[i] = false;
  for (auto t: csc_types)
  {
    if (t >= 0 && t <= CSC_ME42) useCSCChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (csc_types.empty()) useCSCChamberTypes_[CSC_ALL] = true;

  // list of RPC chamber type numbers to use
  std::vector<int> rpc_types = conf().getParameter<std::vector<int> >("useRPCChamberTypes");
  for (int i=0; i <= RPC_MB24n; ++i) useRPCChamberTypes_[i] = false;
  for (auto t: rpc_types)
  {
    if (t >= 0 && t <= RPC_MB24n) useRPCChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (rpc_types.empty()) useRPCChamberTypes_[RPC_ALL] = true;

  // list of DT chamber type numbers to use
  std::vector<int> dt_types = conf().getParameter<std::vector<int> >("useDTChamberTypes");
  for (int i=0; i <= DT_MB24n; ++i) useDTChamberTypes_[i] = false;
  for (auto t: dt_types)
  {
    if (t >= 0 && t <= DT_MB24n) useDTChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (dt_types.empty()) useDTChamberTypes_[DT_ALL] = true;

  // list of GEM chamber type numbers to use
  std::vector<int> gem_types = conf().getParameter<std::vector<int> >("useGEMChamberTypes");
  for (int i=0; i <= GEM_ME21; ++i) useGEMChamberTypes_[i] = false;
  for (auto t: gem_types)
  {
    if (t >= 0 && t <= GEM_ME21) useGEMChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (gem_types.empty()) useGEMChamberTypes_[GEM_ALL] = true;

  // Get the magnetic field
  es.get<IdealMagneticFieldRecord>().get(magfield_);

  // Get the propagators                                                                                  
  es.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);

  /// get the geometry
  hasGEMGeometry_ = true;
  hasRPCGeometry_ = true;
  hasCSCGeometry_ = true;
  hasME0Geometry_ = true;
  hasDTGeometry_ = true;

  try {
    es.get<MuonGeometryRecord>().get(gem_geom_);
    gemGeometry_ = &*gem_geom_;
  } catch (edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    hasGEMGeometry_ = false;
    LogDebug("BaseMatcher") << "+++ Info: GEM geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(me0_geom_);
    me0Geometry_ = &*me0_geom_;
  } catch (edm::eventsetup::NoProxyException<ME0Geometry>& e) {
    hasME0Geometry_ = false;
    LogDebug("BaseMatcher") << "+++ Info: ME0 geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(csc_geom_);
    cscGeometry_ = &*csc_geom_;
  } catch (edm::eventsetup::NoProxyException<CSCGeometry>& e) {
    hasCSCGeometry_ = false;
    LogDebug("BaseMatcher") << "+++ Info: CSC geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(rpc_geom_);
    rpcGeometry_ = &*rpc_geom_;
  } catch (edm::eventsetup::NoProxyException<RPCGeometry>& e) {
    hasRPCGeometry_ = false;
    LogDebug("BaseMatcher") << "+++ Info: RPC geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(dt_geom_);
    dtGeometry_ = &*dt_geom_;
  } catch (edm::eventsetup::NoProxyException<DTGeometry>& e) {
    hasDTGeometry_ = false;
    LogDebug("BaseMatcher") << "+++ Info: DT geometry is unavailable. +++\n";
  }

  simTrackPSet_ = conf().getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrackPSet_.getParameter<int>("verbose");
}


BaseMatcher::~BaseMatcher()
{
}


bool 
BaseMatcher::useGEMChamberType(int gem_type)
{
  if (gem_type < 0 || gem_type > GEM_ME21) return false;
  return useGEMChamberTypes_[gem_type];
}

bool 
BaseMatcher::useRPCChamberType(int rpc_type)
{
  if (rpc_type < 0 || rpc_type > RPC_MB24n) return false;
  return useRPCChamberTypes_[rpc_type];
}

bool 
BaseMatcher::useDTChamberType(int dt_type)
{
  if (dt_type < 0 || dt_type > DT_MB24n) return false;
  return useDTChamberTypes_[dt_type];
}

bool 
BaseMatcher::useCSCChamberType(int csc_type)
{
  if (csc_type < 0 || csc_type > CSC_ME42) return false;
  return useCSCChamberTypes_[csc_type];
}


GlobalPoint
BaseMatcher::propagateToZ(GlobalPoint &inner_point, GlobalVector &inner_vec, float z) const
{
  Plane::PositionType pos(0.f, 0.f, z);
  Plane::RotationType rot;
  Plane::PlanePointer my_plane(Plane::build(pos, rot));

  FreeTrajectoryState state_start(inner_point, inner_vec, trk_.charge(), &*magfield_);

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_plane));
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_plane);

  if (tsos.isValid()) return tsos.globalPosition();
  return GlobalPoint();
}


GlobalPoint
BaseMatcher::propagateToZ(float z) const
{
  GlobalPoint inner_point(vtx_.position().x(), vtx_.position().y(), vtx_.position().z());
  GlobalVector inner_vec (trk_.momentum().x(), trk_.momentum().y(), trk_.momentum().z());
  return propagateToZ(inner_point, inner_vec, z);
}


GlobalPoint
BaseMatcher::propagatedPositionGEM() const
{
  const double eta(trk().momentum().eta());
  const int endcap( (eta > 0.) ? 1 : -1);
  return propagateToZ(endcap*AVERAGE_GEM_Z);
}


double 
BaseMatcher::phiHeavyCorr(double pt, double eta, double phi, double charge) const
{
    // float resEta = eta;
    float etaProp = std::abs(eta);
    if (etaProp< 1.1) etaProp = 1.1;
    float resPhi = phi - 1.464*charge*cosh(1.7)/cosh(etaProp)/pt - M_PI/144.;
    if (resPhi > M_PI) resPhi -= 2.*M_PI;
    if (resPhi < -M_PI) resPhi += 2.*M_PI;
    return resPhi;
}


bool 
BaseMatcher::passDPhicut(CSCDetId id, float dPhi, float pt) const
{
  //  const double GEMdPhi[9][3];
  if (!(id.station()==1 and (id.ring()==1 or id.ring()==4)) &&
	!(id.station()==2 and id.ring()==1))  return true;
   
  auto GEMdPhi( id.station()==1 ? ME11GEMdPhi : ME21GEMdPhi);
   // std::copy(&ME11GEMdPhi[0][0], &ME11GEMdPhi[0][0]+9*3,&GEMdPhi[0][0]);
   //else if (id.station()==2 and id.ring()==1) 
   // std::copy(&ME21GEMdPhi[0][0], &ME21GEMdPhi[0][0]+9*3,&GEMdPhi[0][0]);
   
   bool is_odd(id.chamber()%2==1);
   bool pass = false;

   for (int b = 0; b < 9; b++)
   {
	if (double(pt) >= GEMdPhi[b][0])
	{
		
	    if ((is_odd && GEMdPhi[b][1] > fabs(dPhi)) ||
		(!is_odd && GEMdPhi[b][2] > fabs(dPhi)))
		    pass = true;
	    else    pass = false;
	}
    }
   if (dPhi < -50) pass = true;

   return pass;

}


