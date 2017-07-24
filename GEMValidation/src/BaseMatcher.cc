#include "GEMCode/GEMValidation/interface/BaseMatcher.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "GEMCode/GEMValidation/interface/Helpers.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>
//#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"


BaseMatcher::BaseMatcher(const SimTrack& t, const SimVertex& v,
      const edm::ParameterSet& ps, const edm::Event& ev, const edm::EventSetup& es)
: trk_(t), vtx_(v), conf_(ps), ev_(ev), es_(es), verbose_(0)
{
  // list of CSC chamber type numbers to use
  std::vector<int> csc_types = conf().getParameter<std::vector<int> >("cscStationsToUse");
  for (int i = CSC_ALL; i != CSC_ME42; i++) useCSCChamberTypes_[i] = false;
  for (auto t: csc_types)
  {
    if (t >= 0 && t <= CSC_ME42) useCSCChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (csc_types.empty()) useCSCChamberTypes_[CSC_ALL] = true;

  // list of RPC chamber type numbers to use
  std::vector<int> rpc_types = conf().getParameter<std::vector<int> >("rpcStationsToUse");
  for (int i = RPC_ALL; i != RPC_MB24n; i++) useRPCChamberTypes_[i] = false;
  for (auto t: rpc_types)
  {
    if (t >= 0 && t <= RPC_MB24n) useRPCChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (rpc_types.empty()) useRPCChamberTypes_[RPC_ALL] = true;

  // list of DT chamber type numbers to use
  std::vector<int> dt_types = conf().getParameter<std::vector<int> >("dtStationsToUse");
  for (int i = DT_ALL; i != DT_MB24n; i++) useDTChamberTypes_[i] = false;
  for (auto t: dt_types)
  {
    if (t >= 0 && t <= DT_MB24n) useDTChamberTypes_[t] = true;
  }
  // empty list means use all the chamber types
  if (dt_types.empty()) useDTChamberTypes_[DT_ALL] = true;

  // list of GEM chamber type numbers to use
  std::vector<int> gem_types = conf().getParameter<std::vector<int> >("gemStationsToUse");
  for (int i = GEM_ALL; i != GEM_ME21; i++) useGEMChamberTypes_[i] = false;
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
    std::cout << "+++ Info: GEM geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(me0_geom_);
    me0Geometry_ = &*me0_geom_;
  } catch (edm::eventsetup::NoProxyException<ME0Geometry>& e) {
    hasME0Geometry_ = false;
    std::cout << "+++ Info: ME0 geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(csc_geom_);
    cscGeometry_ = &*csc_geom_;
  } catch (edm::eventsetup::NoProxyException<CSCGeometry>& e) {
    hasCSCGeometry_ = false;
    std::cout << "+++ Info: CSC geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(rpc_geom_);
    rpcGeometry_ = &*rpc_geom_;
  } catch (edm::eventsetup::NoProxyException<RPCGeometry>& e) {
    hasRPCGeometry_ = false;
    std::cout << "+++ Info: RPC geometry is unavailable. +++\n";
  }

  try {
    es.get<MuonGeometryRecord>().get(dt_geom_);
    dtGeometry_ = &*dt_geom_;
  } catch (edm::eventsetup::NoProxyException<DTGeometry>& e) {
    hasDTGeometry_ = false;
    std::cout << "+++ Info: DT geometry is unavailable. +++\n";
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

float
BaseMatcher::lxy() const
{
  // calculate the simTrack dxy from the vertex and the momentum!!

  //Source: https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrackReco/interface/TrackBase.h#119
  const float vx(vtx_.position().x());
  const float vy(vtx_.position().y());

  return sqrt(vx*vx+vy*vy);
}


float
BaseMatcher::dxy() const
{
  // calculate the simTrack dxy from the vertex and the momentum!!

  //Source: https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/TrackReco/interface/TrackBase.h#119
  const float vx(vtx_.position().x());
  const float vy(vtx_.position().y());
  const float px(trk_.momentum().x());
  const float py(trk_.momentum().y());
  const float pt(trk_.momentum().pt());

  return (- vx * py + vy * px ) / pt;
}



GlobalPoint
BaseMatcher::propagateToZCharge(GlobalPoint &inner_point, GlobalVector &inner_vec, float z, int charge) const
{
  Plane::PositionType pos(0.f, 0.f, z);
  Plane::RotationType rot;
  Plane::PlanePointer my_plane(Plane::build(pos, rot));

  FreeTrajectoryState state_start(inner_point, inner_vec, charge, &*magfield_);
  //std::cout <<"state_start  position "<< state_start.position()<<" momentum "<< state_start.momentum()<<" charge "<<state_start.charge() << std::endl;
  if (state_start.hasError()) std::cout <<"state_start has error  "<< std::endl;

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_plane));
  if (!tsos.isValid()) std::cout <<" propagateToZCharge tsos not valid "<< std::endl;
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_plane);
  //broken here  when propagating ME0
  //std::cout <<"propagateToZ GP "<< tsos.globalPosition() <<" eta "<< tsos.globalPosition().eta()<<" phi "<< tsos.globalPosition().phi()<< std::endl;

  if (tsos.isValid()) return tsos.globalPosition();
  return GlobalPoint();
}


GlobalPoint
BaseMatcher::propagateToZ(GlobalPoint &inner_point, GlobalVector &inner_vec, float z) const
{
  Plane::PositionType pos(0.f, 0.f, z);
  Plane::RotationType rot;
  Plane::PlanePointer my_plane(Plane::build(pos, rot));

  FreeTrajectoryState state_start(inner_point, inner_vec, trk_.charge(), &*magfield_);
  //std::cout <<"state_start  position "<< state_start.position()<<" momentum "<< state_start.momentum()<<" charge "<<state_start.charge() << std::endl;
  if (state_start.hasError()) std::cout <<"state_start has error  "<< std::endl;

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_plane));
  if (!tsos.isValid()) std::cout <<" propagateToZ tsos not valid "<< std::endl;
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_plane);
  //broken here  when propagating ME0
  //std::cout <<"propagateToZ GP "<< tsos.globalPosition() <<" eta "<< tsos.globalPosition().eta()<<" phi "<< tsos.globalPosition().phi()<< std::endl;

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
BaseMatcher::propagateToR(GlobalPoint &inner_point, GlobalVector &inner_vec, float R) const
{
  Cylinder::CylinderPointer my_cyl(Cylinder::build(Surface::PositionType(0,0,0), Surface::RotationType(), R));

  FreeTrajectoryState state_start(inner_point, inner_vec, trk_.charge(), &*magfield_);

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_cyl));
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_cyl);

  if (tsos.isValid()) return tsos.globalPosition();
  return GlobalPoint();
}

GlobalPoint
BaseMatcher::propagateToR(float r) const
{
  GlobalPoint inner_point(vtx_.position().x(), vtx_.position().y(), vtx_.position().z());
  GlobalVector inner_vec (trk_.momentum().x(), trk_.momentum().y(), trk_.momentum().z());
  return propagateToR(inner_point, inner_vec, r);
}

GlobalPoint
BaseMatcher::propagatedPositionGEM() const
{
  const double eta(trk().momentum().eta());
  const int endcap( (eta > 0.) ? 1 : -1);
  return propagateToZ(endcap*AVERAGE_GEM_Z);
}

GlobalPoint
BaseMatcher::propagatedPositionSt2() const
{
  const double eta(trk().momentum().eta());
  const int endcap( (eta > 0.) ? 1 : -1);
  float AVERAGE_Z = (AVERAGE_ME21_EVEN_Z+AVERAGE_ME21_ODD_Z)/2.0;
  if (abs(eta)<1.2) return propagateToR(AVERAGE_DT2_R);
  else return  propagateToZ(endcap*AVERAGE_Z);
}


GlobalPoint
BaseMatcher::propagateFromME0ToCSC(ME0Segment segment,float pt, int charge, int st, bool evenodd) const
{
  int chamber = (evenodd? 1:2);
  int ring = 1;
  ME0DetId me0Id(segment.me0DetId());
  int endcap = (me0Id.region()>0) ? 1 : 2;
  auto me0Chamber(me0Geometry_->chamber(me0Id));
  const CSCDetId csclayerId(endcap, st, ring, chamber,  CSCConstants::KEY_CLCT_LAYER);
  //std::cout <<"ME0Id "<< me0Id <<" CSC endcap "<< me0Id.region() <<" st "<< st <<" ring "<< ring <<" chamber "<< chamber<<" cscid "<< csclayerId << std::endl;
  const CSCLayer* csclayer(cscGeometry_->layer(csclayerId));
  const GlobalPoint gp_csc = csclayer->centerOfWireGroup(10);
  //float propagate_z = (evenodd ? AVERAGE_ME11_ODD_Z : AVERAGE_ME11_EVEN_Z) * me0Id.region();
  float propagate_z = gp_csc.z();
  LocalPoint lp(segment.localPosition());
  GlobalPoint SegPos(me0Chamber->toGlobal(lp));
  LocalVector lv(segment.localDirection());
  GlobalVector SegDir(me0Chamber->toGlobal(lv));
  float ratio = pt/SegDir.perp();
  GlobalVector SegVec(SegDir.x()*ratio, SegDir.y()*ratio, SegDir.z()*ratio);
  //std::cout <<"CSC z "<< propagate_z <<" ME0 GP "<< SegPos << " GV "<< SegVec <<" GV.mag "<< SegVec.mag() <<" ratio "<< ratio << std::endl;
  return propagateToZCharge(SegPos, SegVec, propagate_z, charge);

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



bool BaseMatcher::passDPhicut(CSCDetId id, int chargesign, float dphi, float pt) const
//BaseMatcher::passDPhicut(CSCDetId id, float dPhi, float pt) const
{
  //  const double GEMdPhi[9][3];
    //std::cout <<"passdphicut  id " << id << std::endl;
  if (!(id.station()==1 and (id.ring()==1 or id.ring()==4)) &&
	!(id.station()==2 and id.ring()==1))  return true;
  auto GEMdPhi( id.station()==1 ? ME11GEMdPhi : ME21GEMdPhi);
   // std::copy(&ME11GEMdPhi[0][0], &ME11GEMdPhi[0][0]+9*3,&GEMdPhi[0][0]);
   //else if (id.station()==2 and id.ring()==1)
   // std::copy(&ME21GEMdPhi[0][0], &ME21GEMdPhi[0][0]+9*3,&GEMdPhi[0][0]);
  int st = id.station();
  bool is_odd(id.chamber()%2==1);
  bool pass = false;
  unsigned int LUTsize = (st==1)? sizeof(ME11GEMdPhi)/sizeof(ME11GEMdPhi[0]) :sizeof(ME21GEMdPhi)/sizeof(ME21GEMdPhi[0]);
  bool smalldphi = ((is_odd and fabs(dphi)<GEMdPhi[LUTsize-2][1]) || (!is_odd and fabs(dphi)<GEMdPhi[LUTsize-2][2]));
  if (fabs(dphi) < 99 and ((chargesign == 1 and dphi < 0) || (chargesign == 0 and dphi > 0) || smalldphi)){
   for (unsigned int b = 0; b < LUTsize; b++)
   {
//if (st==2) std::cout <<"BaseMatcher LUTpt "<<GEMdPhi[b][0] << " odd " << GEMdPhi[b][1]  <<" even " << GEMdPhi[b][2] <<" dphi "<<dphi <<std::endl;
	if (double(pt) >= GEMdPhi[b][0])
	{

	    if ((is_odd && GEMdPhi[b][1] > fabs(dphi)) ||
		(!is_odd && GEMdPhi[b][2] > fabs(dphi)))
		    pass = true;
	    else    pass = false;
	}
    }
  }
  else pass = false;

 // if (id.station() == 2 and pass) std::cout <<" st2 pass dphicut after comparing with LUT " << std::endl;
 // else if (id.station() == 2 and !pass) std::cout <<" st2 failed to pass dphicut after comparing " << std::endl;

  if (st==2 and pt>=15 and ((is_odd and fabs(dphi)<GEMdPhi[4][1]) || (!is_odd and fabs(dphi)<GEMdPhi[4][2]))) pass = true;

 // if (id.station() == 2 and pass) std::cout <<" st2 pass dphicut after pt>15 " << std::endl;
//  else if (id.station() == 2 and !pass) std::cout <<" st2 failed to pass dphicut after pt>15" << std::endl;


   return pass;

}


const CSCLayerGeometry*
BaseMatcher::retriveCSCKeyLayerGeometry(int rawid) const
{
  const CSCChamber* cscChamber(cscGeometry_->chamber(CSCDetId(rawid)));
  const CSCLayer* cscKeyLayer(cscChamber->layer(3));
  const CSCLayerGeometry* cscKeyLayerGeometry(cscKeyLayer->geometry());
  return cscKeyLayerGeometry;


}
