#ifndef GEMCode_GEMValidation_DisplacedMuonTriggerPtassignment_cc
#define GEMCode_GEMValidation_DisplacedMuonTriggerPtassignment_cc

/**\class DisplacedMuonTriggerPtassignment

  Displaced Muon Trigger Design based on Muon system

  Author: tao.huang@cern.ch, sven.dildick@cern.ch

*/
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "GEMCode/GEMValidation/interface/DisplacedMuonTriggerPtassignment.h"
#include "GEMCode/GEMValidation/interface/PtassignmentHelper.h"
#include "GEMCode/GEMValidation/interface/BarrelTriggerPtAssignmentHelper.h"

using namespace std;

//endcap, we need LCTs and associated cscid, gempads and associated gemid, and all gemometry
//to get position from fitting, we also need all comparator digis
//step0 get LCTs and associated cscids, GEMPads and associated gemids, and geometry.
//step1 get fitting positions from fitting compara digis after assoicating comparator digis to LCTs
//step2 calculate all variables used pt in assignment, requires eta,phi,radius,Z
//step3 assgin L1 pt according to LUTs (in short future)
DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const CSCCorrelatedLCTDigiCollection* lcts, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);
}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const CSCCorrelatedLCTDigiContainer lcts, const CSCDetIdContainer cscids, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);
}

//chamberid_lcts: LCTs matched to simmuon and their associated chamberid, detid_pads: gempads matched to simmuon and their associated detid_pads
DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamberid_lcts, std::map<unsigned int, GEMCSCPadDigiContainer> detid_pads, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  chamberid_lcts_ = chamberid_lcts;
  detid_pads_ = detid_pads;
  ev.getByLabel("simMuonCSCDigis", "MuonCSCComparatorDigi", hCSCComparators);

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);

  initVariables();
  for (auto idlcts : chamberid_lcts_){
    CSCDetId chid(idlcts.first);
    //check the ring number that muon is flying through, 2nd station as reference
    //later use this one to check whether we should use GE21 or ME21only
    if (chid.station()==2) meRing = chid.ring();
    //chamber parity in station3 should be the same as in station2,
    //if there is already qualified one, the skip stubs in station 3
    if (chid.station()==3 and hasStub_st[1] and hasStub_st[2] and
	   ((isEven[1] and isEven[2]) or (not(isEven[1]) and not(isEven[2]))))
	continue;
    if (chid.chamber()%2 == 0) isEven[chid.station()-1] = true;
    if (idlcts.second.size()>0) hasStub_st[chid.station()-1] =true;
    else {
	    std::cout <<" chid "<< chid <<"  number of lcts "<< idlcts.second.size() << std::endl;
	    continue;
    }
    if (idlcts.second.size()>1 and verbose_>0)
	std::cout <<"more than one LCT  available in chamber id "<< chid <<" LCT size "<< idlcts.second.size()<<std::endl;
    globalPositionOfLCT(idlcts.second, chid);
    if (chid.station() == 1 or chid.station()==2){
      //find GEMPads
      for (auto idgempads : detid_pads){
        GEMDetId gemid(idgempads.first);
        if (((chid.station() == 1 and gemid.station() == 1) or (chid.station()==2 and gemid.station() ==3))
            and chid.chamber() == gemid.chamber()){
          //if gp_ge11 or gp_ge21 are taken from GME pad in layer1, then ignore the layer2
          if (hasGEMPad_st1 and gemid.station()==1 and gemid.layer()==2)  continue;
          if (hasGEMPad_st2 and gemid.station()==3 and gemid.layer()==2)  continue;
          if (gemid.station() == 1 ) hasGEMPad_st1 = true;
          else if (gemid.station() == 3) hasGEMPad_st2 = true;
          else if (verbose_>0)
            std::cout <<" gemid "<< gemid <<" CSC chamber id "<< chid << std::endl;
          //maybe also check the dR(csc, gem)
          globalPositionOfGEMPad(idgempads.second[0], gemid);
        }
      }
    }
  }
  nstubs = 0;
  for (int i=0; i<4; i++)
  	if (hasStub_st[i]){
	    nstubs++;
	    radius_st_ME[i] = gp_st_layer3[i].perp();
	}



  //npar>=0 is the flag to do pt assignment
  if (hasStub_st[0] and hasStub_st[1] and hasStub_st[2]){
  	if (not(isEven[0]) and isEven[1] and isEven[2]) npar = 0;
    else if (not(isEven[0]) and not(isEven[1]) and not(isEven[2])) npar = 1;
    else if (isEven[0] and isEven[1] and isEven[2]) npar = 2;
    else if (isEven[0] and not(isEven[1]) and not(isEven[2])) npar = 3;
    else {
	    std::cout <<" hasStub in station 1 2 3  but npar = -1 "<< std::endl;
	    npar= -1;
      for (auto idlcts : chamberid_lcts_){
        CSCDetId chid(idlcts.first);
        std::cout <<"CSC id "<< chid <<" LCT "<< idlcts.second[0] << std::endl;
	    }
    }
  }else

  npar = -1;
  if (nstubs >= 3){
  	//find fitting radius
     fitTrackRadius(gp_st_layer3, radius_st_ME);
  	//reset gp_st_layer3
     for (int i=0; i<4; i++){
  	if (not(hasStub_st[i])) continue;
	float phi = gp_st_layer3[i].phi();
	float z = gp_st_layer3[i].z();
	gp_st_layer3[i] = GlobalPoint(GlobalPoint::Cylindrical(radius_st_ME[i], phi, z));
     }
  }

  //second station
  if (hasStub_st[1]) eta_st2 = gp_st_layer3[1].eta();

  if (hasStub_st[0] and hasStub_st[1])
	xfactor = (gp_st_layer3[1].perp()/gp_st_layer3[0].perp()-1.0)/fabs(gp_st_layer3[0].z()-gp_st_layer3[1].z());

}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const L1CSCTrack& tftrack,
                                                                   const L1CSCTrackCollection& l1Tracks,
                                                                   const CSCCorrelatedLCTDigiCollection& CSCCorrelatedLCTs,
                                                                   bool doStubRecovery,
                                                                   bool matchGEMPads,
                                                                   const edm::EventSetup& es,
                                                                   const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  initVariables();
  setupTriggerScales(es);

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);

  // first step: collect all stubs associated to the CSC TF Track
  std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamberid_lct;
  auto stubCollection = tftrack.second;
  for (auto detUnitIt = stubCollection.begin(); detUnitIt != stubCollection.end(); detUnitIt++) {
    const CSCDetId& ch_id = (*detUnitIt).first;
    const auto range = (*detUnitIt).second;
    // empty vector for stubs
    CSCCorrelatedLCTDigiContainer v;
    for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
      auto stub = *digiIt;
      v.push_back(stub);
    }
    hasStub_st[ch_id.station()-1] = true;
    chamberid_lct[ch_id.rawId()] = v;
  }


  // TF properties
  auto track = tftrack.first;
  const int sign(track.endcap()==1 ? 1 : -1);
  unsigned gpt = 0, quality = 0;
  csc::L1Track::decodeRank(track.rank(), gpt, quality);

  //const float csctf_pt = muPtScale_->getPtScale()->getLowEdge(gpt & 0x1f) + 1.e-6;
  const float csctf_eta = muScales_->getRegionalEtaScale(2)->getCenter(track.eta_packed()) * sign;
  //const float csctf_phi = normalizedPhi(muScales_->getPhiScale()->getLowEdge(phiL1CSCTrack(track)));
  const int csctf_bx = track.bx();

  // second step: stub recovery
  bool stubMissingSt1 = not hasStub_st[0];
  bool stubMissingSt2 = not hasStub_st[1];
  bool stubMissingSt3 = not hasStub_st[2];
  bool stubMissingSt4 = not hasStub_st[3];
  bool atLeast1StubMissing = stubMissingSt1 or stubMissingSt2 or stubMissingSt3 or stubMissingSt4;

  if (doStubRecovery and atLeast1StubMissing){
    int triggerSector = track.sector();

    for (int endcap=1; endcap<=2; endcap++){
      // do not consider stubs in the wrong endcap
      int zendcap(endcap!=1 ? -1 : +1 );
      if (zendcap * csctf_eta < 0) continue;
      for (int station=1; station<=4; station++){

        // ignore station where a L1Mu stub is present!
        if (not stubMissingSt1 and station==1) continue;
        if (not stubMissingSt2 and station==2) continue;
        if (not stubMissingSt3 and station==3) continue;
        if (not stubMissingSt4 and station==4) continue;
        // if(verbose) std::cout << "Recovered stubs in station: " << station << std::endl;
        // temp storage of candidate stubs per station and ring
        CSCCorrelatedLCTDigiId bestMatchingStub;
        int iStub = 0;
        for (int ring=1; ring<=3; ring++){
          if (station!=1 and ring==3) continue;
          //std::cout << "Analyzing ring " << ring << std::endl;

          for (int chamber=1; chamber<=36; chamber++){
            // do not consider invalid detids
            if ( (station==2 or station==3 or station==4) and
                 (ring==1) and chamber>18) continue;
            //std::cout << "Analyzing chamber " << chamber << std::endl;
            // create the detid
            CSCDetId ch_id(endcap, station, ring, chamber);
            //std::cout << "ch_id " <<  ch_id << std::endl;
            // get the stubs in this detid
            auto range = CSCCorrelatedLCTs.get(ch_id);
            for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
              iStub++;

              auto stub(*digiItr);

              // check that this stub is not already part of a CSC TF track
              if (stubInCSCTFTracks(stub, l1Tracks)) continue;

              // trigger sector must be the same
              if (triggerSector != ch_id.triggerSector()) continue;

              // BXs have to match
              int deltaBX = std::abs(stub.getBX() - (6 + csctf_bx));
              if (deltaBX > 1) continue;

              if(verbose()>2) std::cout << ch_id << std::endl;
              if(verbose()>2) std::cout<<"Candidate " << stub << std::endl;
               // bestMatchingStub = pickBestMatchingStub(allxs[ch_id.station()-1], allys[ch_id.station()-1],
               //                                         bestMatchingStub, std::make_pair(ch_id, stub), 6 + csctf_bx);
             }
             // consider the case ME1a
             if (station==1 and ring==1){
               CSCDetId me1a_id(endcap, station, 4, chamber);
               auto range = CSCCorrelatedLCTs.get(me1a_id);
               for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
                 iStub++;
                 auto stub(*digiItr);

                 // check that this stub is not already part of a CSC TF track
                 if (stubInCSCTFTracks(stub, l1Tracks)) continue;

                 // trigger sector must be the same
                 if (triggerSector != me1a_id.triggerSector()) continue;

                 // BXs have to match
                 int deltaBX = std::abs(stub.getBX() - (6 + csctf_bx));
                 if (deltaBX > 1) continue;

                 if(verbose()>2) std::cout << me1a_id << std::endl;
                 if(verbose()>2) std::cout<<"Candidate " << stub << std::endl;
                 // bestMatchingStub = pickBestMatchingStub(allxs[me1a_id.station()-1], allys[me1a_id.station()-1],
                 //                                         bestMatchingStub, std::make_pair(me1a_id, stub), 6 + csctf_bx);
               }
             }
          }
        }
      }
    }
  }

  // third step: add GEM pads
  if (matchGEMPads) {

    // GEMCSCPadDigiId bestPad_GE11_L1;
    // GEMCSCPadDigiId bestPad_GE11_L2;
    // GEMCSCPadDigiId bestPad_GE21_L1;
    // GEMCSCPadDigiId bestPad_GE21_L2;

    // FIXME
  }
}


DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(GlobalPoint gp1, GlobalPoint gp2, GlobalPoint gp3, GlobalPoint gp4, GlobalPoint gp5, GlobalPoint gp6, int npar_in, const edm::EventSetup& es, const edm::Event& ev)
: ev_(ev), es_(es), verbose_(0)
{ //sim level

  setupGeometry(es);
  initVariables();

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);
  gp_st_layer3[0] = GlobalPoint(gp1);
  //use z>100 to make sure this is valid globalpoint
  if (fabs(gp1.z())>100) hasStub_st[0] = true;
  else if (verbose_>0)
      std::cout <<" gp_st1 x "<< gp1.x()<<" y "<< gp1.y()<<" z "<< gp1.z()<< " eta "<< gp1.eta()<<" phi "<< gp1.phi()<< std::endl;
  gp_st_layer3[1] = GlobalPoint(gp2);
  if (fabs(gp2.z())>100) hasStub_st[1] = true;
  else if (verbose_>0)
      std::cout <<" gp_st2 x "<< gp2.x()<<" y "<< gp2.y()<<" z "<< gp2.z()<< " eta "<< gp2.eta()<<" phi "<< gp2.phi()<< std::endl;
  gp_st_layer3[2] = GlobalPoint(gp3);
  if (fabs(gp3.z())>100) hasStub_st[2] = true;
  else if (verbose_>0)
      std::cout <<" gp_st3 x "<< gp3.x()<<" y "<< gp3.y()<<" z "<< gp3.z()<< " eta "<< gp3.eta()<<" phi "<< gp3.phi()<< std::endl;
  gp_st_layer3[3] = GlobalPoint(gp4);
  if (fabs(gp4.z())>100) hasStub_st[3] = true;

  gp_ge11 = GlobalPoint(gp5);
  if (fabs(gp_ge11.z())>100) hasGEMPad_st1 = true;
  else if (verbose_>0)
      std::cout <<" gp_ge11 x "<< gp_ge11.x()<<" y "<< gp_ge11.y()<<" z "<< gp_ge11.z()<< " eta "<< gp_ge11.eta()<<" phi "<< gp_ge11.phi()<< std::endl;
  gp_ge21 = GlobalPoint(gp6);
  if (fabs(gp_ge21.z())>100) hasGEMPad_st2 = true;
  else if (verbose_>0)
      std::cout <<" gp_ge21 x "<< gp_ge21.x()<<" y "<< gp_ge21.y()<<" z "<< gp_ge21.z()<< " eta "<< gp_ge21.eta()<<" phi "<< gp_ge21.phi()<< std::endl;

  if (hasStub_st[0] and hasStub_st[1])
  	npar = npar_in;
  else
      npar = -1;

  if (hasGEMPad_st1)
  	phi_gem[0] = gp_ge11.phi();
  if (hasGEMPad_st2)
  	phi_gem[1] = gp_ge21.phi();

  if (hasStub_st[1] and fabs(gp2.eta())<1.6 and fabs(gp2.eta())>=1.2 and not(hasGEMPad_st2))
      meRing =2;
  else if (hasStub_st[1] and fabs(gp2.eta())<=2.4 and fabs(gp2.eta())>=1.6)
      meRing =1;

  if (hasStub_st[0] and hasStub_st[1])
	xfactor = (gp_st_layer3[1].perp()/gp_st_layer3[0].perp()-1.0)/fabs(gp_st_layer3[0].z()-gp_st_layer3[1].z());
  //second station
  if (hasStub_st[1]) eta_st2 = gp_st_layer3[1].eta();

}


//DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(){ //test constructor
//}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const L1MuDTTrackSegPhiContainer& stubs,
                                                                   const edm::EventSetup& es,
                                                                   const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  initVariables();

  es_.get<IdealMagneticFieldRecord>().get(magfield_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", propagator_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorOpposite", propagatorOpposite_);
  es_.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",      propagatorAny_);

  // check which stubs are available
  for (auto& stub: stubs){
    if (stub.station()==1){
      has_stub_mb1 = true;
      phi_mb1 = stub.phiValue();
      phib_mb1 = stub.phibValue();
      dphi_mb1 = PtassignmentHelper::normalizePhi(phi_mb1 + phib_mb1);
    }
    if (stub.station()==2){
      has_stub_mb2 = true;
      phi_mb2 = stub.phiValue();
      phib_mb2 = stub.phibValue();
      dphi_mb2 = PtassignmentHelper::normalizePhi(phi_mb2 + phib_mb2);
    }
    if (stub.station()==3){
      has_stub_mb3 = true;
      phi_mb3 = stub.phiValue();
      phib_mb3 = stub.phibValue();
      dphi_mb3 = PtassignmentHelper::normalizePhi(phi_mb3 + phib_mb3);
    }
    if (stub.station()==4){
      has_stub_mb4 = true;
      phi_mb4 = stub.phiValue();
      phib_mb4 = stub.phibValue();
      dphi_mb1 = PtassignmentHelper::normalizePhi(phi_mb4 + phib_mb4);
    }
  }

  // calculate bending angles
  if (has_stub_mb1 and has_stub_mb2) dPhi_barrel_dir_12 = deltaPhi(dphi_mb1, dphi_mb2);
  if (has_stub_mb1 and has_stub_mb3) dPhi_barrel_dir_13 = deltaPhi(dphi_mb1, dphi_mb3);
  if (has_stub_mb1 and has_stub_mb4) dPhi_barrel_dir_14 = deltaPhi(dphi_mb1, dphi_mb4);
  if (has_stub_mb2 and has_stub_mb3) dPhi_barrel_dir_23 = deltaPhi(dphi_mb2, dphi_mb3);
  if (has_stub_mb2 and has_stub_mb4) dPhi_barrel_dir_24 = deltaPhi(dphi_mb2, dphi_mb4);
  if (has_stub_mb3 and has_stub_mb4) dPhi_barrel_dir_34 = deltaPhi(dphi_mb3, dphi_mb4);
}


DisplacedMuonTriggerPtassignment::~DisplacedMuonTriggerPtassignment(){


}

void DisplacedMuonTriggerPtassignment::initVariables()
{
  muScalesCacheID_ = 0ULL ;
  muPtScaleCacheID_ = 0ULL ;

  eta_st2 = -9;
  hasGEMPad_st1 = false;
  hasGEMPad_st2 = false;
  npar = -1;
  meRing = -1;
  //position-based
  ddY123 = -99;
  deltaY12 = -99;
  deltaY23 = -99;

  //direction-based
  phiM_st1 = -9;
  phiM_st2 = -9;
  phiM_st12 = -9;
  phiM_st23 = -9;
  dPhi_dir_st1_st2 = -9;
  dPhi_dir_st1_st12 = -9;
  dPhi_dir_st2_st23 = -9;
  dPhi_dir_st12_st23 = -9;

  position_pt = 0.0;
  direction_pt = 0.0;
  hybrid_pt = 0.0;
  /// barrel
  has_stub_mb1 = false;
  has_stub_mb2 = false;
  has_stub_mb3 = false;
  has_stub_mb4 = false;
  phi_mb1 = -9;
  phi_mb2 = -9;
  phi_mb3 = -9;
  phi_mb4 = -9;
  phib_mb1 = -9;
  phib_mb2 = -9;
  phib_mb3 = -9;
  phib_mb4 = -9;
}


void DisplacedMuonTriggerPtassignment::setupGeometry(const edm::EventSetup& es)
{

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


}


void DisplacedMuonTriggerPtassignment::setupTriggerScales(const edm::EventSetup& es)
{
  hasMuScales_ = true;
  hasMuPtScale_ = true;

  try {
    es.get<L1MuTriggerScalesRcd>().get(muScales_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerScalesRcd>& e) {
    hasMuScales_ = false;
    LogDebug("DisplacedMuonTriggerPtassignment") << "+++ Info: L1MuTriggerScalesRcd is unavailable. +++\n";
  }

  try {
    es.get<L1MuTriggerPtScaleRcd>().get(muPtScale_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerPtScaleRcd>& e) {
    hasMuPtScale_ = false;
    LogDebug("DisplacedMuonTriggerPtassignment") << "+++ Info: L1MuTriggerPtScaleRcd is unavailable. +++\n";
  }
}

//void DisplacedMuonTriggerPtassignment::fitComparatorsLCT(const CSCComparatorDigiCollection& hCSCComparators, const CSCCorrelatedLCTDigi& stub,
//	                                  CSCDetId ch_id, float& fit_phi_layer1, float& fit_phi_layer3, float& fit_phi_layer6,
//					  float& fit_z_layer1, float& fit_z_layer3, float& fit_z_layer6, float& perp)
void DisplacedMuonTriggerPtassignment::fitComparatorsLCT(const CSCComparatorDigiCollection& hCSCComparators, const CSCCorrelatedLCTDigi& stub,
	                                  CSCDetId ch_id, float* fit_phi_layers, float* fit_z_layers, float& perp)
{

  auto cscChamber = cscGeometry_->chamber(ch_id);

  // fetch the CSC comparator digis in this chamber
  CSCComparatorDigiContainerIds compDigisIds;
  for (int iLayer=1; iLayer<=6; ++iLayer){
    CSCDetId layerId(ch_id.endcap(), ch_id.station(), ch_id.ring(), ch_id.chamber(), iLayer);
    // get the digis per layer
    auto compRange = hCSCComparators.get(layerId);
    CSCComparatorDigiContainer compDigis;
    for (auto compDigiItr = compRange.first; compDigiItr != compRange.second; compDigiItr++) {
      auto compDigi = *compDigiItr;
      //if (stub.getTimeBin() < 4 or stub.getTimeBin() > 8) continue;
      int stubHalfStrip(getHalfStrip(compDigi));
      // these comparator digis never fit the pattern anyway!
      if (std::abs(stubHalfStrip-stub.getStrip())>5) continue;
      // check if this comparator digi fits the pattern
      //if(verbose) std::cout << "Comparator digi L1Mu " << layerId << " " << compDigi << " HS " << stubHalfStrip << " stubKeyHS " << stub.getStrip() << std::endl;
      if (comparatorInLCTPattern(stub.getStrip(), stub.getPattern(), iLayer, stubHalfStrip)) {
        //if(verbose) std::cout<<"\tACCEPT"<<std::endl;
        compDigis.push_back(compDigi);
      }
      // else{
      //   if(verbose) std::cout<<"\tDECLINE!"<<std::endl;
      // }
    }
    // if(verbose) if (compDigis.size() > 2) std::cout << ">>> INFO: " << compDigis.size() << " matched comp digis in this layer!" << std::endl;
    compDigisIds.push_back(std::make_pair(layerId, compDigis));
  }

  // get the z and phi positions
  perp = 0.0;
  std::vector<float> phis;
  std::vector<float> zs;
  std::vector<float> ephis;
  std::vector<float> ezs;
  std::vector<float> status;
  for (auto p: compDigisIds){
    auto detId = p.first;
    float phi_tmp = 0.0;
    float perp_tmp = 0.0;
    float z_tmp = 0.0;
    if (p.second.size()==0) continue;
    for (auto hit: p.second){
      float fractional_strip = getFractionalStrip(hit);
      auto layer_geo = cscChamber->layer(detId.layer())->geometry();
      float wire = layer_geo->middleWireOfGroup(stub.getKeyWG() + 1);
      LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
      GlobalPoint csc_gp = cscGeometry_->idToDet(detId)->surface().toGlobal(csc_intersect);
      float gpphi = csc_gp.phi();

      if (phis.size()>0 and gpphi>0 and phis[0]<0 and  (gpphi-phis[0])>3.1416)
        phi_tmp += (gpphi-2*3.1415926);
      else if (phis.size()>0 and gpphi<0 and phis[0]>0 and (gpphi-phis[0])<-3.1416)
        phi_tmp += (gpphi+2*3.1415926);
      else
        phi_tmp += (csc_gp.phi());

      z_tmp = csc_gp.z();
      perp_tmp += csc_gp.perp();
    }
    //in case there are more than one comparator digis in one layer
    perp_tmp = perp_tmp/(p.second).size();
    phi_tmp = phi_tmp/(p.second).size();
    if (verbose_>0)
    	std::cout <<"detid "<< detId <<" perp "<< perp_tmp <<" phi "<< phi_tmp <<" z "<< z_tmp << std::endl;
    perp += perp_tmp;
    phis.push_back(phi_tmp);
    zs.push_back(z_tmp);
    ezs.push_back(0);
    // phis.push_back(csc_gp.phi());
    ephis.push_back(gemvalidation::cscHalfStripWidth(detId)/sqrt(12));
  }


  CSCDetId key_id(ch_id.endcap(), ch_id.station(), ch_id.ring(), ch_id.chamber(), CSCConstants::KEY_CLCT_LAYER);
  float fractional_strip = 0.5 * (stub.getStrip() + 1) - 0.25;
  auto layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();
  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(stub.getKeyWG() + 1);
  LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  GlobalPoint csc_gp = cscGeometry_->idToDet(key_id)->surface().toGlobal(csc_intersect);
  perp = csc_gp.perp();
  // use average perp
  //perp = perp/phis.size();
  // do a fit to the comparator digis
  float alpha = -99., beta = 0.;
  PtassignmentHelper::calculateAlphaBeta(zs, phis, ezs, ephis, status, alpha, beta);
  if (phis.size() <= 2 or fabs(alpha)>=99){
      if (verbose_>0)
      	std::cout <<"warning, falied to fit comparator digis,num of digis: "<< phis.size()<<" alpha "<< alpha <<" beta "<< beta << std::endl;
      alpha = csc_gp.phi();
      beta = 0.0;
  }
  if (verbose_>0)
      std::cout <<"fitting results: alpha "<< alpha <<" beta "<< beta << std::endl;
  for (int i=0; i<6; i++){
      fit_z_layers[i] = cscChamber->layer(i+1)->centerOfStrip(20).z();
      fit_phi_layers[i] = PtassignmentHelper::normalizePhi(alpha + beta * fit_z_layers[i]);
      if (verbose_>0)
      	std::cout <<"i "<< i <<" fit_z "<< fit_z_layers[i]<< " fit_phi "<< fit_phi_layers[i]<<" perp "<< perp << std::endl;
  }

}


void DisplacedMuonTriggerPtassignment::fitTrackRadius(GlobalPoint* gps, float* radius)
{

  std::vector<float> gps_r;
  std::vector<float> gps_z;
  std::vector<float> gps_er;
  std::vector<float> gps_ez;
  std::vector<float> status;
  for (int i=0; i<4; i++){
  	if (not(hasStub_st[i])) continue;
	gps_r.push_back(gps[i].perp());
	gps_z.push_back(gps[i].z());
	gps_er.push_back(1.5);//how to set error on r?
	gps_ez.push_back(0.);

  }
  float alpha = 0., beta = 0.;
  PtassignmentHelper::calculateAlphaBeta(gps_z, gps_r, gps_ez, gps_er, status, alpha, beta);
  for (int i=0; i<4; i++){
  	if (hasStub_st[i])
	    radius[i] = alpha + beta*gps[i].z();
  	else
	    radius[i] = 0.0;
	if (fabs(radius[i]-gps[i].perp())>=.02*gps[i].perp() and hasStub_st[i]){
	    std::cout <<" warning!!! difference bewteen before fitting and after fitting is large "<< std::endl;
	//if (verbose_>=0)
	    std::cout <<"station "<< i+1 <<" z "<< gps[i].z() <<" radius from gp "<< gps[i].perp()<<" from fit "<< radius[i]<< std::endl;
	}
  }


}

void DisplacedMuonTriggerPtassignment::globalPositionOfLCT(const CSCCorrelatedLCTDigi stub, CSCDetId chid)
{
  float perp;
  int st = chid.station();
  fitComparatorsLCT(*hCSCComparators.product(), stub, chid, phi_st_layers[st-1], z_st_layers[st-1], perp);
  //gp calculated here can have negative Z!!! later use fabs() to get distance

  gp_st_layer1[st-1] = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][0], z_st_layers[st-1][0]));
  gp_st_layer3[st-1] = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][2], z_st_layers[st-1][2]));
  gp_st_layer6[st-1] = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][5], z_st_layers[st-1][5]));
  if (verbose_>0)
      std::cout <<"LCT position, chid "<< chid<<" hs "<<stub.getStrip()+1<<" wg "<< stub.getKeyWG()+1 <<" gp eta "<< gp_st_layer3[st-1].eta()<<" phi "<<gp_st_layer3[st-1].phi()<<" perp "<< gp_st_layer3[st-1].perp() << std::endl;

}


void DisplacedMuonTriggerPtassignment::globalPositionOfLCT(CSCCorrelatedLCTDigiContainer stubs, CSCDetId chid)
{
  float dR = 99;
  int st = chid.station();
  GlobalPoint gp_ref;
  bool hasRefStub = false;
  unsigned int beststub=0;
  //find closest stub as ref
  for (int i=st-1; i>0; i--)
      if (hasStub_st[i-1]){
  	gp_ref = GlobalPoint(gp_st_layer3[i-1]);
	hasRefStub = true;
	break;
  	}
  if (hasRefStub and stubs.size()>1){
        unsigned int istub = -1;
  	for (auto stub : stubs){
	    	istub++;
  		globalPositionOfLCT(stub, chid);
		//assign higher weight to phi comparison
    		float dphi = 10.*deltaPhi(gp_st_layer3[st-1].phi(), gp_ref.phi());
    		float deta = gp_st_layer3[st-1].eta() - gp_ref.eta();
    		float curr_dr2 = dphi*dphi + deta*deta;
		if (curr_dr2<dR){
		    dR = curr_dr2;
		    beststub = istub;
		}
  	}
   }
  if (beststub >= stubs.size())
      std::cout <<"error beststub >= stubs.size() , beststub "<< beststub <<" stubs size "<< stubs.size() << std::endl;
  globalPositionOfLCT(stubs[beststub], chid);
  if (verbose_>0)
      std::cout <<"LCT position, chid "<< chid<<" hs "<<stubs[beststub].getStrip()+1<<" wg "<< stubs[beststub].getKeyWG()+1 <<" gp eta "<< gp_st_layer3[st-1].eta()<<" phi "<<gp_st_layer3[st-1].phi()<<" perp "<< gp_st_layer3[st-1].perp() << std::endl;
}


void DisplacedMuonTriggerPtassignment::globalPositionOfGEMPad(const GEMCSCPadDigi gempad, GEMDetId gemid)
{
  GEMDetId ch_id(gemid.region(), gemid.ring(), gemid.station(), gemid.layer(), gemid.chamber(), 0);
  const GEMChamber* gemChamber(gemGeometry_->chamber(ch_id.rawId()));
  auto gemRoll(gemChamber->etaPartition(gemid.roll()));//any roll
  /*ignore this since Geometry of GE21 is wrong(should be 2strip-pad)?
  const int nGEMPads(gemRoll->npads());
  if (gempad.pad() > nGEMPads or gempad.pad() < 0){
      std::cout <<" gempad.pad() is within pad range gempad "<< gempad <<" npad "<< nGEMPads << std::endl;
      return;
  }*/

  if (gemid.station() == 1){
  	const LocalPoint lpGEM(gemRoll->centreOfPad(gempad.pad()));
  	gp_ge11 = GlobalPoint(gemRoll->toGlobal(lpGEM));
	phi_gem[0] = gp_ge11.phi();
	if (verbose_>0) std::cout <<" gempad in GE11 id " << gemid <<" gp eta "<< gp_ge11.eta()<<" phi "<< gp_ge11.phi()<<" pad "<<gempad.pad()<< std::endl;
  }else if (gemid.station() == 3){
  	const LocalPoint lpGEM(gemRoll->centreOfStrip(float(gempad.pad()*2.0-1.0)));
  	gp_ge21 = GlobalPoint(gemRoll->toGlobal(lpGEM));
	phi_gem[1] = gp_ge21.phi();
	if (verbose_>0) std::cout <<" gempad in GE21 id "<< gemid <<" gp eta "<< gp_ge21.eta()<<" phi "<< gp_ge21.phi()<<" pad "<<gempad.pad()<< std::endl;
  }else if (verbose_>0)
      std::cout <<" gemid "<< gemid  <<" not in station 1 or 3" << std::endl;

}


void DisplacedMuonTriggerPtassignment::globalPositionOfGEMPad(GEMCSCPadDigiContainer gempads, GEMDetId gemid)
{
  GEMDetId ch_id(gemid.region(), gemid.ring(), gemid.station(), gemid.layer(), gemid.chamber(), 0);
  const GEMChamber* gemChamber(gemGeometry_->chamber(ch_id.rawId()));
  auto gemRoll(gemChamber->etaPartition(gemid.roll()));//any roll
  const int nGEMPads(gemRoll->npads());
  int st = gemid.station();
  if (st==3) st=2;//use CSC station as reference
  for (auto gempad : gempads){
  	if (gempad.pad() > nGEMPads or gempad.pad() < 0){
      		std::cout <<" gempad.pad() is within pad range gempad "<< gempad <<" npad "<< nGEMPads << std::endl;
      		return;
	}
  	const LocalPoint lpGEM(gemRoll->centreOfPad(gempad.pad()));
  	GlobalPoint gp_pad = GlobalPoint(gemRoll->toGlobal(lpGEM));
	if (hasStub_st[st-1] and fabs(deltaPhi(gp_pad.phi(), gp_st_layer3[st-1].phi()))<fabs(dphi_gemcsc_st[st-1]) ){
		gp_ge11 = GlobalPoint(gp_pad);
		dphi_gemcsc_st[st-1] = fabs(deltaPhi(gp_pad.phi(), gp_st_layer3[st-1].phi()));
	}

  }

}

float DisplacedMuonTriggerPtassignment::deltaYcalculation(GlobalPoint gp1, GlobalPoint gp2) const
{
   float anglea = gp2.phi();
   float newyst1 = -gp1.x()*sin(anglea) + gp1.y()*cos(anglea);
   float newyst2 = -gp2.x()*sin(anglea) + gp2.y()*cos(anglea);
   return (newyst2-newyst1);

}


float DisplacedMuonTriggerPtassignment::deltadeltaYcalculation(GlobalPoint gp1, GlobalPoint gp2, GlobalPoint gp3, float eta, int par) const
{

   float anglea = gp2.phi();
   float newyst1 = -gp1.x()*sin(anglea) + gp1.y()*cos(anglea);
   float newyst2 = -gp2.x()*sin(anglea) + gp2.y()*cos(anglea);
	//float newxst3 = gp3.x()*cos(anglea) + gp3.y()*sin(anglea);
   float newyst3 = -gp3.x()*sin(anglea) + gp3.y()*cos(anglea);
   float deltay12 = newyst2-newyst1;
   float deltay23 = newyst3-newyst2;
   //std::cout <<" angle in st2 "<< anglea <<" newyst1 "<< newyst1 <<" newyst2 "<< newyst2 << " newyst3 "<< newyst3 << std::endl;
   int neta = PtassignmentHelper::GetEtaPartition(eta);

   if (par<0 or par>3 or neta==-1) return -99;
   return (deltay23-PtassignmentHelper::PositionEpLUT[par][neta][0]*deltay12);

}

float DisplacedMuonTriggerPtassignment::phiMomentum_Xfactor(float phi_CSC, float phi_GEM, float xfactor) const
{


   if (fabs(phi_CSC) > M_PI or fabs(phi_GEM) > M_PI) return -9;
   float dphi = deltaPhi(phi_CSC,phi_GEM);
   float y = 1.0-cos(dphi)- xfactor;

   float phi_diff = 0.0;
   if (fabs(y) > 0.0) phi_diff = atan(sin(dphi)/y);
   else phi_diff = M_PI/2.0;

   if (phi_diff <= -M_PI) phi_diff = phi_diff+2*M_PI;
   else if (phi_diff > M_PI) phi_diff = phi_diff-2*M_PI;

   float phiM = phi_GEM-phi_diff;
   if (phiM <= -M_PI) phiM = phiM+2*M_PI;
   else if (phiM > M_PI) phiM = phiM-2*M_PI;

   //std::cout <<"PhiMomentum_Xfactor: dphi "<< dphi <<" phi_poistion1 "<< phi_GEM <<" phi_position2 "<< phi_CSC <<" Xfactor "<<X <<" phi_diff "<< phi_diff <<" phiM "<< phiM << std::endl;

   return phiM;
}


bool DisplacedMuonTriggerPtassignment::runPositionbased()
{
   if (npar<0 or npar>=4 or not(hasStub_st[2])){
        std::cout <<" failed to runPositionbased  npar "<< npar << std::endl;
   	return false;
   }
   ddY123 = deltadeltaYcalculation(gp_st_layer3[0], gp_st_layer3[1], gp_st_layer3[2], gp_st_layer3[1].eta(), npar);
   deltaY12 = deltaYcalculation(gp_st_layer3[0], gp_st_layer3[1]);
   deltaY23 = -deltaYcalculation(gp_st_layer3[2], gp_st_layer3[1]);
   if (npar>=0 and npar<=3){
        position_pt = 2.0;
   	int neta = PtassignmentHelper::GetEtaPartition(eta_st2);
   	for (int i=0; i<PtassignmentHelper::NPt2; i++){
	    if (fabs(ddY123) <= PtassignmentHelper::PositionbasedDDYLUT[i][neta][npar])
		position_pt = float(PtassignmentHelper::PtBins2[i]);
	    else
		break;
	    if (verbose_>0)
		std::cout <<"eta "<< eta_st2 <<" neta "<< neta <<" npar "<< npar <<" fabs ddY123 "<< fabs(ddY123) <<" cut "<< PtassignmentHelper::PositionbasedDDYLUT[i][neta][npar] <<" position pt "<< position_pt<<std::endl;
	}
   }
   return true;
}

// run the direction based algorithm. Option to include GE21 hits or not
bool DisplacedMuonTriggerPtassignment::runDirectionbased(bool useGE21)
{
  if (useGE21 and meRing==1) return runDirectionbasedGE21();
  else return runDirectionbasedCSConly();
}

//use GE21 if GE21 pads are available. use GE11 if GE11 pads are available
bool DisplacedMuonTriggerPtassignment::runDirectionbasedGE21()
{
   //if (not (npar<4 and npar>=0 and hasGEMPad_st1 and hasGEMPad_st2)) return false;
   if (not (npar<4 and npar>=0)) return false;
   //if (fabs(phi_gem[1])>4) return false;//check this because we want to use setPhiGE21() to set phi_gem[1](using 2strips-pad)

   float xfactor_st1 = 0.0;
   float xfactor_st2 = 0.0;
   if (meRing==1 and hasGEMPad_st1){
	xfactor_st1 = xfactor*fabs(gp_ge11.z() - gp_st_layer3[0].z());
   	phiM_st1 = phiMomentum_Xfactor(gp_st_layer3[0].phi(), gp_ge11.phi(), xfactor_st1);//
   }else{
        xfactor_st1 = xfactor*fabs(z_st_layers[0][0] - z_st_layers[0][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[0][5])+1);
   	phiM_st1 = phiMomentum_Xfactor(gp_st_layer6[0].phi(), gp_st_layer1[0].phi(), xfactor_st1);//
   }

   if (meRing==1 and hasGEMPad_st2){
	xfactor_st2 = xfactor*fabs(gp_ge21.z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   	phiM_st2 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), phi_gem[1], xfactor_st2);
   }else if(meRing==2){
   	xfactor_st2 = xfactor*fabs(z_st_layers[1][0] - z_st_layers[1][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[1][5])+1);
   	phiM_st2 = phiMomentum_Xfactor(gp_st_layer6[1].phi(), gp_st_layer1[1].phi(), xfactor_st2);//
   }else return false;

   float xfactor_st12 = xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   float xfactor_st23 = xfactor*fabs(gp_st_layer3[1].z() - gp_st_layer3[2].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[2].z())+1);
   phiM_st12 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), gp_st_layer3[0].phi(), xfactor_st12);
   phiM_st23 = phiMomentum_Xfactor(gp_st_layer3[2].phi(), gp_st_layer3[1].phi(), xfactor_st23);
   if (verbose_>0)  std::cout <<"DisplacedMuonTrigger, direction with GE21, meRing "<< meRing <<" xfactor_st1 "<< xfactor_st1 <<" phiM_st1 "<< phiM_st1
       			<<" xfactor_st2 "<< xfactor_st2 <<" phiM_st2 "<< phiM_st2 << std::endl;

   //make sure both phiM_st1 and phiM_st2 are reasonable, 4 can be changed into M_PI later
   dPhi_dir_st1_st2 = (fabs(phiM_st1)<4 and fabs(phiM_st2)<4)? deltaPhi(phiM_st1, phiM_st2):-9;
   dPhi_dir_st1_st12 = (fabs(phiM_st1)<4 and fabs(phiM_st12)<4)? deltaPhi(phiM_st1, phiM_st12):-9;
   dPhi_dir_st2_st23 = (fabs(phiM_st2)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st2, phiM_st23):-9;
   dPhi_dir_st12_st23 = (fabs(phiM_st12)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st12, phiM_st23):-9;

   if (npar>=0 and npar<=3){
        direction_pt = 2.0;
   	int neta = PtassignmentHelper::GetEtaPartition(eta_st2);
   	for (int i=0; i<PtassignmentHelper::NPt2; i++){
	    if (fabs(dPhi_dir_st1_st2) <= PtassignmentHelper::DirectionbasedDeltaPhiLUT[i][neta][npar])
		direction_pt = float(PtassignmentHelper::PtBins2[i]);
	    else
		break;
	    if (verbose_>0)
		std::cout <<"eta "<< eta_st2 <<" neta "<< neta <<" npar "<< npar <<" fabs dphi "<< fabs(dPhi_dir_st1_st2) <<" cut "<< PtassignmentHelper::DirectionbasedDeltaPhiLUT[i][neta][npar] <<" direction pt "<< direction_pt<<std::endl;
	}
   }
   return true;
}


bool DisplacedMuonTriggerPtassignment::runDirectionbasedCSConly()
{

   //z_st_layers should be used at sim level, set Z and phi for layer1 and layer6 at sim level, or rebuild constructor?
   if (not ( npar<4 and npar>=0)) return false;
   float xfactor_st1 = 0.0;
   float xfactor_st2 = 0.0;
   if (meRing==1 and hasGEMPad_st1){
   	xfactor_st1 = xfactor*fabs(gp_ge11.z() - gp_st_layer3[0].z());
   	phiM_st1 = phiMomentum_Xfactor(gp_st_layer3[0].phi(), gp_ge11.phi(), xfactor_st1);//
   }else{
   	//xfactor_st1 = xfactor*fabs(gp_st1_layer1.z() - gp_st1_layer6.z())/(xfactor*fabs(gp_st1.z() - gp_st1_layer6.z())+1);
   	xfactor_st1 = xfactor*fabs(z_st_layers[0][0] - z_st_layers[0][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[0][5])+1);
   	phiM_st1 = phiMomentum_Xfactor(phi_st_layers[0][5], phi_st_layers[0][0], xfactor_st1);//
   }
   xfactor_st2 = xfactor*fabs(z_st_layers[1][0] - z_st_layers[1][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[1][5])+1);
   float xfactor_st12 = xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   float xfactor_st23 = xfactor*fabs(gp_st_layer3[1].z() - gp_st_layer3[2].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[2].z())+1);
   phiM_st2 = phiMomentum_Xfactor(phi_st_layers[1][5], phi_st_layers[1][0], xfactor_st2);
   phiM_st12 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), gp_st_layer3[0].phi(), xfactor_st12);
   phiM_st23 = phiMomentum_Xfactor(gp_st_layer3[2].phi(), gp_st_layer3[1].phi(), xfactor_st23);
   //if phi in layer1 and layer6 in station1 and 2 are not set, then here phiM return -9
   if (verbose_>0)  std::cout <<"DisplacedMuonTrigger CSConly direction, meRing "<< meRing <<" xfactor_st1 "<< xfactor_st1 <<" phiM_st1 "<< phiM_st1
       			<<" xfactor_st2 "<< xfactor_st2 <<" phiM_st2 "<< phiM_st2 << std::endl;

   dPhi_dir_st1_st2 = (fabs(phiM_st1)<4 and fabs(phiM_st2)<4)? deltaPhi(phiM_st1, phiM_st2):-9;
   dPhi_dir_st1_st12 = (fabs(phiM_st1)<4 and fabs(phiM_st12)<4)? deltaPhi(phiM_st1, phiM_st12):-9;
   dPhi_dir_st2_st23 = (fabs(phiM_st2)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st2, phiM_st23):-9;
   dPhi_dir_st12_st23 = (fabs(phiM_st12)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st12, phiM_st23):-9;

   if (npar>=0 and npar<=3){
        direction_pt = 2.0;
   	int neta = PtassignmentHelper::GetEtaPartition(eta_st2);
   	for (int i=0; i<PtassignmentHelper::NPt2; i++){
	    if (fabs(dPhi_dir_st1_st2) <= PtassignmentHelper::DirectionbasedDeltaPhiME21CSConlyLUT[i][neta][npar])
		direction_pt = float(PtassignmentHelper::PtBins2[i]);
	    else
		break;
	    if (verbose_>0)
		std::cout <<"eta "<< eta_st2 <<" neta "<< neta <<" npar "<< npar <<" fabs dphi "<< fabs(dPhi_dir_st1_st2) <<" cut "<< PtassignmentHelper::DirectionbasedDeltaPhiME21CSConlyLUT[i][neta][npar] <<" direction pt "<< direction_pt<<std::endl;
	}

   }
   return true;
}


bool DisplacedMuonTriggerPtassignment::runHybrid(bool useGE21)
{

   //firstly to run through position-based and direction-based
   bool checkPosition = runPositionbased();
   bool checkDirection = false;
   if (useGE21)
   	checkDirection = runDirectionbasedGE21();
   else
 	checkDirection = runDirectionbasedCSConly();
   if (not (checkPosition and checkDirection))
       return false;
   hybrid_pt = 2.0;
   if (npar>=0 and npar<=3){
   	int neta = PtassignmentHelper::GetEtaPartition(eta_st2);
	if (fabs(ddY123)>=40 or fabs(dPhi_dir_st1_st2)>=1.0){//rejected by hybrid
	    return true;
	}
	//ignore pt=40
   	for (int i=0; i<PtassignmentHelper::NPt-1; i++){
           if(useGE21 and PtassignmentHelper::ellipse(PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][0],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][1],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][2],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][3],
					//PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][4], ddY123*charge, dPhi_dir_st1_st2*charge) <=1)
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUT[i][neta][npar][4], ddY123, dPhi_dir_st1_st2) <=1)
		hybrid_pt = PtassignmentHelper::PtBins[i];
	   else if(not(useGE21) and PtassignmentHelper::ellipse(PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][0],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][1],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][2],
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][3],
		  		//PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][4], ddY123*charge, dPhi_dir_st1_st2*charge) <=1)
		  			  PtassignmentHelper::HybridDDYAndDeltaPhiLUTME21CSConly[i][neta][npar][4], ddY123, dPhi_dir_st1_st2) <=1)
		hybrid_pt = PtassignmentHelper::PtBins[i];
	   else//make sure LUT is consitent
	   	break;
	   if (verbose_>0)
   		std::cout <<"eta_st2 "<< eta_st2 <<" npar "<< npar <<" charge "<< charge <<" ddY123 "<< ddY123 << " dphi_dir "<< dPhi_dir_st1_st2 <<" hybrid_pt "<< hybrid_pt << std::endl;

	}
   }
   return true;
}

void DisplacedMuonTriggerPtassignment::runDirectionBasedBarrel()
{
  // check case
  int dt_stub_case = getBarrelStubCase(has_stub_mb1, has_stub_mb2, has_stub_mb3, has_stub_mb4);

  barrel_direction_pt = -1;

  switch(dt_stub_case){
  case 0:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_12, "DT1_DT2");
    break;
  case 1:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_13, "DT1_DT3");
    break;
  case 2:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_14, "DT1_DT4");
   break;
  case 3:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_23, "DT2_DT3");
    break;
  case 4:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_24, "DT2_DT4");
    break;
  case 5:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(dPhi_barrel_dir_34, "DT3_DT4");
    break;
  case 6:
    // first dphi is x value, second dphi is y value!!!!
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(dPhi_barrel_dir_12, dPhi_barrel_dir_13, "DT1_DT2__DT1_DT3");
    break;
  case 7:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(dPhi_barrel_dir_12, dPhi_barrel_dir_14, "DT1_DT2__DT1_DT4");
    break;
  case 8:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(dPhi_barrel_dir_13, dPhi_barrel_dir_14, "DT1_DT3__DT1_DT4");
    break;
  case 9:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(dPhi_barrel_dir_23, dPhi_barrel_dir_24, "DT2_DT3__DT2_DT4");
    break;
  case 10:
    barrel_direction_pt = BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(dPhi_barrel_dir_14, dPhi_barrel_dir_23, "DT1_DT4__DT2_DT3");
    break;
  default:
    // all else fails; assign lowest possible pT
    barrel_direction_pt = 2;
    break;
  };
}

int DisplacedMuonTriggerPtassignment::getBarrelStubCase(bool MB1, bool MB2, bool MB3, bool MB4)
{
  if (    MB1 and     MB2 and not MB3 and not MB4) return 0;
  if (    MB1 and not MB2 and     MB3 and not MB4) return 1;
  if (    MB1 and not MB2 and not MB3 and     MB4) return 2;
  if (not MB1 and     MB2 and     MB3 and not MB4) return 3;
  if (not MB1 and     MB2 and not MB3 and     MB4) return 4;
  if (not MB1 and not MB2 and     MB3 and     MB4) return 5;

  if (    MB1 and     MB2 and     MB3 and not MB4) return 6;
  if (    MB1 and     MB2 and not MB3 and     MB4) return 7;
  if (    MB1 and not MB2 and     MB3 and     MB4) return 8;
  if (not MB1 and     MB2 and     MB3 and     MB4) return 9;

  if (MB1 and MB2 and MB3 and MB4) return 10;

  return -1;
}


void DisplacedMuonTriggerPtassignment::runPositionBasedBarrel(){}
void DisplacedMuonTriggerPtassignment::runHybridBasedBarrel(){}


float DisplacedMuonTriggerPtassignment::getlocalPhiDirection(int st) const
{
    //st =1 :station1 , st=2: station2
    //st = 12 : between station1 and station2; st = 23 : between station2 and station3
   if (st==1 and hasStub_st[0]) return phiM_st1;
   else if (st==2 and hasStub_st[1]) return phiM_st2;
   else if (st == 12 and hasStub_st[0] and hasStub_st[1]) return phiM_st12;
   else if (st == 23 and hasStub_st[1] and hasStub_st[2]) return phiM_st23;
   else{
   	std::cout <<" error in getlocalPhiDirection, st  "<<st <<" not in range or not not have stub or GEMpad" << std::endl;
	return -99;
   }
}

float DisplacedMuonTriggerPtassignment::getdeltaPhiDirection(int st1, int st2) const
{
   if (((st1 == 1 and st2 == 2) or (st1 == 1 and st2 == 2)) and hasStub_st[0] and hasStub_st[1]) return dPhi_dir_st1_st2;
   else if (((st1 == 1 and st2 == 12) or (st1 == 12 and st2 == 1)) and hasStub_st[0] and hasStub_st[1]) return dPhi_dir_st1_st12;
   else if (((st1 == 2 and st2 == 23) or (st1 == 23 and st2 == 2)) and hasStub_st[1] and hasStub_st[2]) return dPhi_dir_st2_st23;
   else if (((st1 == 12 and st2 == 23) or (st1 == 23 and st2 == 12)) and hasStub_st[0] and hasStub_st[1] and hasStub_st[2]) return dPhi_dir_st12_st23;
   else{
   	std::cout <<" error in getdeltaPhiDirection, st1 "<< st1 <<" st2 "<< st2 <<" not in range or not not have stub or GEMpad" << std::endl;
	return -99;
   }

}

int
DisplacedMuonTriggerPtassignment::getHalfStrip(const CSCComparatorDigi& digi)
{
  return (digi.getStrip() - 1) * 2 + digi.getComparator();
}

float
DisplacedMuonTriggerPtassignment::getFractionalStrip(const CSCComparatorDigi&d)
{
  return d.getStrip() + d.getComparator()/2. - 3/4.;
}

bool
DisplacedMuonTriggerPtassignment::stubInDTTFTracks(const L1MuDTTrackSegPhi& candidateStub,
                                                   const L1MuDTTrackCollection& l1Tracks) const
{
  bool isMatched = false;
  for (auto tftrack: l1Tracks){
    auto stubCollection = tftrack.second;
    for (auto stub: stubCollection) {
      if (candidateStub == stub) {
        isMatched = true;
        break;
      }
    }
  }
  return isMatched;
}

bool
DisplacedMuonTriggerPtassignment::stubInCSCTFTracks(const CSCCorrelatedLCTDigi& candidateStub,
                                                    const L1CSCTrackCollection& l1Tracks) const
{
  bool isMatched = false;
  for (auto tftrack: l1Tracks){
    auto stubCollection = tftrack.second;
    for (auto detUnitIt = stubCollection.begin(); detUnitIt != stubCollection.end(); detUnitIt++) {
      const auto range = (*detUnitIt).second;
      for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
        //if (!(*digiIt).isValid()) continue;
        auto stub = *digiIt;
        if (candidateStub == stub) {
          isMatched = true;
          break;
        }
      }
    }
  }
  return isMatched;
}

double
phiL1DTTrack(const L1MuDTTrack& track)
{
  int phi_local = track.phi_packed(); //range: 0 < phi_local < 31
  if ( phi_local > 15 ) phi_local -= 32; //range: -16 < phi_local < 15
  double dttrk_phi_global = normalizedPhi((phi_local*(M_PI/72.))+((M_PI/6.)*track.spid().sector()));// + 12*i->scNum(); //range: -16 < phi_global < 147
  // if(dttrk_phi_global < 0) dttrk_phi_global+=2*M_PI; //range: 0 < phi_global < 147
  // if(dttrk_phi_global > 2*M_PI) dttrk_phi_global-=2*M_PI; //range: 0 < phi_global < 143
  return dttrk_phi_global;
}

double
phiL1CSCTrack(const csc::L1Track& track)
{
  unsigned gbl_phi(track.localPhi() + ((track.sector() - 1)*24) + 6);
  if(gbl_phi > 143) gbl_phi -= 143;
  double phi_packed = gbl_phi & 0xff;
  return phi_packed;
}

void DisplacedMuonTriggerPtassignment::calculateTTIsolation()
{
  for (unsigned int j=0; j<tttracks_.size(); ++j) {
    auto l1Tk = tttracks_[j];
    const double l1Tk_pt = l1Tk.getMomentum().perp();
    // const double l1Tk_eta = l1Tk.getMomentum().eta();
    // const double l1Tk_phi = normalizedPhi(l1Tk.getMomentum().phi());
    // const double l1Tk_charge = l1Tk.getRInv()>0? 1: -1;

    // if(false) {
    //   cout << "l1Tk " << j << endl;
    //   cout << "l1Tk_pt " << l1Tk_pt << endl;
    //   cout << "l1Tk_eta " << l1Tk_eta << endl;
    //   cout << "l1Tk_phi " << l1Tk_phi << endl;
    //   cout << "l1Tk_phi_corr " << l1Tk_phi_corr << endl;
    //   cout << "l1Tk_charge " << l1Tk_charge << endl;
    // }

    double l1Tk_eta_prop = -99;
    double l1Tk_phi_prop = -99;
    GlobalPoint ex_point(extrapolateGP(l1Tk));
    if (!(ex_point == GlobalPoint())) {
      l1Tk_eta_prop = ex_point.eta();
      l1Tk_phi_prop = ex_point.phi();
      if(false) {
        cout << "l1Tk_eta_prop " << l1Tk_eta_prop << endl;
        cout << "l1Tk_phi_prop " << l1Tk_phi_prop << endl;
      }
      const double dR_l1Mu_l1Tk_prop = reco::deltaR(l1Tk_eta_prop, l1Tk_phi_prop,
                                                    getTrackEta(), getTrackPhi(2));
      if (dR_l1Mu_l1Tk_prop < L1Mu_L1Tk_dR_min_) {
        L1Mu_L1Tk_dR_min_ = dR_l1Mu_l1Tk_prop;
        L1Mu_L1Tk_pt_min_ = l1Tk_pt;
      }
    }
  }
  // end of loop on TTTracks

  if (L1Mu_L1Tk_dR_min_ <= 0.12 and L1Mu_L1Tk_pt_min_ >= 4) isLooseVeto_ = true;
  if (L1Mu_L1Tk_dR_min_ <= 0.12 and L1Mu_L1Tk_pt_min_ >= 3) isMediumVeto_ = true;
  if (L1Mu_L1Tk_dR_min_ <= 0.12 and L1Mu_L1Tk_pt_min_ >= 2) isTightVeto_ = true;
}

GlobalPoint
DisplacedMuonTriggerPtassignment::extrapolateGP(const TTTrack< Ref_PixelDigi_ > &tk, int station)
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
DisplacedMuonTriggerPtassignment::propagateToZ(const GlobalPoint &inner_point, const GlobalVector &inner_vec, double charge, double z) const
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
DisplacedMuonTriggerPtassignment::propagateToR(const GlobalPoint &inner_point, const GlobalVector &inner_vec, double charge, double R) const
{
  Cylinder::CylinderPointer my_cyl(Cylinder::build(Surface::PositionType(0,0,0), Surface::RotationType(), R));

  FreeTrajectoryState state_start(inner_point, inner_vec, charge, &*magfield_);

  TrajectoryStateOnSurface tsos(propagator_->propagate(state_start, *my_cyl));
  if (!tsos.isValid()) tsos = propagatorOpposite_->propagate(state_start, *my_cyl);
  return tsos;
}

#endif

