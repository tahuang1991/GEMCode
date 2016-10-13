#ifndef GEMCode_GEMValidation_DisplacedMuonTriggerPtassignment_cc
#define GEMCode_GEMValidation_DisplacedMuonTriggerPtassignment_cc

/**\class DisplacedMuonTriggerPtassignment

  Displaced Muon Trigger Design based on Muon system
  
  Author: tao.huang@cern.ch, sven.dildick@cern.ch

*/
#include "DataFormats/Math/interface/deltaPhi.h"
#include "GEMCode/GEMValidation/interface/DisplacedMuonTriggerPtassignment.h"
#include "GEMCode/GEMValidation/interface/PtassignmentHelper.h"


//endcap, we need LCTs and associated cscid, gempads and associated gemid, and all gemometry
//to get position from fitting, we also need all comparator digis
//step0 get LCTs and associated cscids, GEMPads and associated gemids, and geometry. 
//step1 get fitting positions from fitting compara digis after assoicating comparator digis to LCTs
//step2 calculate all variables used pt assignment, requires eta,phi,radius,Z
//step3 assgin L1 pt according to LUTs (in short future)
DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const CSCCorrelatedLCTDigiCollection* lcts, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{  
  setupGeometry(es);
}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const CSCCorrelatedLCTDigiContainer lcts, const CSCDetIdContainer cscids, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
}

//chamberid_lcts: LCTs matched to simmuon and their associated chamberid, detid_pads: gempads matched to simmuon and their associated detid_pads
DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(std::map<unsigned int, CSCCorrelatedLCTDigiContainer> chamberid_lcts, std::map<unsigned int, GEMCSCPadDigiContainer> detid_pads, const edm::EventSetup& es, const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  chamberid_lcts_ = chamberid_lcts;
  detid_pads_ = detid_pads;
  ev.getByLabel("simMuonCSCDigis", "MuonCSCComparatorDigi", hCSCComparators);
  
  initVariables(); 
  for (auto idlcts : chamberid_lcts_){
  	CSCDetId chid(idlcts.first);
/*<<<<<<< Updated upstream
    //check the ring number that muon is flying through, 2nd station as reference
    //later use this one to check whether we should use GE21 or ME21only
    if (chid.station()==2) meRing = chid.ring();
    if (chid.chamber()%2 == 0) isEven[chid.station()-1] = true;
    if (chid.station() == 1 and idlcts.second.size()>0 ) hasStub_st1 = true;
    else if (chid.station() == 2 and idlcts.second.size()>0 ) hasStub_st2 = true;
    else if (chid.station() == 3 and idlcts.second.size()>0 ) hasStub_st3 = true;
    else if (chid.station() == 4 and idlcts.second.size()>0 ) hasStub_st4 = true;
    else {
=======*/
	//check the ring number that muon is flying through, 2nd station as reference
	//later use this one to check whether we should use GE21 or ME21only
	if (chid.station()==2) meRing = chid.ring();
	if (chid.chamber()%2 == 0) isEven[chid.station()-1] = true;
	if (idlcts.second.size()>0) hasStub_st[chid.station()-1] =true;
	/*
	if (chid.station() == 1 and idlcts.second.size()>0 ) hasStub_st1 = true;
	else if (chid.station() == 2 and idlcts.second.size()>0 ) hasStub_st2 = true;
	else if (chid.station() == 3 and idlcts.second.size()>0 ) hasStub_st3 = true;
	else if (chid.station() == 4 and idlcts.second.size()>0 ) hasStub_st4 = true;
	*/
	else {
	    std::cout <<" chid "<< chid <<"  number of lcts "<< idlcts.second.size() << std::endl;
	    continue;
    }
    globalPositionOfLCT(idlcts.second[0], chid);
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
          else if (verbose_>=0)
            std::cout <<" gemid "<< gemid <<" CSC chamber id "<< chid << std::endl;
          //maybe also check the dR(csc, gem)
          globalPositionOfGEMPad(idgempads.second[0], gemid);
        }
      }
    }
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
}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const L1CSCTrack& tftrack,
                                                                   const CSCCorrelatedLCTDigiCollection& lcts,
                                                                   bool doStubRecovery,
                                                                   const edm::EventSetup& es, 
                                                                   const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  initVariables(); 
  
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
    /*
    if (ch_id.station() == 1) hasStub_st1 = true;
    if (ch_id.station() == 2) hasStub_st2 = true;
    if (ch_id.station() == 3) hasStub_st3 = true;
    if (ch_id.station() == 4) hasStub_st4 = true;
    */

    chamberid_lct[ch_id.rawId()] = v;
  }

  // TF eta

  // second step: stub recovery
  const bool atLeast1StubMissing( (not hasStub_st[0]) or
                                  (not hasStub_st[0]) or
                                  (not hasStub_st[0]) or
                                  (not hasStub_st[0]) ); 
  
  if (doStubRecovery and atLeast1StubMissing){
    //int triggerSector = tftrack.sector();

    // for (int endcap=1; endcap<=2; endcap++){
    //   //do not consider stubs in the wrong endcap
    //   int zendcap(endcap!=1 ? -1 : +1 );
    //   if (zendcap * event_.CSCTF_eta[j] < 0) continue;
    //   for (int station=1; station<=4; station++){
          
    //     // ignore station where a L1Mu stub is present!
    //     if (not stubMissingSt1 and station==1) continue;
    //     if (not stubMissingSt2 and station==2) continue;
    //     if (not stubMissingSt3 and station==3) continue;
    //     if (not stubMissingSt4 and station==4) continue;
    //     if(verbose) std::cout << "Recovered stubs in station: " << station << std::endl;
    //     // temp storage of candidate stubs per station and ring
    //     CSCCorrelatedLCTDigiId bestMatchingStub;
    //     int iStub = 0;
    //     for (int ring=1; ring<=3; ring++){
    //       if (station!=1 and ring==3) continue;
    //       //std::cout << "Analyzing ring " << ring << std::endl;
            
    //       for (int chamber=1; chamber<=36; chamber++){
    //         // do not consider invalid detids
    //         if ( (station==2 or station==3 or station==4) and 
    //              (ring==1) and chamber>18) continue;
    //         //std::cout << "Analyzing chamber " << chamber << std::endl; 
    //         // create the detid
    //         CSCDetId ch_id(endcap, station, ring, chamber);
    //         //std::cout << "ch_id " <<  ch_id << std::endl;
    //         // get the stubs in this detid
    //         auto range = CSCCorrelatedLCTs.get(ch_id);
    //         for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
    //           iStub++; 

    //           auto stub(*digiItr);

    //           // check that this stub is not already part of a CSC TF track
    //           if (stubInCSCTFTracks(stub, l1Tracks)) continue;

    //           // trigger sector must be the same
    //           if (triggerSector != ch_id.triggerSector()) continue;
                
    //           // BXs have to match
    //           int deltaBX = std::abs(stub.getBX() - (6 + event_.CSCTF_bx[j]));
    //           if (deltaBX > 1) continue;

    //           if(verbose) std::cout << ch_id << std::endl;
    //           if(verbose) std::cout<<"Candidate " << stub << std::endl;
    //           bestMatchingStub = pickBestMatchingStub(allxs[ch_id.station()-1], allys[ch_id.station()-1], 
    //                                                   bestMatchingStub, std::make_pair(ch_id, stub), 6 + event_.CSCTF_bx[j]);
    //         }
    //         // consider the case ME1a
    //         if (station==1 and ring==1){
    //           CSCDetId me1a_id(endcap, station, 4, chamber);
    //           auto range = CSCCorrelatedLCTs.get(me1a_id);
    //           for (auto digiItr = range.first; digiItr != range.second; ++digiItr){
    //             iStub++;
    //             auto stub(*digiItr);

    //             // check that this stub is not already part of a CSC TF track
    //             if (stubInCSCTFTracks(stub, l1Tracks)) continue;
                  
    //             // trigger sector must be the same
    //             if (triggerSector != me1a_id.triggerSector()) continue;
                  
    //             // BXs have to match
    //             int deltaBX = std::abs(stub.getBX() - (6 + event_.CSCTF_bx[j]));
    //             if (deltaBX > 1) continue; 

    //             if(verbose) std::cout << me1a_id << std::endl;
    //             if(verbose) std::cout<<"Candidate " << stub << std::endl;
    //             bestMatchingStub = pickBestMatchingStub(allxs[me1a_id.station()-1], allys[me1a_id.station()-1], 
    //                                                     bestMatchingStub, std::make_pair(me1a_id, stub), 6 + event_.CSCTF_bx[j]);
                  
    //           }
    //         }
    //       }
    //     }

  }
}


DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(GlobalPoint gp1, GlobalPoint gp2, GlobalPoint gp3, GlobalPoint gp4, GlobalPoint gp5, GlobalPoint gp6, int npar_in, const edm::EventSetup& es, const edm::Event& ev)
: ev_(ev), es_(es), verbose_(0)
{ //sim level

  setupGeometry(es);
  initVariables(); 
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
  if (hasGEMPad_st2)
  	phi_ge21 = gp_ge21.phi();

  if (hasStub_st[1] and fabs(gp2.eta())<1.6 and fabs(gp2.eta())>=1.2 and not(hasGEMPad_st2))
      meRing =2;
  else if (hasStub_st[1] and fabs(gp2.eta())<=2.4 and fabs(gp2.eta())>=1.6)
      meRing =1;
        
  if (hasStub_st[0] and hasStub_st[1])
	xfactor = (gp_st_layer3[1].perp()/gp_st_layer3[0].perp()-1.0)/fabs(gp_st_layer3[0].z()-gp_st_layer3[1].z());
        
}


//DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(){ //test constructor
//}

DisplacedMuonTriggerPtassignment::DisplacedMuonTriggerPtassignment(const L1MuDTTrackSegPhiContainer& tracks,
                                                                   const edm::EventSetup& es, 
                                                                   const edm::Event& ev)
  : ev_(ev), es_(es), verbose_(0)
{
  setupGeometry(es);
  initVariables();
}


DisplacedMuonTriggerPtassignment::~DisplacedMuonTriggerPtassignment(){


}

void DisplacedMuonTriggerPtassignment::initVariables()
{
/*<<<<<<< Updated upstream
  /// endcap
  hasStub_st1 = false; 
  hasStub_st2 = false; 
  hasStub_st3 = false; 
  hasStub_st4 = false; 
=======

>>>>>>> Stashed changes*/
  hasGEMPad_st1 = false; 
  hasGEMPad_st2 = false; 
  npar = -1;
  meRing = -1;
  //position-based
  ddY123 = -99;
  deltaY12 = -99;
  deltaY23 = -99;

  //direction-based
  phi_ge21 = -9;
  phiM_st1 = -9;
  phiM_st2 = -9;
  phiM_st12 = -9;
  phiM_st23 = -9;
  dPhi_dir_st1_st2 = -9;
  dPhi_dir_st1_st12 = -9;
  dPhi_dir_st2_st23 = -9;
  dPhi_dir_st12_st23 = -9;

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
  if (verbose_>=0) std::cout <<"detid "<< ch_id<<" perp "<< perp << std::endl;
  // use average perp  
  //perp = perp/phis.size();
  // do a fit to the comparator digis
  float alpha = 0., beta = 0.;
  calculateAlphaBeta(zs, phis, ezs, ephis, status, alpha, beta);
  //std::cout <<" alpha "<< alpha <<" beta "<< beta << std::endl;
  for (int i=0; i<6; i++){
      fit_z_layers[i] = cscChamber->layer(i+1)->centerOfStrip(20).z();
      fit_phi_layers[i] = normalizePhi(alpha + beta * fit_z_layers[i]);
      if (verbose_>0)
      	std::cout <<"i "<< i <<" fit_z "<< fit_z_layers[i]<< " fit_phi "<< fit_phi_layers[i]<< std::endl;
  }
  /*
  fit_z_layer3 = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->centerOfStrip(20).z();
  
  fit_phi_layer3 = normalizePhi(alpha + beta * fit_z_layer3);
  
  if(verbose_>0) {
    std::cout << "Number of comparator digis used in the fit " << ezs.size() << std::endl;
    std::cout << "best CSC stub fit phi position (L1Only) " << fit_z_layer3 << " " << fit_phi_layer3 << std::endl;
  }

  // calculate the z position in L1 and L6
  float l1_z = cscChamber->layer(1)->centerOfStrip(20).z();
  float l6_z = cscChamber->layer(6)->centerOfStrip(20).z();
  fit_z_layer1 = l1_z;
  fit_z_layer6 = l6_z;
  fit_phi_layer1 = normalizePhi(alpha + beta * l1_z);
  fit_phi_layer6 = normalizePhi(alpha + beta * l6_z);
  */
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
  if (verbose_>=0)
      std::cout <<"LCT position, chid "<< chid <<" gp eta "<< gp_st_layer3[st-1].eta()<<" phi "<<gp_st_layer3[st-1].phi()<<" perp "<< gp_st_layer3[st-1].perp() << std::endl;
 /* 
  if (chid.station() == 1){
      gp_st1_layer1 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][0], z_st_layers[st-1][0]));
      gp_st1 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][2], z_st_layers[st-1][2]));
      gp_st1_layer6 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][5], z_st_layers[st-1][5]));
      if (verbose_>=0)
      	std::cout <<"LCT position st1 chid "<< chid <<" gp eta "<< gp_st1.eta()<<" phi "<<gp_st1.phi()<<" perp "<< gp_st1.perp() << std::endl;
  }else if (chid.station() == 2){
      gp_st2_layer1 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][0], z_st_layers[st-1][0]));
      gp_st2 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][2], z_st_layers[st-1][2]));
      gp_st2_layer6 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][5], z_st_layers[st-1][5]));
      if (verbose_>=0)
      	std::cout <<"LCT position st2 chid "<< chid <<" gp eta "<< gp_st2.eta()<<" phi "<<gp_st2.phi() <<" perp "<< gp_st2.perp() << std::endl;
  }else if (chid.station() == 3){
      gp_st3 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][2], z_st_layers[st-1][2]));
      if (verbose_>=0)
      	std::cout <<"LCT position st3 chid "<< chid <<" gp eta "<< gp_st3.eta()<<" phi "<<gp_st3.phi() <<" perp "<< gp_st3.perp() << std::endl;
  }else if (chid.station() == 4){
      gp_st4 = GlobalPoint(GlobalPoint::Cylindrical(perp, phi_st_layers[st-1][2], z_st_layers[st-1][2]));
      if (verbose_>=0)
      	std::cout <<"LCT position st4 chid "<< chid <<" gp eta "<< gp_st3.eta()<<" phi "<<gp_st3.phi() <<" perp "<< gp_st3.perp() << std::endl;
  }else if (verbose_>0) 
      std::cout <<" not in CSC station 1 , 2 ,3 , 4, chamber id  "<< chid << std::endl;
 */ 

}

void DisplacedMuonTriggerPtassignment::globalPositionOfGEMPad(const GEMCSCPadDigi gempad, GEMDetId gemid)
{
  GEMDetId ch_id(gemid.region(), gemid.ring(), gemid.station(), gemid.layer(), gemid.chamber(), 0);
  const GEMChamber* gemChamber(gemGeometry_->chamber(ch_id.rawId()));
  auto gemRoll(gemChamber->etaPartition(gemid.roll()));//any roll
  const int nGEMPads(gemRoll->npads());
  if (gempad.pad() > nGEMPads or gempad.pad() < 0){
      std::cout <<" gempad.pad() is within pad range gempad "<< gempad <<" npad "<< nGEMPads << std::endl;
      return;
  }

  const LocalPoint lpGEM(gemRoll->centreOfPad(gempad.pad()));
  if (gemid.station() == 1){
  	gp_ge11 = GlobalPoint(gemRoll->toGlobal(lpGEM));
	if (verbose_>0) std::cout <<" gempad in GE11 id " << gemid <<" gp eta "<< gp_ge11.eta()<<" phi "<< gp_ge11.phi()<<" pad "<<gempad.pad()<< std::endl;
  }else if (gemid.station() == 3){
  	gp_ge21 = GlobalPoint(gemRoll->toGlobal(lpGEM));
	if (verbose_>0) std::cout <<" gempad in GE21 id "<< gemid <<" gp eta "<< gp_ge21.eta()<<" phi "<< gp_ge21.phi()<<" pad "<<gempad.pad()<< std::endl;
  }else if (verbose_>0) 
      std::cout <<" gemid "<< gemid  <<" not in station 1 or 3" << std::endl;

}



int DisplacedMuonTriggerPtassignment::getEtaPartition(float eta) const
{
    int neta=-1;
    if (fabs(eta)>=1.2 and fabs(eta)<1.4)
	neta=0;
    else if (fabs(eta)>=1.4 and fabs(eta)<1.6)
	neta=1;
    else if (fabs(eta)>=1.6 and fabs(eta)<1.8)
	neta=2;
    else if (fabs(eta)>=1.8 and fabs(eta)<2.0)
	neta=3;
    else if (fabs(eta)>=2.0 and fabs(eta)<2.2)
	neta=4;
    else if (fabs(eta)>=2.2 and fabs(eta)<2.4)
	neta=5;

    return neta;

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
   int neta = getEtaPartition(eta);
   
   if (par<0 or par>3 or neta==-1) return -99;
   return (deltay23-PositionEpLUT[par][neta][0]*deltay12);
   
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
   return true;
}


bool DisplacedMuonTriggerPtassignment::runDirectionbasedGE21()
{
   if (not (npar<4 and npar>=0 and hasGEMPad_st1 and hasGEMPad_st2)) return false; 
   if (fabs(phi_ge21)>4) return false;//check this because we want to use setPhiGE21() to set phi_ge21 (using 2strips-pad)
   /*
   float xfactor_st1 = xfactor*fabs(gp_ge11.z() - gp_st1.z());
   float xfactor_st2 = xfactor*fabs(gp_ge21.z() - gp_st2.z())/(xfactor*fabs(gp_st1.z() - gp_st2.z())+1);
   float xfactor_st12 = xfactor*fabs(gp_st1.z() - gp_st2.z())/(xfactor*fabs(gp_st1.z() - gp_st2.z())+1);
   float xfactor_st23 = xfactor*fabs(gp_st2.z() - gp_st3.z())/(xfactor*fabs(gp_st1.z() - gp_st3.z())+1);
   if (verbose_>0) std::cout <<"DisplacedMuonTrigger meRing "<< meRing <<" xfactor st1 "<< xfactor_st1 <<" xfactor st2 "<< xfactor_st2 << std::endl;
   phiM_st1 = phiMomentum_Xfactor(gp_st1.phi(), gp_ge11.phi(), xfactor_st1);//
   phiM_st2 = phiMomentum_Xfactor(gp_st2.phi(), phi_ge21, xfactor_st2);
   phiM_st12 = phiMomentum_Xfactor(gp_st2.phi(), gp_st1.phi(), xfactor_st12);
   phiM_st23 = phiMomentum_Xfactor(gp_st3.phi(), gp_st2.phi(), xfactor_st23);
   */

   float xfactor_st1 = xfactor*fabs(gp_ge11.z() - gp_st_layer3[0].z());
   float xfactor_st2 = xfactor*fabs(gp_ge21.z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   float xfactor_st12 = xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   float xfactor_st23 = xfactor*fabs(gp_st_layer3[1].z() - gp_st_layer3[2].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[2].z())+1);
   phiM_st1 = phiMomentum_Xfactor(gp_st_layer3[0].phi(), gp_ge11.phi(), xfactor_st1);//
   phiM_st2 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), phi_ge21, xfactor_st2);
   phiM_st12 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), gp_st_layer3[0].phi(), xfactor_st12);
   phiM_st23 = phiMomentum_Xfactor(gp_st_layer3[2].phi(), gp_st_layer3[1].phi(), xfactor_st23);
   if (verbose_>=0)  std::cout <<"DisplacedMuonTrigger, direction with GE21, meRing "<< meRing <<" xfactor_st1 "<< xfactor_st1 <<" phiM_st1 "<< phiM_st1
       			<<" xfactor_st2 "<< xfactor_st2 <<" phiM_st2 "<< phiM_st2 << std::endl;


   dPhi_dir_st1_st2 = (fabs(phiM_st1)<4 and fabs(phiM_st2)<4)? deltaPhi(phiM_st1, phiM_st2):-9;
   dPhi_dir_st1_st12 = (fabs(phiM_st1)<4 and fabs(phiM_st12)<4)? deltaPhi(phiM_st1, phiM_st12):-9;
   dPhi_dir_st2_st23 = (fabs(phiM_st2)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st2, phiM_st23):-9;
   dPhi_dir_st12_st23 = (fabs(phiM_st12)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st12, phiM_st23):-9;
   return true;
}


bool DisplacedMuonTriggerPtassignment::runDirectionbasedCSConly() 
{

   //z_st_layers should be used at sim level, set Z and phi for layer1 and layer6 at sim level, or rebuild constructor?
   if (not ( npar<4 and npar>=0 and ((meRing==1 and hasGEMPad_st1) or meRing==2))) return false; 
   float xfactor_st1 = 0;
   if (meRing==1){
   	xfactor_st1 = xfactor*fabs(gp_ge11.z() - gp_st_layer3[0].z());
   	phiM_st1 = phiMomentum_Xfactor(gp_st_layer3[0].phi(), gp_ge11.phi(), xfactor_st1);//
   }else if (meRing == 2){
   	//xfactor_st1 = xfactor*fabs(gp_st1_layer1.z() - gp_st1_layer6.z())/(xfactor*fabs(gp_st1.z() - gp_st1_layer6.z())+1);
   	xfactor_st1 = xfactor*fabs(z_st_layers[0][0] - z_st_layers[0][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[0][5])+1);
   	phiM_st1 = phiMomentum_Xfactor(phi_st_layers[0][5], phi_st_layers[0][0], xfactor_st1);//
   }
   float xfactor_st2 = xfactor*fabs(z_st_layers[1][0] - z_st_layers[1][5])/(xfactor*fabs(gp_st_layer3[0].z() - z_st_layers[1][5])+1);
   float xfactor_st12 = xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[1].z())+1);
   float xfactor_st23 = xfactor*fabs(gp_st_layer3[1].z() - gp_st_layer3[2].z())/(xfactor*fabs(gp_st_layer3[0].z() - gp_st_layer3[2].z())+1);
   phiM_st2 = phiMomentum_Xfactor(phi_st_layers[1][5], phi_st_layers[1][0], xfactor_st2);
   phiM_st12 = phiMomentum_Xfactor(gp_st_layer3[1].phi(), gp_st_layer3[0].phi(), xfactor_st12);
   phiM_st23 = phiMomentum_Xfactor(gp_st_layer3[2].phi(), gp_st_layer3[1].phi(), xfactor_st23);
   //if phi in layer1 and layer6 in station1 and 2 are not set, then here phiM return -9
   if (verbose_>=0)  std::cout <<"DisplacedMuonTrigger CSConly direction, meRing "<< meRing <<" xfactor_st1 "<< xfactor_st1 <<" phiM_st1 "<< phiM_st1
       			<<" xfactor_st2 "<< xfactor_st2 <<" phiM_st2 "<< phiM_st2 << std::endl;
   
   dPhi_dir_st1_st2 = (fabs(phiM_st1)<4 and fabs(phiM_st2)<4)? deltaPhi(phiM_st1, phiM_st2):-9;
   dPhi_dir_st1_st12 = (fabs(phiM_st1)<4 and fabs(phiM_st12)<4)? deltaPhi(phiM_st1, phiM_st12):-9;
   dPhi_dir_st2_st23 = (fabs(phiM_st2)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st2, phiM_st23):-9;
   dPhi_dir_st12_st23 = (fabs(phiM_st12)<4 and fabs(phiM_st23)<4)? deltaPhi(phiM_st12, phiM_st23):-9;
   return true;
}



float DisplacedMuonTriggerPtassignment::getlocalPhiDirection(int st) const
{
    //st =1 :station1 , st=2: station2 
    //st = 12 : between station1 and station2; st = 23 : between station2 and station3
   bool hasGEMPad1((meRing==1 and hasGEMPad_st1) or meRing==2);
   bool hasGEMPad2((meRing==1 and hasGEMPad_st2) or meRing==2);
   if (st==1 and hasStub_st[0] and hasGEMPad1) return phiM_st1;
   else if (st==2 and hasStub_st[1] and hasGEMPad2) return phiM_st2;
   else if (st == 12 and hasStub_st[0] and hasStub_st[1]) return phiM_st12;
   else if (st == 23 and hasStub_st[1] and hasStub_st[2]) return phiM_st23;
   else{
   	std::cout <<" error in getlocalPhiDirection, st  "<<st <<" not in range or not not have stub or GEMpad" << std::endl;
	return -99; 
   }
}

float DisplacedMuonTriggerPtassignment::getdeltaPhiDirection(int st1, int st2) const
{
   bool hasGEMPad1((meRing==1 and hasGEMPad_st1) or meRing==2);
   bool hasGEMPad2((meRing==1 and hasGEMPad_st2) or meRing==2);
   if (((st1 == 1 and st2 == 2) or (st1 == 1 and st2 == 2)) and hasStub_st[0] and hasGEMPad1 and hasStub_st[1] and hasGEMPad2) return dPhi_dir_st1_st2;
   else if (((st1 == 1 and st2 == 12) or (st1 == 12 and st2 == 1)) and hasStub_st[0] and hasGEMPad1 and hasStub_st[1]) return dPhi_dir_st1_st12;
   else if (((st1 == 2 and st2 == 23) or (st1 == 23 and st2 == 2)) and hasStub_st[1] and hasGEMPad2 and hasStub_st[2]) return dPhi_dir_st2_st23;
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

#endif
