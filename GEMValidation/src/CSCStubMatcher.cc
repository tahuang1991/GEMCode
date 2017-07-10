#include "GEMCode/GEMValidation/interface/CSCStubMatcher.h"
#include "GEMCode/GEMValidation/interface/SimHitMatcher.h"
#include "L1Trigger/CSCTriggerPrimitives/src/CSCCathodeLCTProcessor.h"

#include <algorithm>

using namespace std;
using namespace matching;


CSCStubMatcher::CSCStubMatcher(SimHitMatcher& sh, CSCDigiMatcher& dg, GEMDigiMatcher& gem_dg, RPCDigiMatcher& rpc_dg,
                               edm::EDGetTokenT<CSCCLCTDigiCollection>& clctInputs_, 
                               edm::EDGetTokenT<CSCALCTDigiCollection>& alctInputs_, 
                               edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection>& lctInputs_, 
                               edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection>& mplctInputs_)
: DigiMatcher(sh)
, digi_matcher_(&dg)
, gem_digi_matcher_(&gem_dg)
, rpc_digi_matcher_(&rpc_dg)
, sh_matcher_(&sh)
{
  auto cscCLCT_ = conf().getParameter<edm::ParameterSet>("cscCLCT");
  minBXCLCT_ = cscCLCT_.getParameter<int>("minBX");
  maxBXCLCT_ = cscCLCT_.getParameter<int>("maxBX");
  verboseCLCT_ = cscCLCT_.getParameter<int>("verbose");
  minNHitsChamberCLCT_ = cscCLCT_.getParameter<int>("minNHitsChamber");
  runCLCT_ = cscCLCT_.getParameter<bool>("run");

  auto cscALCT_ = conf().getParameter<edm::ParameterSet>("cscALCT");
  minBXALCT_ = cscALCT_.getParameter<int>("minBX");
  maxBXALCT_ = cscALCT_.getParameter<int>("maxBX");
  verboseALCT_ = cscALCT_.getParameter<int>("verbose");
  minNHitsChamberALCT_ = cscALCT_.getParameter<int>("minNHitsChamber");
  runALCT_ = cscCLCT_.getParameter<bool>("run");

  auto cscLCT_ = conf().getParameter<edm::ParameterSet>("cscLCT");
  minBXLCT_ = cscLCT_.getParameter<int>("minBX");
  maxBXLCT_ = cscLCT_.getParameter<int>("maxBX");
  verboseLCT_ = cscLCT_.getParameter<int>("verbose");
  minNHitsChamberLCT_ = cscLCT_.getParameter<int>("minNHitsChamber");
  addGhostLCTs_ = cscLCT_.getParameter<bool>("addGhosts");
  matchAlctGemME11_ = cscLCT_.getParameter<bool>("matchAlctGemME11");
  matchAlctGemME21_ = cscLCT_.getParameter<bool>("matchAlctGemME21");
  matchClctGemME11_ = cscLCT_.getParameter<bool>("matchClctGemME11");
  matchClctGemME21_ = cscLCT_.getParameter<bool>("matchClctGemME21");
  matchAlctRpc_ = cscLCT_.getParameter<bool>("matchAlctRpc");
  matchClctRpc_ = cscLCT_.getParameter<bool>("matchClctRpc");
  hsFromSimHitMean_ = cscLCT_.getParameter<bool>("hsFromSimHitMean");
  runLCT_ = cscLCT_.getParameter<bool>("run");

  auto cscMPLCT_ = conf().getParameter<edm::ParameterSet>("cscMPLCT");
  minBXMPLCT_ = cscMPLCT_.getParameter<int>("minBX");
  maxBXMPLCT_ = cscMPLCT_.getParameter<int>("maxBX");
  verboseMPLCT_ = cscMPLCT_.getParameter<int>("verbose");
  minNHitsChamberMPLCT_ = cscMPLCT_.getParameter<int>("minNHitsChamber");
  addGhostMPLCTs_ = cscMPLCT_.getParameter<bool>("addGhosts");
  runMPLCT_ = cscMPLCT_.getParameter<bool>("run");

  minNHitsChamber_ = conf().getUntrackedParameter<int>("minNHitsChamber", 4);

  if (hasCSCGeometry_) {
    edm::Handle<CSCCLCTDigiCollection> clcts;
    if (gemvalidation::getByToken(clctInputs_, clcts, event())) if (runCLCT_) matchCLCTsToSimTrack(*clcts.product());    
    
    edm::Handle<CSCALCTDigiCollection> alcts;
    if (gemvalidation::getByToken(alctInputs_, alcts, event())) if (runALCT_) matchALCTsToSimTrack(*alcts.product());    
    
    edm::Handle<CSCCorrelatedLCTDigiCollection> lcts;
    if (gemvalidation::getByToken(lctInputs_, lcts, event())) if (runLCT_) matchLCTsToSimTrack(*lcts.product());    
    
    edm::Handle<CSCCorrelatedLCTDigiCollection> mplcts;
    if (gemvalidation::getByToken(mplctInputs_, mplcts, event())) if (runMPLCT_) matchMPLCTsToSimTrack(*mplcts.product());    
  }
}


CSCStubMatcher::~CSCStubMatcher() {}


void
CSCStubMatcher::matchCLCTsToSimTrack(const CSCCLCTDigiCollection& clcts)
{
  // only look for stub in chambers that have digis matching to this track
  setVerbose(verboseCLCT_);

  auto cathode_ids = digi_matcher_->chamberIdsStrip(0);
  int n_minLayers = 0;
  for (auto id: cathode_ids)
  {
    CSCDetId ch_id(id);
    if (digi_matcher_->nLayersWithStripInChamber(id) >= minNHitsChamberCLCT_) ++n_minLayers;

    // fill 1 half-strip wide gaps
    auto digi_strips = digi_matcher_->stripsInChamber(id, 1);
    if (verbose())
    {
      cout<<"clct: digi_strips "<<ch_id<<" ";
      copy(digi_strips.begin(), digi_strips.end(), ostream_iterator<int>(cout, " ")); cout<<endl;
    }

    auto clcts_in_det = clcts.get(ch_id);
    for (auto c = clcts_in_det.first; c != clcts_in_det.second; ++c)
    {
      if (!c->isValid()) continue;

      if (verbose()) cout<<"clct "<<ch_id<<" "<<*c<<endl;

      // check that the BX for this stub wasn't too early or too late
      if (c->getBX() < minBXCLCT_ || c->getBX() > maxBXCLCT_) continue;

      int half_strip = c->getKeyStrip() + 1; // CLCT halfstrip numbers start from 0
      auto mydigi = make_digi(id, half_strip, c->getBX(), CSC_CLCT, c->getQuality(), c->getPattern());

      // store all CLCTs in this chamber
      chamber_to_clcts_all_[id].push_back(mydigi);
      chamber_to_cscClcts_all_[id].push_back(*c);

      // match by half-strip with the digis
      if (digi_strips.find(half_strip) == digi_strips.end())
      {
        if (verbose()) cout<<"clctBAD"<<endl;
        continue;
      }
      if (verbose()) cout<<"clctGOOD"<<endl;


      // store matching CLCTs in this chamber
      chamber_to_clcts_[id].push_back(mydigi);
      chamber_to_cscClcts_[id].push_back(*c);

      if (chamber_to_clct_.find(id) != chamber_to_clct_.end())
      {
        //cout<<"WARNING!!! there already was matching CLCT "<<chamber_to_clct_[id]<<endl;
        //cout<<"   new digi: "<<mydigi<<endl;
	//sh_matcher_.	
	//GlobalPoint keygp(sh_matcher_.simHitPositionKeyLayer(id));
  	//LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  	//GlobalPoint csc_gp = cscGeometry_->idToDet(key_id)->surface().toGlobal(csc_intersect);
        // decide which one to choose
        int q_old = digi_quality(chamber_to_clct_[id]);
        int q_new = digi_quality(mydigi);

        if (q_old > q_new) continue; // keep old
        else if (q_old == q_new)
        {
          int p_old = digi_pattern(chamber_to_clct_[id]);
          int p_new = digi_pattern(mydigi);
          if (p_old > p_new) continue; // keep old
        }
        //cout<<"   new chosen"<<endl;
      }

      chamber_to_clct_[id] = mydigi;

    }
    if (chamber_to_clcts_[id].size() > 2)
    {
      //cout<<"WARNING!!! too many CLCTs "<<chamber_to_clcts_[id].size()<<" in "<<ch_id<<endl;
      //for (auto &c: chamber_to_clcts_[id]) cout<<"  "<<c<<endl;
    }
  }

  if (verbose() and n_minLayers > 0)
  {
    if (chamber_to_clct_.size() == 0)
    {
      cout<<"effNoCLCT"<<endl;
      for (const auto &it: clcts)
      {
        CSCDetId id(it.first);
        if (useCSCChamberType(id.iChamberType())) continue;
        auto clcts_in_det = clcts.get(id);
        for (auto c = clcts_in_det.first; c != clcts_in_det.second; ++c)
        {
          if (!c->isValid()) continue;
          if (verbose()) cout<<" clct: "<<id<<"  "<<*c<<endl;
        }
      }
    }
    else cout<<"effYesCLCT"<<endl;
  }
}


void
CSCStubMatcher::matchALCTsToSimTrack(const CSCALCTDigiCollection& alcts)
{
  setVerbose(verboseALCT_);
  // only look for stub in chambers that have digis matching to this track

  auto anode_ids = digi_matcher_->chamberIdsWire(0);
  int n_minLayers = 0;
  for (auto id: anode_ids)
  {
    if (digi_matcher_->nLayersWithWireInChamber(id) >= minNHitsChamberALCT_) ++n_minLayers;
    CSCDetId ch_id(id);

    // fill 1 WG wide gaps
    auto digi_wgs = digi_matcher_->wiregroupsInChamber(id, 1);
    if (verbose())
    {
      cout<<"alct: digi_wgs "<<ch_id<<" ";
      copy(digi_wgs.begin(), digi_wgs.end(), ostream_iterator<int>(cout, " ")); cout<<endl;
    }

    auto alcts_in_det = alcts.get(ch_id);
    for (auto a = alcts_in_det.first; a != alcts_in_det.second; ++a)
    {
      if (!a->isValid()) continue;

      if (verbose()) 
	  cout<<"alct "<<ch_id<<" "<<*a<<endl;

      // check that the BX for stub wasn't too early or too late
      if (a->getBX() < minBXALCT_ || a->getBX() > maxBXALCT_) continue;

      int wg = a->getKeyWG() + 1; // as ALCT wiregroups numbers start from 0
      auto mydigi = make_digi(id, wg, a->getBX(), CSC_ALCT, a->getQuality());

      // store all ALCTs in this chamber
      chamber_to_alcts_all_[id].push_back(mydigi);
      chamber_to_cscAlcts_all_[id].push_back(*a);

      // match by wiregroup with the digis
      if (digi_wgs.find(wg) == digi_wgs.end())
      {
        if (verbose()) cout<<"alctBAD"<<endl;
        continue;
      }
      if (verbose()) cout<<"alctGOOD"<<endl;

      // store matching ALCTs in this chamber
      chamber_to_alcts_[id].push_back(mydigi);
      chamber_to_cscAlcts_[id].push_back(*a);

      if (chamber_to_alct_.find(id) != chamber_to_alct_.end())
      {
        //cout<<"WARNING!!! there already was matching ALCT "<<chamber_to_alct_[id]<<endl;
        //cout<<"   new digi: "<<mydigi<<endl;

        // decide which one to choose
        int q_old = digi_quality(chamber_to_alct_[id]);
        int q_new = digi_quality(mydigi);
        if (q_old > q_new) continue; // keep old
        //cout<<"   new chosen"<<endl;
      }

      chamber_to_alct_[id] = mydigi;

    }
    if (chamber_to_alcts_[id].size() > 2)
    {
      //cout<<"WARNING!!! too many ALCTs "<<chamber_to_alcts_[id].size()<<" in "<<ch_id<<endl;
      //for (auto &a: chamber_to_alcts_[id]) cout<<"  "<<a<<endl;
    }
  }

  if (verbose() and n_minLayers > 0)
  {
    if (chamber_to_alct_.size() == 0)
    {
      cout<<"effNoALCT"<<endl;
      for (const auto &it: alcts)
      {
        CSCDetId id(it.first);
        if (useCSCChamberType(id.iChamberType())) continue;
        auto alcts_in_det = alcts.get(id);
        for (auto a = alcts_in_det.first; a != alcts_in_det.second; ++a)
        {
          if (!a->isValid()) continue;
          if (verbose()) cout<<" alct: "<<id<<"  "<<*a<<endl;
        }
      }
    }
    else cout<<"effYesALCT"<<endl;
  }
}


void
CSCStubMatcher::matchLCTsToSimTrack(const CSCCorrelatedLCTDigiCollection& lcts)
{
  setVerbose(verboseLCT_);
  // only look for stubs in chambers that already have CLCT and ALCT
  auto cathode_ids = chamberIdsAllCLCT(0);
  auto anode_ids = chamberIdsAllALCT(0);

  std::set<int> cathode_and_anode_ids;
  std::set_union(
      cathode_ids.begin(), cathode_ids.end(),
      anode_ids.begin(), anode_ids.end(),
      std::inserter(cathode_and_anode_ids, cathode_and_anode_ids.end())
  );

  int n_minLayers = 0;
  for (auto id: cathode_and_anode_ids)
  {
    if (digi_matcher_->nLayersWithStripInChamber(id) >= minNHitsChamberCLCT_ and digi_matcher_->nLayersWithWireInChamber(id) >= minNHitsChamberALCT_) ++n_minLayers;
    CSCDetId ch_id(id);

    auto lcts_in_det = lcts.get(ch_id);
    DigiContainer lcts_tmp;
    CSCCorrelatedLCTDigiContainer cscLcts_tmp;
    map<int, DigiContainer> bx_to_lcts;
    for (auto lct = lcts_in_det.first; lct != lcts_in_det.second; ++lct)
    {
      if (!lct->isValid()) continue;

      if (verbose()) cout<<"\n lct in detId "<<ch_id<<" "<<*lct<<endl;

      int bx = lct->getBX();

      // check that the BX for stub wasn't too early or too late
      if (bx < minBXLCT_ || bx > maxBXLCT_) continue;

      int hs = lct->getStrip() + 1; // LCT halfstrip and wiregoup numbers start from 0
      int wg = lct->getKeyWG() + 1;

      //fixedME
      //float dphi = lct->getGEMDPhi();
      float dphi = -9;

      auto mydigi = make_digi(id, hs, bx, CSC_LCT, lct->getQuality(), lct->getPattern(), wg, dphi);
      lcts_tmp.push_back(mydigi);
      cscLcts_tmp.push_back(*lct);
      bx_to_lcts[bx].push_back(mydigi);

      // Add ghost LCTs when there are two in bx
      // and the two don't share half-strip or wiregroup
      // TODO: when GEMs would be used to resolve this, there might ned to be an option to turn this off!
      if (bx_to_lcts[bx].size() == 2 and addGhostLCTs_)
      {
        auto lct11 = bx_to_lcts[bx][0];
        auto lct22 = bx_to_lcts[bx][1];
        int wg1 = digi_wg(lct11);
        int wg2 = digi_wg(lct22);
        int hs1 = digi_channel(lct11);
        int hs2 = digi_channel(lct22);

        if ( ! (wg1 == wg2 || hs1 == hs2) )
        {
          auto lct12 = lct11;
          digi_wg(lct12) = wg2;
          lcts_tmp.push_back(lct12);
	  //CSCCorrelatedLCTDigi LCT12(0, const int valid, const int quality, const int keywire, const int strip, const int pattern, const int bend, const int bx, const int mpclink=0, const uint16_t bx0=0, const uint16_t syncErr=0, const uint16_t cscID=0)
	  CSCCorrelatedLCTDigi LCT12(0, 1, digi_quality(lct12), digi_wg(lct12), digi_channel(lct12), digi_pattern(lct12), 0, digi_bx(lct12));
	  cscLcts_tmp.push_back(LCT12);

          auto lct21 = lct22;
          digi_wg(lct21) = wg1;
          lcts_tmp.push_back(lct21);
	  CSCCorrelatedLCTDigi LCT21(0, 1, digi_quality(lct21), digi_wg(lct21), digi_channel(lct21), digi_pattern(lct21), 0, digi_bx(lct21));
	  cscLcts_tmp.push_back(LCT21);
	  if (verbose())
	      cout<<"added ghosts"<<endl<<lct11<<"    "<<lct22<<endl <<lct12<<"    "<<lct21<<endl;
        }
      }
    } // lcts_in_det

    size_t n_lct = lcts_tmp.size();
    if (verbose()) cout<< "number of lcts = "<<n_lct <<" from cscLcts_tmp "<< cscLcts_tmp.size() <<endl;
    if (n_lct == 0) continue; // no LCTs in this chamber

    // assign the non necessarily matching LCTs
    chamber_to_lcts_all_[id] = lcts_tmp;
    chamber_to_cscLcts_all_[id] = cscLcts_tmp;

    if (verbose() and !(n_lct == 1 || n_lct == 2 || n_lct == 4 ) )
    {
      cout<<"WARNING!!! weird #LCTs="<<n_lct;
      for (auto &s: lcts_tmp) cout<<"  "<<s<<endl;
      //continue;
    }

    auto clct(clctsInChamber(id));
    auto alct(alctsInChamber(id));
    const auto hits = sh_matcher_->hitsInChamber(id);
    bool hasPad = false;
    bool hasDigis = false;
    if ((ch_id.station()==1 or ch_id.station()==2) and ch_id.ring()==1){
	// find a matching LCT
	const GEMDetId gemDetId(GEMDetId(ch_id.zendcap(),1,ch_id.station(),0,ch_id.chamber(),0));
	auto pads(gem_digi_matcher_->coPadsInSuperChamber(gemDetId));
	hasPad = (pads.size()>0);
	if (verbose())
	    std::cout <<"GEMDetId "<< gemDetId << (hasPad ? " has CoPad ":" has NO CoPad") << std::endl;
    }
    if ((ch_id.station()==3 or ch_id.station()==4) and ch_id.ring()==1){
	// find matching rpc chamber (only valid for me31 and me41)
	const int csc_trig_sect(CSCTriggerNumbering::triggerSectorFromLabels(ch_id));
	const int csc_trig_id( CSCTriggerNumbering::triggerCscIdFromLabels(ch_id));
	const int csc_trig_chid((3*(csc_trig_sect-1)+csc_trig_id)%18 +1);
	const int rpc_trig_sect((csc_trig_chid-1)/3+1);
	const int rpc_trig_subsect((csc_trig_chid-1)%3+1);
	const RPCDetId rpcDetId(RPCDetId(ch_id.zendcap(),1,ch_id.station(),rpc_trig_sect,1,rpc_trig_subsect,0));
	auto rpcDigis(rpc_digi_matcher_->digisInChamber(rpcDetId));
	hasDigis = (rpcDigis.size()>0);
    }


    
    
    int iLct = -1;
    for (auto &lct: lcts_tmp)
    {
      iLct++;
      bool lct_matched(false);
      for (unsigned int j=0; j<alct.size();j++){
        for (unsigned int i=0; i<clct.size()+1;i++){
          
          if (verbose()) 
            std::cout <<"lct size "<< lcts_tmp.size()<<" ilct "<< iLct << " available LCT " << lct << std::endl;
          if (verbose())
            std::cout <<"alcts size "<< alct.size() <<" jth "<<j<< " available ALCT " << alct[j] << std::endl;
          if (verbose() and i<clct.size()) 
            std::cout<<"clcts size"<< clct.size()<<" ith "<<i << " available CLCT " << clct[i] << std::endl; 
          
          const bool caseAlctClct(j<alct.size() && i<clct.size());
          const bool caseAlctGem(is_valid(alct[j]) and i==clct.size() and (ch_id.station() == 1 or ch_id.station() == 2) and ch_id.ring()==1);
          //const bool caseClctGem(is_valid(clct[i]) and hasPad and !is_valid(alct[j]) and (ch_id.station() == 1 or ch_id.station() == 2));
          const bool caseAlctRpc(is_valid(alct[j]) and i==clct.size() and (ch_id.station() == 3 or ch_id.station() == 4) and ch_id.ring()==1);
          //const bool caseClctRpc(is_valid(clct[i]) and hasDigis and !is_valid(alct[j]) and (ch_id.station() == 3 or ch_id.station() == 4));
          if(verbose()){
            if (caseAlctGem and hasPad) std::cout <<" caseAlctGem and hasPad "<< std::endl;
            else if (caseAlctGem and not(hasPad)) cout <<" caseAlct Gem and not hasPad "<< std::endl;
            if (caseAlctRpc and hasDigis) std::cout <<" caseAlctRpc and hasDigis "<< std::endl;
            else if (caseAlctRpc and not(hasDigis)) cout <<" caseAlctRpc and not hasDigis "<< std::endl;
          }
          const CSCChamber* cscChamber(cscGeometry_->chamber(CSCDetId(id)));
          
          const CSCLayer* cscKeyLayer(cscChamber->layer(3));
          const CSCLayerGeometry* cscKeyLayerGeometry(cscKeyLayer->geometry());
          const int nStrips(cscKeyLayerGeometry->numberOfStrips());
          const float averageZ((cscKeyLayer->centerOfStrip(0)).z());
          auto GpME(propagateToZ(averageZ));
          auto lpME(cscKeyLayer->toLocal(GpME));
          
          const int my_hs_gem_propagate((nStrips-cscKeyLayerGeometry->nearestStrip(lpME))*2);
          const auto hits = sh_matcher_->hitsInChamber(id);
          const float my_hs_gemrpc_mean(sh_matcher_->simHitsMeanStrip(hits));
          //const float my_wg_gemrpc_mean(sh_matcher_->simHitsMeanStrip(hits));
          bool matchAlctGem_((matchAlctGemME11_ and ch_id.station()==1 )||(matchAlctGemME21_ and ch_id.station()==2));
          //bool matchClctGem_((matchClctGemME11_ and cscid.station()==1 )||(matchClctGemME21_ and cscid.station()==2));
          /* Printing Stuff
             if (caseAlctClct[i][j]) std::cout << "caseAlctClct" << std::endl;
             else if(matchAlctGem_)std::cout << "caseAlctGem" << std::endl;
             std::cout << "mean half strip from simhits " << sh_matcher_->simHitsMeanStrip(hits)
             <<"   half strip by propagating track " << my_hs_gem_propagate << std::endl;
          */
          float my_hs_gemrpc;
          if (hsFromSimHitMean_)  my_hs_gemrpc= my_hs_gemrpc_mean;
          else my_hs_gemrpc = my_hs_gem_propagate;
          
          //All the matching information can go here.
          
          auto digi_strips = digi_matcher_->stripsInChamber(id, 1);
          int my_hs = -1 ;
          if ( i < clct.size() ) my_hs = digi_channel(clct[i]);
          
          const int my_wg(digi_wg(alct[j]));
          const int my_bx(digi_bx(alct[j]));
          auto digi_wgs = digi_matcher_->wiregroupsInChamber(id,1);
          
	  if (not(caseAlctClct) and not(caseAlctGem) and not(caseAlctRpc)){
	      if (verbose())
	      	std::cout <<"this LCT can not be match AlctClct, AlctGem, AlctRpc, skip it: "<< lct<< std::endl;
	      continue;
	  }
            
          if ( caseAlctClct and !((my_bx == digi_bx(lct) || 6 == digi_bx(lct)) and my_hs == digi_channel(lct) and my_wg == digi_wg(lct))){
            if (verbose()) cout<<"  BAD LCT in AlctClct case"<<endl;
            continue;
          }
          
          //add hadPad here
          if (matchAlctGem_ and caseAlctGem and !(hasPad and  my_bx == digi_bx(lct) and std::abs(my_hs_gemrpc - digi_channel(lct))<3 and my_wg == digi_wg(lct) ) ){
            if (verbose()) cout<<"  BAD LCT in AlctGem case"<<endl;
            continue;
          }
          else if (caseAlctGem and !matchAlctGem_) continue;
          
          /* if (matchClctGem_ and caseClctGem and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and std::abs(my_wg_gemrpc_mean - digi_wg(lct))<3) ){
             if (verbose()) cout<<"  BAD"<<endl;
             continue;
             }*/
          
          //if (matchAlctRpc_ and caseAlctRpc and !(my_bx == digi_bx(lct) and std::abs(my_hs_gemrpc - digi_channel(lct))<3 and my_wg == digi_wg(lct) ) ){
          if (matchAlctRpc_ and caseAlctRpc and !(hasDigis and  std::fabs(my_hs_gemrpc - digi_channel(lct))<3.0 and my_wg == digi_wg(lct) ) ){
            if (verbose()) cout<<"  BAD LCT in AlctRpc case"<<endl;
            //std::cout << "my_hs_gemrpc " << my_hs_gemrpc << "  hs from LCT " << digi_channel(lct) 
            //        <<" my_wg " << my_wg << " wg in LCT "<<digi_wg(lct) <<std::endl;
            continue;
          }
          else if (caseAlctRpc and !matchAlctRpc_) continue;
          /* if (matchClctRpc_ and caseClctRpc and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and std::abs(my_wg_gemrpc_mean - digi_wg(lct))<3) ){
             if (verbose()) cout<<"  BAD"<<endl;
             continue;
             }*/
          
          if (verbose() and ch_id.ring()==1 and (ch_id.station()==1 or ch_id.station()==2) and !caseAlctClct and !caseAlctGem) 
            std::cout <<"ME11 or ME21: jth alct "<< j<<" ith clct "<<i <<"  not caseAlctClct caseAlctGem, LCT "<< lct << std::endl; 
          if (chamber_to_lct_.find(id) == chamber_to_lct_.end())   chamber_to_lct_[id] = lct;
          else if (chamber_to_lct_.find(id) != chamber_to_lct_.end() and 
                   fabs(digi_channel(chamber_to_lct_[id])-my_hs_gemrpc) > fabs(digi_channel(lct)-my_hs_gemrpc)){
            if (verbose()){
              cout<<"ALARM!!! here already was matching LCT "<<chamber_to_lct_[id]<<endl;
              cout<<"   new digi: "<<lct<<endl;
            }
            chamber_to_lct_[id] = lct;
          }
          
          chamber_to_lcts_[id].push_back(lct);
          chamber_to_cscLcts_[id].push_back(cscLcts_tmp[iLct]);
          lct_matched = true;
          
          if (verbose() and caseAlctClct)  std::cout << " this LCT is matched to simtrack in AlctClct case" << std::endl;
          if (verbose() and caseAlctGem)  std::cout << " this LCT is matched to simtrack in AlctGem case" << std::endl;
          if (verbose() and caseAlctRpc)  std::cout << " this LCT is matched to simtrack in AlctRpc case" << std::endl;
	 
          break;
        } //clct loop over
        if (lct_matched) break;
      }// alct loop over
      
    } // lct loop over
    
  }
  
  if (verbose() and n_minLayers > 0)
  {
    if (chamber_to_lct_.size() == 0)
    {
      cout<<"No Matched LCT"<<endl;
      for (const auto &it: lcts)
      {
        CSCDetId id(it.first);
        if (useCSCChamberType(id.iChamberType())) continue;
        auto lcts_in_det = lcts.get(id);
        for (auto a = lcts_in_det.first; a != lcts_in_det.second; ++a)
        {
          if (!a->isValid()) continue;
          if (verbose()) cout<<" lct: "<<id<<"  "<<*a<<endl;
        }
      }

    }
    else cout<<"at least one matched LCT" << std::endl;
  }
}


void
CSCStubMatcher::matchMPLCTsToSimTrack(const CSCCorrelatedLCTDigiCollection& mplcts)
{
  setVerbose(verboseMPLCT_);
  // only look for stubs in chambers that already have CLCT and ALCT
  auto cathode_ids = chamberIdsAllCLCT(0);
  auto anode_ids = chamberIdsAllALCT(0);

  std::set<int> cathode_and_anode_ids;
  std::set_union(
      cathode_ids.begin(), cathode_ids.end(),
      anode_ids.begin(), anode_ids.end(),
      std::inserter(cathode_and_anode_ids, cathode_and_anode_ids.end())
  );

  int n_minLayers = 0;
  for (auto id: cathode_and_anode_ids)
  {
    if (digi_matcher_->nLayersWithStripInChamber(id) >= minNHitsChamberCLCT_ and digi_matcher_->nLayersWithWireInChamber(id) >= minNHitsChamberALCT_) ++n_minLayers;
    CSCDetId ch_id(id);

    auto mplcts_in_det = mplcts.get(ch_id);
    DigiContainer mplcts_tmp;
    CSCCorrelatedLCTDigiContainer cscMplcts_tmp;
    map<int, DigiContainer> bx_to_mplcts;
    for (auto lct = mplcts_in_det.first; lct != mplcts_in_det.second; ++lct)
    {
      if (!lct->isValid()) continue;

      if (verbose()) cout<<"mplct in detId"<<ch_id<<" "<<*lct<<endl;

      int bx = lct->getBX();

      // check that the BX for stub wasn't too early or too late
      if (bx < minBXLCT_ || bx > maxBXLCT_) continue;

      int hs = lct->getStrip() + 1; // LCT halfstrip and wiregoup numbers start from 0
      int wg = lct->getKeyWG() + 1;

      //fixedME
      //float dphi = lct->getGEMDPhi();
      float dphi = -9;

      auto mydigi = make_digi(id, hs, bx, CSC_LCT, lct->getQuality(), lct->getPattern(), wg, dphi);
      mplcts_tmp.push_back(mydigi);
      cscMplcts_tmp.push_back(*lct);
      bx_to_mplcts[bx].push_back(mydigi);

      // Add ghost mplcts when there are two in bx
      // and the two don't share half-strip or wiregroup
      // TODO: when GEMs would be used to resolve this, there might ned to be an option to turn this off!
      if (bx_to_mplcts[bx].size() == 2 and addGhostMPLCTs_)
      {
        auto lct11 = bx_to_mplcts[bx][0];
        auto lct22 = bx_to_mplcts[bx][1];
        int wg1 = digi_wg(lct11);
        int wg2 = digi_wg(lct22);
        int hs1 = digi_channel(lct11);
        int hs2 = digi_channel(lct22);

        if ( ! (wg1 == wg2 || hs1 == hs2) )
        {
          auto lct12 = lct11;
          digi_wg(lct12) = wg2;
          mplcts_tmp.push_back(lct12);

          auto lct21 = lct22;
          digi_wg(lct21) = wg1;
          mplcts_tmp.push_back(lct21);
          //cout<<"added ghosts"<<endl<<lct11<<"    "<<lct22<<endl <<lct12<<"    "<<lct21<<endl;
        }
      }
    } // mplcts_in_det

    size_t n_lct = mplcts_tmp.size();
    if (verbose()) cout<<"number of mplct = "<<n_lct<<endl;
    if (n_lct == 0) continue; // no mplcts in this chamber

    // assign the non necessarily matching Mplcts
    chamber_to_mplcts_all_[id] = mplcts_tmp;
    chamber_to_cscMplcts_all_[id] = cscMplcts_tmp;

    if (verbose() and !(n_lct == 1 || n_lct == 2 || n_lct == 4 ) )
    {
      cout<<"WARNING!!! weird #Mplcts="<<n_lct;
      for (auto &s: mplcts_tmp) cout<<"  "<<s<<endl;
      //continue;
    }

    // find a matching LCT

    auto clct(clctsInChamber(id));
    auto alct(alctsInChamber(id));

    for (unsigned int i=0; i<clct.size();i++){

        if (!is_valid(clct[i])) continue;


        for (unsigned int j=0; j<alct.size();j++){
            if(!is_valid(alct[j])) continue;


            int my_hs = digi_channel(clct[i]);
            int my_wg = digi_wg(alct[j]);
            int my_bx = digi_bx(alct[j]);

            if (verbose()) cout<<"will match hs"<<my_hs<<" wg"<<my_wg<<" bx"<<my_bx<<" to #lct "<<n_lct<<endl;
            for (auto &lct: mplcts_tmp)
            {
              if (verbose()) cout<<" corlct "<<lct;
              if ( is_valid(alct[j]) and is_valid(clct[i]) and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and my_wg == digi_wg(lct)) ){
              if (verbose()) cout<<"  BAD"<<endl;
                continue;
                 }  
              if (verbose()) cout<<"  GOOD"<<endl;

              if (chamber_to_mplct_.find(id) != chamber_to_mplct_.end())
                {
                //cout<<"ALARM!!! there already was matching LCT "<<chamber_to_mplct_[id]<<endl;
            //cout<<"   new digi: "<<lct<<endl;
                }
            chamber_to_mplct_[id] = lct;

      // assign the matching Mplcts
              chamber_to_mplcts_[id].push_back(lct);
            }
        }//End of ALCT loop 
    } // End of CLCT loop
    
  }

  if (verbose() and n_minLayers > 0)
  {
    if (chamber_to_mplct_.size() == 0)
    {
      cout<<"effNoLCT"<<endl;
      for (const auto &it: mplcts)
      {
        CSCDetId id(it.first);
        if (useCSCChamberType(id.iChamberType())) continue;
        auto mplcts_in_det = mplcts.get(id);
        for (auto a = mplcts_in_det.first; a != mplcts_in_det.second; ++a)
        {
          if (!a->isValid()) continue;
          if (verbose()) cout<<" lct: "<<id<<"  "<<*a<<endl;
        }
      }

    }
    else cout<<"effYesLCT" << std::endl;
  }
}


std::set<unsigned int>
CSCStubMatcher::chamberIdsAllCLCT(int csc_type) const
{
  return selectDetIds(chamber_to_clcts_all_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsAllALCT(int csc_type) const
{
  return selectDetIds(chamber_to_alcts_all_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsAllLCT(int csc_type) const
{
  return selectDetIds(chamber_to_lcts_all_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsAllMPLCT(int csc_type) const
{
  return selectDetIds(chamber_to_mplcts_all_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsCLCT(int csc_type) const
{
  return selectDetIds(chamber_to_clct_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsALCT(int csc_type) const
{
  return selectDetIds(chamber_to_alct_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsLCT(int csc_type) const
{
  return selectDetIds(chamber_to_lct_, csc_type);
}

std::set<unsigned int>
CSCStubMatcher::chamberIdsMPLCT(int csc_type) const
{
  return selectDetIds(chamber_to_mplct_, csc_type);
}


matching::Digi
CSCStubMatcher::clctInChamber(unsigned int detid) const
{
  if (chamber_to_clct_.find(detid) == chamber_to_clct_.end()) return make_digi();
  return chamber_to_clct_.at(detid);
}

matching::Digi
CSCStubMatcher::alctInChamber(unsigned int detid) const
{
  if (chamber_to_alct_.find(detid) == chamber_to_alct_.end()) return make_digi();
  return chamber_to_alct_.at(detid);
}

matching::Digi
CSCStubMatcher::lctInChamber(unsigned int detid) const
{
  if (chamber_to_lct_.find(detid) == chamber_to_lct_.end()) return make_digi();
  return chamber_to_lct_.at(detid);
}

matching::Digi
CSCStubMatcher::mplctInChamber(unsigned int detid) const
{
  if (chamber_to_mplct_.find(detid) == chamber_to_mplct_.end()) return make_digi();
  return chamber_to_mplct_.at(detid);
}


const matching::DigiContainer&
CSCStubMatcher::allCLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_clcts_all_.find(detid) == chamber_to_clcts_all_.end()) return no_digis_;
  return chamber_to_clcts_all_.at(detid);
}

const matching::DigiContainer&
CSCStubMatcher::allALCTsInChamber(unsigned int detid) const
{
  if (chamber_to_alcts_all_.find(detid) == chamber_to_alcts_all_.end()) return no_digis_;
  return chamber_to_alcts_all_.at(detid);
}

const matching::DigiContainer&
CSCStubMatcher::allLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_lcts_all_.find(detid) == chamber_to_lcts_all_.end()) return no_digis_;
  return chamber_to_lcts_all_.at(detid);
}


const matching::DigiContainer&
CSCStubMatcher::allMPLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_mplcts_all_.find(detid) == chamber_to_mplcts_all_.end()) return no_digis_;
  return chamber_to_mplcts_all_.at(detid);
}

const matching::DigiContainer&
CSCStubMatcher::clctsInChamber(unsigned int detid) const
{
  if (chamber_to_clcts_.find(detid) == chamber_to_clcts_.end()) return no_digis_;
  return chamber_to_clcts_.at(detid);
}

const matching::DigiContainer&
CSCStubMatcher::alctsInChamber(unsigned int detid) const
{
  if (chamber_to_alcts_.find(detid) == chamber_to_alcts_.end()) return no_digis_;
  return chamber_to_alcts_.at(detid);
}

const matching::DigiContainer&
CSCStubMatcher::lctsInChamber(unsigned int detid) const
{
  if (chamber_to_lcts_.find(detid) == chamber_to_lcts_.end()) return no_digis_;
  return chamber_to_lcts_.at(detid);
}


const matching::DigiContainer&
CSCStubMatcher::mplctsInChamber(unsigned int detid) const
{
  if (chamber_to_mplcts_.find(detid) == chamber_to_mplcts_.end()) return no_digis_;
  return chamber_to_mplcts_.at(detid);
}


const matching::DigiContainer
CSCStubMatcher::lctsInStation(int st) const
{
  DigiContainer lcts;
  for (auto p : chamber_to_lcts_) 
  {
    CSCDetId id(p.first);
    if (id.station() == st) lcts.insert(lcts.end(), (p.second).begin(), (p.second).end());
    else continue; 
  }
  return lcts;
}


const CSCCLCTDigiContainer&
CSCStubMatcher::allCscCLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_cscClcts_all_.find(detid) == chamber_to_cscClcts_all_.end()) return no_csc_clcts_;
  return chamber_to_cscClcts_all_.at(detid);
}

const CSCALCTDigiContainer&
CSCStubMatcher::allCscALCTsInChamber(unsigned int detid) const
{
  if (chamber_to_cscAlcts_all_.find(detid) == chamber_to_cscAlcts_all_.end()) return no_csc_alcts_;
  return chamber_to_cscAlcts_all_.at(detid);
}

const CSCCorrelatedLCTDigiContainer&
CSCStubMatcher::allCscLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_cscLcts_all_.find(detid) == chamber_to_cscLcts_all_.end()) return no_csc_lcts_;
  return chamber_to_cscLcts_all_.at(detid);
}

const CSCCorrelatedLCTDigiContainer&
CSCStubMatcher::allCscMPLCTsInChamber(unsigned int detid) const
{
  if (chamber_to_cscMplcts_all_.find(detid) == chamber_to_cscMplcts_all_.end()) return no_csc_mplcts_;
  return chamber_to_cscMplcts_all_.at(detid);
}

const CSCCLCTDigiContainer&
CSCStubMatcher::cscClctsInChamber(unsigned int detid) const
{
  if (chamber_to_cscClcts_.find(detid) == chamber_to_cscClcts_.end()) return no_csc_clcts_;
  return chamber_to_cscClcts_.at(detid);
}

const CSCALCTDigiContainer&
CSCStubMatcher::cscAlctsInChamber(unsigned int detid) const
{
  if (chamber_to_cscAlcts_.find(detid) == chamber_to_cscAlcts_.end()) return no_csc_alcts_;
  return chamber_to_cscAlcts_.at(detid);
}

const CSCCorrelatedLCTDigiContainer&
CSCStubMatcher::cscLctsInChamber(unsigned int detid) const
{
  if (chamber_to_cscLcts_.find(detid) == chamber_to_cscLcts_.end()) return no_csc_lcts_;
  return chamber_to_cscLcts_.at(detid);
}

const CSCCorrelatedLCTDigiContainer&
CSCStubMatcher::cscMplctsInChamber(unsigned int detid) const
{
  if (chamber_to_cscMplcts_.find(detid) == chamber_to_cscMplcts_.end()) return no_csc_mplcts_;
  return chamber_to_cscMplcts_.at(detid);
}

CSCCLCTDigi 
CSCStubMatcher::bestCscClctInChamber(unsigned int detid) const
{
  //sort stubs based on quality
  auto input(cscClctsInChamber(detid));
  int bestQ = 0;
  int index = -1;
  for (unsigned int i=0; i<input.size(); ++i){
    int quality = input[i].getQuality();
    if (quality>bestQ){
      bestQ = quality;
      index = i;
    }
  }
  if (index != -1) return input[index];
  return CSCCLCTDigi();
}

CSCALCTDigi 
CSCStubMatcher::bestCscAlctInChamber(unsigned int detid) const
{
  //sort stubs based on quality
  auto input(cscAlctsInChamber(detid));
  int bestQ = 0;
  int index = -1;
  for (unsigned int i=0; i<input.size(); ++i){
    int quality = input[i].getQuality();
    if (quality>bestQ){
      bestQ = quality;
      index = i;
    }
  }
  if (index != -1) return input[index];
  return CSCALCTDigi();
}

CSCCorrelatedLCTDigi 
CSCStubMatcher::bestCscLctInChamber(unsigned int detid) const
{
 
  //sort stubs based on quality
  auto input(cscLctsInChamber(detid));
  int bestQ = 0;
  int index = -1;
  for (unsigned int i=0; i<input.size(); ++i){
    int quality = input[i].getQuality();
    if (quality>bestQ){
      bestQ = quality;
      index = i;
    }
  }
  if (index != -1) return input[index];
  return CSCCorrelatedLCTDigi();
  /* 
 // get the detid of the keylayer of this chamber
  auto chamberId = CSCDetId(detid);
  
  auto layer(cscGeometry_->chamber(chamberId)->layer(3));
  auto layerId(layer->id());
  auto hit_strips(sh_matcher_->hitStripsInDetId(layerId,0));
  std::cout << "In function bestCscLctInChamber " << std::endl;
  std::cout << "hit CSC strips (delta0) " << std::endl;
  copy(hit_strips.begin(), hit_strips.end(), ostream_iterator<int>(cout, " "));
  auto input(cscLctsInChamber(detid));
  int minDeltaStrip = 99;
  int index = -1;
  for (unsigned int i=0; i<input.size(); ++i){
    std::cout << input[i] << std::endl;
    int quality = input[i].getQuality();
    if (quality<minDeltaStrip){
      minDeltaStrip = quality;
      index = i;
    }
  }
  if (index != -1) return input[index];
  return CSCCorrelatedLCTDigi();
  */
}

CSCComparatorDigiDetIdContainer 
CSCStubMatcher::matchingComparatorDigisLCT(unsigned int detid, const CSCCorrelatedLCTDigi& stub) const
{
  //std::cout << "In function matchingComparatorDigisLCT" << std::endl;
  CSCComparatorDigiDetIdContainer output;
  //1) get the keystrip of the stub
  int keyStrip = stub.getStrip();
  int patternNumber = stub.getPattern();
  //2) loop on the layers
  auto chamberId = CSCDetId(detid);
  if (chamberId.ring()==4 and verboseLCT_)
  	std::cout << chamberId << std::endl;
  for (int ilayer=1; ilayer<=6; ilayer++){
    auto layer(cscGeometry_->chamber(chamberId)->layer(ilayer));
    auto layerId(layer->id());
    //std::cout << "\t" << layerId << std::endl;
    // get the comparator digis in this layer
    auto comps(digi_matcher_->cscComparatorDigisInDetId(layerId));
    // loop on the comparator digis
    for (auto comp: comps){
      //check if they match the LCT pattern
      int halfStrip = digi_matcher_->getHalfStrip(comp);
      if (chamberId.ring()==4 and verboseLCT_)
      	std::cout <<"layerid "<< layerId << "\t"<<comp<<" "<<halfStrip<<endl;
      // auto pattern = CSCCathodeLCTProcessor::pattern2007[patternNumber];
      // int actualLayer = ilayer - 1;
      // auto subPatternL0 = 
      if (this->comparatorInCLCTPattern(keyStrip, patternNumber, ilayer, halfStrip)){
        output.push_back(std::make_pair(layerId.rawId(), comp));
        //std::cout << "\t\tIt's in!" << std::endl;
      }
      else{
        //std::cout << "\t\tIt's out!" << std::endl;
      }
    }
  }
  return output;
}

bool 
CSCStubMatcher::comparatorInCLCTPattern(int keyStrip, int pattern, int layer, int halfStrip) const
{
  /*
  // first, get the (sub)pattern
  bool returnValue;
  std::vector<int> subpat = patIndexToPatternDelta[pattern].at(layer-1);
  for (auto p: subpat) std::cout << "\t" << p << std::endl;
  int halfStripDelta = halfStrip - keyStrip;
  std::cout << "\tkeyStrip pattern layer halfstrip " << keyStrip << " " <<pattern << " " <<layer << " " <<halfStrip <<std::endl <<std::endl;
  returnValue = std::find(subpat.begin(), subpat.end(), halfStripDelta) != subpat.end();
  */
  return true;//returnValue;
}


void 
CSCStubMatcher::positionsOfComparatorInLCT(unsigned int detid, const CSCCorrelatedLCTDigi& stub, std::vector<GlobalPoint>& gps) const
{
  auto digis = matchingComparatorDigisLCT (detid, stub);
  auto id = CSCDetId(detid);
  for (int l=1; l<=6; l++){
     CSCDetId l_id(id.endcap(), id.station(), id.ring(), id.chamber(), l); 
     float sum_x =0.0; 
     float sum_y = 0.0;
     float sum_z = 0.0;
     int n=0;
     //float z_pos = cscGeometry_->chamber(id)->layer(l)->centerOfStrip(20).z();
     for (auto p: digis){
	if (p.first == l_id.rawId()){
	    auto comp = p.second;
	    float fractional_strip = digi_matcher_->getFractionalStrip(comp);
	    auto layer_geo(cscGeometry_->chamber(id)->layer(l)->geometry());
	    float wire = layer_geo->middleWireOfGroup(stub.getKeyWG() + 1);
	    LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
	    GlobalPoint csc_gp = cscGeometry_->idToDet(l_id)->surface().toGlobal(csc_intersect);
	    //std::cout <<"In positionsOfComparatorInLCT local Point "<< csc_intersect <<" globalPoint "<< csc_gp << std::endl;
	    sum_x += csc_gp.x();
	    sum_y += csc_gp.y();
	    sum_z += csc_gp.z();
	    n++;
	    //std::cout <<"n "<< n <<" layerid "<< l_id <<" comparadigi hs "<< comp << " fractional strip "<< fractional_strip << " wg "<< wire <<" csc_phi "<< csc_gp.phi()<<" lp "<< csc_intersect.phi() << std::endl;
	}
     }
    // if (n>0 and verboseLCT_) std::cout <<"Digi cscid "<< l_id <<" n "<< n <<" x "<< sum_x/n <<" y "<< sum_y/n <<" z "<< sum_z/n << std::endl;
     if (n>0) gps.push_back(GlobalPoint(sum_x/n, sum_y/n, sum_z/n));
  }
  
}

//to do 
float 
CSCStubMatcher::zpositionOfLayer(unsigned int detid, int layer) const{

  auto id = CSCDetId(detid);
  auto cscChamber(cscGeometry_->chamber(id));
  return fabs(cscChamber->layer(layer)->centerOfStrip(20).z());

}

int
CSCStubMatcher::nChambersWithCLCT(int min_quality) const
{
  int result = 0;
  auto chamber_ids = chamberIdsCLCT();
  for (auto id: chamber_ids)
  {
    auto clct = clctInChamber(id);
    if (!is_valid(clct)) continue;
    if (digi_quality(clct) >= min_quality) ++result;
  }
  return result;
}

int
CSCStubMatcher::nChambersWithALCT(int min_quality) const
{
  int result = 0;
  auto chamber_ids = chamberIdsALCT();
  for (auto id: chamber_ids)
  {
    auto alct = alctInChamber(id);
    if (!is_valid(alct)) continue;
    if (digi_quality(alct) >= min_quality) ++result;
  }
  return result;
}

int
CSCStubMatcher::nChambersWithLCT(int min_quality) const
{
  int result = 0;
  auto chamber_ids = chamberIdsLCT();
  for (auto id: chamber_ids)
  {
    auto lct = lctInChamber(id);
    if (!is_valid(lct)) continue;
    if (digi_quality(lct) >= min_quality) ++result;
  }
  return result;
}

int
CSCStubMatcher::nChambersWithMPLCT(int min_quality) const
{
  int result = 0;
  auto chamber_ids = chamberIdsMPLCT();
  for (auto id: chamber_ids)
  {
    auto mplct = mplctInChamber(id);
    if (!is_valid(mplct)) continue;
    if (digi_quality(mplct) >= min_quality) ++result;
  }
  return result;
}


bool 
CSCStubMatcher::checkStubInChamber(CSCDetId id, CSCCorrelatedLCTDigi lct) const
{
  
  //    int hs = lct->getStrip() + 1; // LCT halfstrip and wiregoup numbers start from 0
  //    int wg = lct->getKeyWG() + 1;
  //    float dphi = lct->getGEMDPhi();
//  std::cout << " checks Stubs, DetId " << id << std::endl;
      //fixedME
      float dphi = -9;
  auto mydigi = make_digi(id.rawId(), lct.getStrip()+1, lct.getBX(), CSC_LCT, lct.getQuality(), lct.getPattern(),lct.getKeyWG()+1, dphi);
//      std::cout << "                   mydigi " << mydigi << std::endl;
  try{
  auto alldigis(chamber_to_lcts_all_.at(id.rawId()));
  for (auto p : alldigis) 
  {   
//      std::cout << " digi matched to simtrack " << p << std::endl;
      if (p==mydigi) return true;
   }
  }
  catch (exception& e)
    {
//	std::cout << " exception happens! " << std::endl;
    }
// std::cout << " matching is failed " << std::endl;
 return false;
}


bool 
CSCStubMatcher::wasStubMatchedInChamber(CSCDetId id, CSCCorrelatedLCTDigi lct) const
{
  auto stubs(cscLctsInChamber(id.rawId()));
  for (auto stub: stubs){
    if (stub==lct) return true;
  }
  return false;
}


GlobalPoint 
CSCStubMatcher::getGlobalPosition(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const
{
  // taken from https://github.com/cms-sw/cmssw/blob/dc9f78b6af4ad56c9342cf14041b6485a60b0691/L1Trigger/CSCTriggerPrimitives/src/CSCMotherboardME11GEM.cc
  CSCDetId cscId = CSCDetId(rawId);
  CSCDetId key_id(cscId.endcap(), cscId.station(), cscId.ring(), 
                  cscId.chamber(), CSCConstants::KEY_CLCT_LAYER);
  auto cscChamber = getCSCGeometry()->chamber(cscId);
  float fractional_strip = getFractionalStrip(lct);
  auto layer_geo = cscChamber->layer(CSCConstants::KEY_CLCT_LAYER)->geometry();
  // LCT::getKeyWG() also starts from 0
  float wire = layer_geo->middleWireOfGroup(lct.getKeyWG() + 1);
  LocalPoint csc_intersect = layer_geo->intersectionOfStripAndWire(fractional_strip, wire);
  GlobalPoint csc_gp = cscGeometry_->idToDet(key_id)->surface().toGlobal(csc_intersect);
  return csc_gp;

}


float 
CSCStubMatcher::getFractionalStrip(const CSCCorrelatedLCTDigi& lct) const
{
  return 0.5 * (lct.getStrip() + 1) - 0.25;
}


float 
CSCStubMatcher::getAverageBendingLCT(unsigned int rawId, const CSCCorrelatedLCTDigi& lct) const
{
  CSCDetId cscId = CSCDetId(rawId);
  float hsWidth(gemvalidation::cscHalfStripWidth(cscId));
  
  // map the pattern number to min int bending, max int bending
  int maxBendList[CSCConstants::NUM_CLCT_PATTERNS] = {
    0, 11, 11, -11, 9, -9, 7, -7, 5, -5, 3};
  int minBendList[CSCConstants::NUM_CLCT_PATTERNS] = {
    0, 0,   7, -7,  5, -5, 3, -3, 0, 0, -3};
  
  int averageBend = 0.5 * (minBendList[lct.getPattern()] + maxBendList[lct.getPattern()] );
                   
  // return the average bending
  return averageBend * hsWidth;
}
