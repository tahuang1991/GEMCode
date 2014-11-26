#include "GEMCode/GEMValidation/src/CSCStubMatcher.h"
#include "GEMCode/GEMValidation/src/SimHitMatcher.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include <DataFormats/MuonDetId/interface/CSCTriggerNumbering.h>

#include <algorithm>

using namespace std;
using namespace matching;


CSCStubMatcher::CSCStubMatcher(SimHitMatcher& sh, CSCDigiMatcher& dg, GEMDigiMatcher& gem_dg, RPCDigiMatcher& rpc_dg)
: DigiMatcher(sh)
, digi_matcher_(&dg)
, gem_digi_matcher_(&gem_dg)
, rpc_digi_matcher_(&rpc_dg)
, sh_matcher_(&sh)
{
  auto cscCLCT_ = conf().getParameter<edm::ParameterSet>("cscCLCT");
  clctInput_ = cscCLCT_.getParameter<edm::InputTag>("input");
  minBXCLCT_ = cscCLCT_.getParameter<int>("minBX");
  maxBXCLCT_ = cscCLCT_.getParameter<int>("maxBX");
  verboseCLCT_ = cscCLCT_.getParameter<int>("verbose");
  minNHitsChamberCLCT_ = cscCLCT_.getParameter<int>("minNHitsChamber");
  runCLCT_ = cscCLCT_.getParameter<bool>("run");

  auto cscALCT_ = conf().getParameter<edm::ParameterSet>("cscALCT");
  alctInput_ = cscALCT_.getParameter<edm::InputTag>("input");
  minBXALCT_ = cscALCT_.getParameter<int>("minBX");
  maxBXALCT_ = cscALCT_.getParameter<int>("maxBX");
  verboseALCT_ = cscALCT_.getParameter<int>("verbose");
  minNHitsChamberALCT_ = cscALCT_.getParameter<int>("minNHitsChamber");
  runALCT_ = cscCLCT_.getParameter<bool>("run");

  auto cscLCT_ = conf().getParameter<edm::ParameterSet>("cscLCT");
  lctInput_ = cscLCT_.getParameter<edm::InputTag>("input");
  minBXLCT_ = cscLCT_.getParameter<int>("minBX");
  maxBXLCT_ = cscLCT_.getParameter<int>("maxBX");
  verboseLCT_ = cscLCT_.getParameter<int>("verbose");
  minNHitsChamberLCT_ = cscLCT_.getParameter<int>("minNHitsChamber");
  addGhostLCTs_ = cscLCT_.getParameter<bool>("addGhosts");
  matchAlctGem_ = cscLCT_.getParameter<bool>("matchAlctGem");
  matchClctGem_ = cscLCT_.getParameter<bool>("matchClctGem");
  matchAlctRpc_ = cscLCT_.getParameter<bool>("matchAlctRpc");
  matchClctRpc_ = cscLCT_.getParameter<bool>("matchClctRpc");
  hsFromSimHitMean_ = cscLCT_.getParameter<bool>("hsFromSimHitMean");
  runLCT_ = cscLCT_.getParameter<bool>("run");

  auto cscMPLCT_ = conf().getParameter<edm::ParameterSet>("cscMPLCT");
  mplctInput_ = cscMPLCT_.getParameter<edm::InputTag>("input");
  minBXMPLCT_ = cscMPLCT_.getParameter<int>("minBX");
  maxBXMPLCT_ = cscMPLCT_.getParameter<int>("maxBX");
  verboseMPLCT_ = cscMPLCT_.getParameter<int>("verbose");
  minNHitsChamberMPLCT_ = cscMPLCT_.getParameter<int>("minNHitsChamber");
  addGhostMPLCTs_ = cscMPLCT_.getParameter<bool>("addGhosts");
  runMPLCT_ = cscMPLCT_.getParameter<bool>("run");

  minNHitsChamber_ = conf().getUntrackedParameter<int>("minNHitsChamber", 4);

  setVerbose(conf().getUntrackedParameter<int>("verboseCSCStub", 0));

  if (! ( clctInput_.label().empty() || alctInput_.label().empty() ||
          lctInput_.label().empty() || mplctInput_.label().empty() )
      )
  {
   init();
  }
}


CSCStubMatcher::~CSCStubMatcher() {}


void CSCStubMatcher::init()
{
  edm::Handle<CSCCLCTDigiCollection> clcts;
  event().getByLabel(clctInput_, clcts);

  edm::Handle<CSCALCTDigiCollection> alcts;
  event().getByLabel(alctInput_, alcts);

  edm::Handle<CSCCorrelatedLCTDigiCollection> lcts;
  event().getByLabel(lctInput_, lcts);

  edm::Handle<CSCCorrelatedLCTDigiCollection> mplcts;
  event().getByLabel(mplctInput_, mplcts);

  if (runCLCT_) matchCLCTsToSimTrack(*clcts.product());
  if (runALCT_) matchALCTsToSimTrack(*alcts.product());
  if (runLCT_) matchLCTsToSimTrack(*lcts.product());
  if (runMPLCT_) matchMPLCTsToSimTrack(*mplcts.product());
}


void
CSCStubMatcher::matchCLCTsToSimTrack(const CSCCLCTDigiCollection& clcts)
{
  // only look for stub in chambers that have digis matching to this track

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

      // match by half-strip with the digis
      if (digi_strips.find(half_strip) == digi_strips.end())
      {
        if (verbose()) cout<<"clctBAD"<<endl;
        continue;
      }
      if (verbose()) cout<<"clctGOOD"<<endl;


      // store matching CLCTs in this chamber
      chamber_to_clcts_[id].push_back(mydigi);

      if (chamber_to_clct_.find(id) != chamber_to_clct_.end())
      {
        //cout<<"WARNING!!! there already was matching CLCT "<<chamber_to_clct_[id]<<endl;
        //cout<<"   new digi: "<<mydigi<<endl;

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

      if (verbose()) cout<<"alct "<<ch_id<<" "<<*a<<endl;

      // check that the BX for stub wasn't too early or too late
      if (a->getBX() < minBXALCT_ || a->getBX() > maxBXALCT_) continue;

      int wg = a->getKeyWG() + 1; // as ALCT wiregroups numbers start from 0
      auto mydigi = make_digi(id, wg, a->getBX(), CSC_ALCT, a->getQuality());

      // store all ALCTs in this chamber
      chamber_to_alcts_all_[id].push_back(mydigi);

      // match by wiregroup with the digis
      if (digi_wgs.find(wg) == digi_wgs.end())
      {
        if (verbose()) cout<<"alctBAD"<<endl;
        continue;
      }
      if (verbose()) cout<<"alctGOOD"<<endl;

      // store matching ALCTs in this chamber
      chamber_to_alcts_[id].push_back(mydigi);

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
    map<int, DigiContainer> bx_to_lcts;
    for (auto lct = lcts_in_det.first; lct != lcts_in_det.second; ++lct)
    {
      if (!lct->isValid()) continue;

      if (verbose()) cout<<"lct in detId "<<ch_id<<" "<<*lct<<endl;

      int bx = lct->getBX();

      // check that the BX for stub wasn't too early or too late
      if (bx < minBXLCT_ || bx > maxBXLCT_) continue;

      int hs = lct->getStrip() + 1; // LCT halfstrip and wiregoup numbers start from 0
      int wg = lct->getKeyWG() + 1;

      float dphi = lct->getGEMDPhi();

      auto mydigi = make_digi(id, hs, bx, CSC_LCT, lct->getQuality(), lct->getPattern(), wg, dphi);
      lcts_tmp.push_back(mydigi);
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

          auto lct21 = lct22;
          digi_wg(lct21) = wg1;
          lcts_tmp.push_back(lct21);
          //cout<<"added ghosts"<<endl<<lct11<<"    "<<lct22<<endl <<lct12<<"    "<<lct21<<endl;
        }
      }
    } // lcts_in_det

    size_t n_lct = lcts_tmp.size();
    if (verbose()) cout<< "number of lcts = "<<n_lct<<endl;
    if (n_lct == 0) continue; // no LCTs in this chamber

    // assign the non necessarily matching LCTs
    chamber_to_lcts_all_[id] = lcts_tmp;

    if (verbose() and !(n_lct == 1 || n_lct == 2 || n_lct == 4 ) )
    {
      cout<<"WARNING!!! weird #LCTs="<<n_lct;
      for (auto &s: lcts_tmp) cout<<"  "<<s<<endl;
      //continue;
    }

    // find a matching LCT
    const GEMDetId gemDetId(GEMDetId(ch_id.zendcap(),1,ch_id.station(),1,ch_id.chamber(),0));
    // find matching rpc chamber (only valid for me31 and me41)
    const int csc_trig_sect(CSCTriggerNumbering::triggerSectorFromLabels(ch_id));
    const int csc_trig_id( CSCTriggerNumbering::triggerCscIdFromLabels(ch_id));
    const int csc_trig_chid((3*(csc_trig_sect-1)+csc_trig_id)%18 +1);
    const int rpc_trig_sect((csc_trig_chid-1)/3+1);
    const int rpc_trig_subsect((csc_trig_chid-1)%3+1);
    const RPCDetId rpcDetId(RPCDetId(ch_id.zendcap(),1,ch_id.station(),rpc_trig_sect,1,rpc_trig_subsect,0));

    auto clct(clctsInChamber(id));
    auto alct(alctsInChamber(id));
    auto pads(gem_digi_matcher_->coPadsInSuperChamber(gemDetId));
    auto rpcDigis(rpc_digi_matcher_->digisInChamber(rpcDetId));
    const auto hits = sh_matcher_->hitsInChamber(id);
    
    for (auto &lct: lcts_tmp)
     {  
	 //std::cout << " LCT " << lct << std::endl;
        for (unsigned int j=0; j<alct.size();j++){
        for (unsigned int i=0; i<clct.size()+1;i++){
           //if (i < clct.size()) 
	   //std::cout << " ALCT " << alct[j] << " CLCT " << clct[i] << std::endl;
	   //else std::cout << " NO CLCT case " << " ALCT " << alct[j] << std::endl;
           auto hasPad(pads.size()!=0);
           auto hasDigis(rpcDigis.size()!=0);

           const bool caseAlctClct(j<alct.size() && i<clct.size());
           const bool caseAlctGem(is_valid(alct[j]) and hasPad and i==clct.size() and (ch_id.station() == 1 or ch_id.station() == 2));
           //const bool caseClctGem(is_valid(clct[i]) and hasPad and !is_valid(alct[j]) and (ch_id.station() == 1 or ch_id.station() == 2));
           const bool caseAlctRpc(is_valid(alct[j]) and hasDigis and i==clct.size() and (ch_id.station() == 3 or ch_id.station() == 4));
           //const bool caseClctRpc(is_valid(clct[i]) and hasDigis and !is_valid(alct[j]) and (ch_id.station() == 3 or ch_id.station() == 4));

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


               if ( caseAlctClct and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and my_wg == digi_wg(lct))){
                    if (verbose()) cout<<"  BAD LCT in caseAlctClct"<<endl;
                    continue;
                }


               if (matchAlctGem_ and caseAlctGem and !(my_bx == digi_bx(lct) and std::abs(my_hs_gemrpc - digi_channel(lct))<3 and my_wg == digi_wg(lct) ) ){
                    if (verbose()) cout<<"  BAD LCT in matchAlctGem"<<endl;
                    continue;
               }
	       else if (caseAlctGem and !matchAlctGem_) continue;

              /* if (matchClctGem_ and caseClctGem and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and std::abs(my_wg_gemrpc_mean - digi_wg(lct))<3) ){
                    if (verbose()) cout<<"  BAD"<<endl;
                    continue;
               }*/

              if (matchAlctRpc_ and caseAlctRpc and !(my_bx == digi_bx(lct) and std::abs(my_hs_gemrpc - digi_channel(lct))<3 and my_wg == digi_wg(lct) ) ){
                   if (verbose()) cout<<"  BAD LCT in matchAlctRpc"<<endl;
                   continue;
              }
	      else if (caseAlctRpc and !matchAlctRpc_) continue;
             /* if (matchClctRpc_ and caseClctRpc and !(my_bx == digi_bx(lct) and my_hs == digi_channel(lct) and std::abs(my_wg_gemrpc_mean - digi_wg(lct))<3) ){
                    if (verbose()) cout<<"  BAD"<<endl;
                    continue;
              }*/


             if (chamber_to_lct_.find(id) != chamber_to_lct_.end()){
                       //  cout<<"ALARM!!! here already was matching LCT "<<chamber_to_lct_[id]<<endl;
                       //  cout<<"   new digi: "<<lct<<endl;
             }
             
              chamber_to_lct_[id] = lct;
              chamber_to_lcts_[id].push_back(lct);
              //std::cout << " this LCT is matched to simtrack!! " << std::endl;
              break;
            } //clct loop over
        }// alct loop over

    } // lct loop over

  }

  if (verbose() and n_minLayers > 0)
  {
    if (chamber_to_lct_.size() == 0)
    {
      cout<<"effNoLCT"<<endl;
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
    else cout<<"effYesLCT" << std::endl;
  }
}


void
CSCStubMatcher::matchMPLCTsToSimTrack(const CSCCorrelatedLCTDigiCollection& mplcts)
{
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

      float dphi = lct->getGEMDPhi();

      auto mydigi = make_digi(id, hs, bx, CSC_LCT, lct->getQuality(), lct->getPattern(), wg, dphi);
      mplcts_tmp.push_back(mydigi);
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


const matching::DigiContainer&
CSCStubMatcher::mplctsInChamber(unsigned int detid) const
{
  if (chamber_to_mplcts_.find(detid) == chamber_to_mplcts_.end()) return no_digis_;
  return chamber_to_mplcts_.at(detid);
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

  auto mydigi = make_digi(id.rawId(), lct.getStrip()+1, lct.getBX(), CSC_LCT, lct.getQuality(), lct.getPattern(),lct.getKeyWG()+1,lct.getGEMDPhi());
  try{
  auto alldigis(chamber_to_lcts_.at(id.rawId()));
//  std::cout << "in checkstub, mydigi  " << mydigi << std::endl;
  for (auto p : alldigis) 
  {   
 //     std::cout << " digi matched to simtrack " << p << std::endl;
      if (p==mydigi) return true;
   }
  }
  catch (exception& e)
    {
//	return false;
    }
  //std::cout << "   matching failed   " << std::endl;
 return false;
}

