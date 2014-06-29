#include "GEMCode/GEMValidation/src/TrackMatcher.h"

TrackMatcher::TrackMatcher(SimHitMatcher& sh, CSCDigiMatcher& csc_dg, 
                           GEMDigiMatcher& gem_dg, RPCDigiMatcher& rpc_dg, 
                           CSCStubMatcher& csc_st)
: CSCStubMatcher(sh, csc_dg, gem_dg, rpc_dg)
, sh_matcher_(&sh)
, gem_digi_matcher_(&gem_dg)
, csc_digi_matcher_(&csc_dg)
, rpc_digi_matcher_(&rpc_dg)                 
{
  auto tfTrack = conf().getParameter<edm::ParameterSet>("cscTfTrack");
  cscTfTrackInputLabel_ = tfTrack.getParameter<edm::InputTag>("input");
  minBXTFTrack_ = tfTrack.getParameter<int>("minBX");
  maxBXTFTrack_ = tfTrack.getParameter<int>("minBX");
  verboseTFTrack_ = tfTrack.getParameter<int>("verbose");
  
  auto tfCand = conf().getParameter<edm::ParameterSet>("cscTfCand");
  cscTfCandInputLabel_ = tfCand.getParameter<edm::InputTag>("input");
  minBXTFCand_ = tfCand.getParameter<int>("minBX");
  maxBXTFCand_ = tfCand.getParameter<int>("minBX");
  verboseTFCand_ = tfCand.getParameter<int>("verbose");
  
  auto gmtRegCand = conf().getParameter<edm::ParameterSet>("gmtRegCand");
  gmtRegCandInputLabel_ = gmtRegCand.getParameter<edm::InputTag>("input");
  minBXGMTRegCand_ = gmtRegCand.getParameter<int>("minBX");
  maxBXGMTRegCand_ = gmtRegCand.getParameter<int>("minBX");
  verboseGMTRegCand_ = gmtRegCand.getParameter<int>("verbose");
  
  auto gmtCand = conf().getParameter<edm::ParameterSet>("gmtCand");
  gmtCandInputLabel_ = gmtCand.getParameter<edm::InputTag>("input");
  minBXGMTCand_ = gmtCand.getParameter<int>("minBX");
  maxBXGMTCand_ = gmtCand.getParameter<int>("minBX");
  verboseGMTCand_ = gmtCand.getParameter<int>("verbose");
  
  auto l1Extra = conf().getParameter<edm::ParameterSet>("l1Extra");
  l1ExtraInputLabel_ = l1Extra.getParameter<edm::InputTag>("input");
  minBXL1Extra_ = l1Extra.getParameter<int>("minBX");
  maxBXL1Extra_ = l1Extra.getParameter<int>("minBX");
  verboseL1Extra_ = l1Extra.getParameter<int>("verbose");

  CSCTFSPset_ = conf().getParameter<edm::ParameterSet>("sectorProcessor");
  ptLUTset_ = CSCTFSPset_.getParameter<edm::ParameterSet>("PTLUT");
  
  clear();
  init();
}

TrackMatcher::~TrackMatcher() 
{
}

void 
TrackMatcher::clear()
{
  tfTracks_.clear();
  tfCands_.clear();
  gmtRegCands_.clear();
  gmtCands_.clear();
  l1Extras_.clear();
}

void 
TrackMatcher::init()
{  
  // Trigger scales
  eventSetup().get<L1MuTriggerScalesRcd>().get(muScales_);
  eventSetup().get<L1MuTriggerPtScaleRcd>().get(muPtScale_);
  
  if (ptLUT_) delete ptLUT_;  
  ptLUT_ = new CSCTFPtLUT(ptLUTset_, muScales_.product(), muPtScale_.product());
  
  for(int e=0; e<2; e++) {
    for (int s=0; s<6; s++){
      if  (my_SPs_[e][s]) delete my_SPs_[e][s];
      my_SPs_[e][s] = new CSCTFSectorProcessor(e+1, s+1, CSCTFSPset_, true, muScales_.product(), muPtScale_.product());
      my_SPs_[e][s]->initialize(eventSetup());
    }
  }

  // tracks produced by TF
  edm::Handle<L1CSCTrackCollection> hl1Tracks;
  event().getByLabel(cscTfTrackInputLabel_,hl1Tracks);
  matchTfTrackToSimTrack(*hl1Tracks.product());
  
  // L1 muon candidates after CSC sorter
  edm::Handle<std::vector<L1MuRegionalCand> > hl1TfCands;
  event().getByLabel(cscTfCandInputLabel_, hl1TfCands);
  matchTfCandToSimTrack(*hl1TfCands.product());
}

void 
TrackMatcher::matchTfTrackToSimTrack(const L1CSCTrackCollection& tracks)
{
  for (L1CSCTrackCollection::const_iterator trk = tracks.begin(); trk != tracks.end(); trk++) {
    TFTrack track(&trk->first);
    track.init(ptLUT_, muScales_, muPtScale_);
    // calculate the DR
    track.setDR(this->trk());

    // debugging
    if (verboseTFTrack_){
      std::cout << "L1CSC TFTrack "<< trk-tracks.begin() << " information:" << std::endl;
      std::cout << "pt (GeV/c) = " << track.pt() << ", eta = " << track.eta() 
                << ", phi = " << track.phi() << ", dR(sim,L1) = " << track.dr() << std::endl;      
    }

    // check the stubs
    if (verboseTFTrack_){
      std::cout << "Stubs:" << std::endl;
    }
    for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator detUnitIt = trk->second.begin(); 
         detUnitIt != trk->second.end(); detUnitIt++) {
      const CSCDetId& id = (*detUnitIt).first;
      if (verboseTFTrack_){
        std::cout << "DetId: " << id << std::endl;
      }
      const CSCCorrelatedLCTDigiCollection::Range& range = (*detUnitIt).second;
      for (CSCCorrelatedLCTDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; digiIt++) {
        auto lct(*digiIt);
        // check for valid stubs
        if (!(lct.isValid())) continue;
        track.addTriggerDigi(&lct);
        track.addTriggerDigiId(id);
        track.addTriggerEtaPhi(intersectionEtaPhi(id, lct.getKeyWG(), lct.getStrip()));
        track.addTriggerStub(buildTrackStub(lct, id));
        if (verboseTFTrack_){
          std::cout << "LCT" << digiIt-range.first << ": " << lct << std::endl;
        }
      }
    }
  }
}

void 
TrackMatcher::matchTfCandToSimTrack(const std::vector<L1MuRegionalCand>& tracks)
{
}

void 
TrackMatcher::matchGmtRegCandToSimTrack(const L1MuRegionalCand& tracks)
{
}

void 
TrackMatcher::matchGmtCandToSimTrack(const L1MuGMTExtendedCand& tracks)
{
}

TFTrack* 
TrackMatcher::bestTFTrack(bool sortPtFirst)
{
  return 0;
  /*
  if (tfTracks_.size()==0) return NULL;
  
  // determine max # of matched stubs in the TFTrack collection
  int maxNMatchedMPC = 0;
  for (unsigned i=0; i<tfTracks_.size(); i++) {
    int nst=0;
    for (size_t s=0; s<tfTracks_.at(i).ids.size(); s++) 
      if (tfTracks_.at(i).mplcts[s]->deltaOk) nst++;
    if (nst>maxNMatchedMPC) maxNMatchedMPC = nst;
  }
  // collect tracks with max # of matched stubs
  std::vector<unsigned> bestMatchedTracks;
  for (unsigned i=0; i<tfTracks_.size(); i++) {
    int nst=0;
    for (size_t s=0; s<tfTracks_.at(i).ids.size(); s++) 
      if (tfTracks_.at(i).mplcts[s]->deltaOk) nst++;
    if (nst==maxNMatchedMPC) bestMatchedTracks.push_back(i);
  }
  
  // already found the best TFTrack
  if (bestMatchedTracks.size()==1) return &(tfTracks_[bestMatchedTracks[0]]);
  
  // case when you have more than 1 best TFTrack
  // first sort by quality
  double qBase  = 1000000.;
  // then sort by Pt inside the cone (if sortPtFirst), then sort by DR
  double ptBase = 0.;
  if (sortPtFirst) ptBase = 1000.;
  unsigned maxI = 99;
  double maxRank = -999999.;
  for (unsigned i=0; i<tfTracks_.size(); i++) {
    if (bestMatchedTracks.size()) {
      bool gotit=0;
      for (unsigned m=0;m<bestMatchedTracks.size();m++) if (bestMatchedTracks[m]==i) gotit=1;
      if (!gotit) continue;
    }
    double rr = qBase*tfTracks_.at(i).q_packed + ptBase*tfTracks_.at(i).pt_packed + 1./(0.01 + tfTracks_.at(i).dr);
    if (rr > maxRank) { maxRank = rr; maxI = i;}
  }
  if (maxI==99) return NULL;
  return &(tfTracks_.at(maxI));
  */
}


TFCand* 
TrackMatcher::bestTFCand(bool sortPtFirst)
{
  /*
    if (cands.size()==0) return NULL;

    // determine max # of matched stubs
    int maxNMatchedMPC = 0;
    for (unsigned i=0; i<cands.size(); i++) 
    {
        int nst=0;
        if (cands[i].tftrack==0) continue;
        for (size_t s=0; s<cands[i].tftrack->ids.size(); s++) 
            if (cands[i].tftrack->mplcts[s]->deltaOk) nst++;
        if (nst>maxNMatchedMPC) maxNMatchedMPC = nst;
    }

    // collect tracks with max # of matched stubs
    std::vector<unsigned> bestMatchedTracks;
    if (maxNMatchedMPC>0) {
        for (unsigned i=0; i<cands.size(); i++) 
        {
            int nst=0;
            if (cands[i].tftrack==0) continue;
            for (size_t s=0; s<cands[i].tftrack->ids.size(); s++) 
                if (cands[i].tftrack->mplcts[s]->deltaOk) nst++;
            if (nst==maxNMatchedMPC) bestMatchedTracks.push_back(i);
        }
        if (bestMatchedTracks.size()==1) return &(cands[bestMatchedTracks[0]]);
    }

    // first sort by quality
    double qBase  = 1000000.;
    // then sort by Pt inside the cone (if sortPtFirst), then sort by DR
    double ptBase = 0.;
    if (sortPtFirst) ptBase = 1000.;
    unsigned maxI = 99;
    double maxRank = -999999.;
    for (unsigned i=0; i<cands.size(); i++) 
    {
        if (bestMatchedTracks.size()) {
            bool gotit=0;
            for (unsigned m=0;m<bestMatchedTracks.size();m++) if (bestMatchedTracks[m]==i) gotit=1;
            if (!gotit) continue;
        }
        // quality criterium you apply to get the best TFCand
        double rr = qBase*cands[i].l1cand->quality_packed() + ptBase*cands[i].l1cand->pt_packed() + 1./(0.01 + cands[i].dr);
        if (rr > maxRank) { maxRank = rr; maxI = i;}
    }
OB    if (maxI==99) return NULL;
    return &(cands[maxI]);
  */
  return 0;
}

GMTRegCand* 
TrackMatcher::bestGMTRegCand(bool sortPtFirst)
{
  // first sort by Pt inside the cone (if sortPtFirst), then sort by DR
  if (gmtRegCands_.size()==0) return nullptr;
  const double ptBase(sortPtFirst ? 1000. : 0.);
  unsigned maxI = 99;
  double maxRank = -999999.;
  for (unsigned i=0; i<gmtRegCands_.size(); i++) {
    // quality criterium to sort the GMT Regional candidates
    const double rank(ptBase*gmtRegCands_.at(i)->pt() + 1./(0.01 + gmtRegCands_.at(i)->dr()));
    if (rank > maxRank) { 
      maxRank = rank; 
      maxI = i;
    }
  }
  if (maxI==99) return nullptr;
  return gmtRegCands_.at(maxI);
}

GMTCand* 
TrackMatcher::bestGMTCand(bool sortPtFirst)
{
  // first sort by Pt inside the cone (if sortPtFirst), then sort by DR
  if (gmtCands_.size()==0) return nullptr;
  const double ptBase(sortPtFirst ? 1000. : 0.);
  unsigned maxI = 99;
  double maxRank = -999999.;
  for (unsigned i=0; i<gmtCands_.size(); i++) {
    // quality criterium to sort the GMT candidates
    const double rank(ptBase*gmtCands_.at(i)->pt() + 1./(0.01 + gmtCands_.at(i)->dr()));
    if (rank > maxRank) { 
      maxRank = rank; 
      maxI = i;
    }
  }
  if (maxI==99) return nullptr;
  return gmtCands_.at(maxI);
}

L1Extra* 
TrackMatcher::bestL1Extra(bool sortPtFirst)
{
  return 0;
}

csctf::TrackStub 
TrackMatcher::buildTrackStub(const CSCCorrelatedLCTDigi &d, CSCDetId id)
{
  const unsigned fpga((id.station() == 1) ? CSCTriggerNumbering::triggerSubSectorFromLabels(id) - 1 : id.station());
  CSCSectorReceiverLUT* srLUT = srLUTs_[fpga][id.triggerSector()-1][id.endcap()-1];
  const unsigned cscid(CSCTriggerNumbering::triggerCscIdFromLabels(id));
  const unsigned cscid_special((id.station()==1 && id.ring()==4) ? cscid + 9 : cscid);
  const lclphidat lclPhi(srLUT->localPhi(d.getStrip(), d.getPattern(), d.getQuality(), d.getBend()));
  const gblphidat gblPhi(srLUT->globalPhiME(lclPhi.phi_local, d.getKeyWG(), cscid_special));
  const gbletadat gblEta(srLUT->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, d.getKeyWG(), cscid));

  return csctf::TrackStub(d, id, gblPhi.global_phi, gblEta.global_eta);
}



std::pair<float, float> 
TrackMatcher::intersectionEtaPhi(CSCDetId id, int wg, int hs)
{
  const CSCDetId layerId(id.endcap(), id.station(), id.ring(), id.chamber(), CSCConstants::KEY_CLCT_LAYER);
  const CSCLayer* csclayer(cscGeometry_->layer(layerId));
  const CSCLayerGeometry* layer_geo(csclayer->geometry());

  // LCT::getKeyWG() starts from 0
  const float wire(layer_geo->middleWireOfGroup(wg + 1));

  // half-strip to strip
  // note that LCT's HS starts from 0, but in geometry strips start from 1
  const float fractional_strip(0.5 * (hs + 1) - 0.25);
  const LocalPoint csc_intersect(layer_geo->intersectionOfStripAndWire(fractional_strip, wire));
  const GlobalPoint csc_gp(cscGeometry_->idToDet(layerId)->surface().toGlobal(csc_intersect));

  return std::make_pair(csc_gp.eta(), csc_gp.phi());
}
