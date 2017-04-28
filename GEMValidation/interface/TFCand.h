#ifndef GEMCode_GEMValidation_TFCand_h
#define GEMCode_GEMValidation_TFCand_h

//#include "GEMCode/GEMValidation/interface/PtassignmentHelper.h"
#include "GEMCode/GEMValidation/interface/TFTrack.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

class TFCand
{
 public:
  /// constructor
  TFCand(const L1MuRegionalCand* t);
  TFCand(const l1t::RegionalMuonCand* t);
  /// copy constructor
  TFCand(const TFCand&);
  /// destructor
  ~TFCand();

  void init(CSCTFPtLUT* ptLUT, 
	    edm::ESHandle< L1MuTriggerScales > &muScales, 
 	    edm::ESHandle< L1MuTriggerPtScale > &muPtScale); 

  void setDR(double);
  void setGlobalPhi(double x) { phi_ = x ; }
  void setBx(int x ) { bx_= x;}
  void setMatchedTFTrack(TFTrack* trk)  { nTFStubs = trk->nStubs(); tftrack_ = trk; }
  void print();

  const L1MuRegionalCand * l1Cand() const {return l1Cand_;}

  TFTrack* tftrack() const {return tftrack_;}

  std::vector < CSCDetId > ids() const {return ids_;}

  double pt() const {return pt_;}
  double eta() const {return eta_;}
  double phi() const {return phi_;}
  double phi_local() const {return phi_local_;}
  int quality() const {return quality_;}
  int charge() const {return charge_;}
  int bx() const {return bx_;}
  double dr() const {return dr_;}
  int tracktype() const {return trackType_; }
  unsigned int nStubs() const { return  nTFStubs; }
  
 private:

  const L1MuRegionalCand* l1Cand_;
  const l1t::RegionalMuonCand* gmtCand;
  TFTrack* tftrack_;
  std::vector<CSCDetId> ids_;
  double pt_;
  double eta_;
  double phi_;
  double phi_local_;
  double dr_;
  int quality_;
  int charge_;
  int bx_;
  int trackType_; 
  unsigned int nTFStubs;
};

#endif
