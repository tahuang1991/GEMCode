#ifndef GEMCode_GEMValidation_Helpers_h
#define GEMCode_GEMValidation_Helpers_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"

static const float AVERAGE_GEM_Z(568.6); // [cm]

static const float AVERAGE_GE11_ODD_Z(568.6); // [cm]
static const float AVERAGE_GE11_EVEN_Z(568.6); // [cm]

static const float AVERAGE_GE21_LONG_Z(568.6); // [cm]
static const float AVERAGE_GE21_SHORT_Z(568.6); // [cm]

static const float AVERAGE_ME11_EVEN_Z(585); // [cm]
static const float AVERAGE_ME11_ODD_Z(615); // [cm]

static const float AVERAGE_ME21_EVEN_Z(820); // [cm]
static const float AVERAGE_ME21_ODD_Z(835); // [cm]

static const float AVERAGE_ME0_Z(568.6); // [cm]

/// Muon Subsystem
enum MuonType {DT= 1, CSC=2, RPC=3, GEM=4, ME0=5};

/// CSC chamber types, according to CSCDetId::iChamberType()
enum CSCType {CSC_ALL = 0, CSC_ME11,
	      CSC_ME1a, CSC_ME1b, CSC_ME12, CSC_ME13, CSC_ME21, 
	      CSC_ME22, CSC_ME31, CSC_ME32, CSC_ME41, CSC_ME42};

/// GEM chamber types
enum GEMType {GEM_ALL = 0, 
	      GEM_ME11, GEM_ME21};

/// RPC endcap chamber types -- FIXME
enum RPCType {RPC_ALL = 0, 
	      RPC_ME12, RPC_ME13, RPC_ME22, RPC_ME23, RPC_ME31, 
	      RPC_ME32, RPC_ME33, RPC_ME41, RPC_ME42, RPC_ME43,
	      RPC_MB01,  RPC_MB02,  RPC_MB03,  RPC_MB04, 
	      RPC_MB11p, RPC_MB12p, RPC_MB13p, RPC_MB14p, 
	      RPC_MB21p, RPC_MB22p, RPC_MB23p, RPC_MB24p, 
	      RPC_MB11n, RPC_MB12n, RPC_MB13n, RPC_MB14n, 
	      RPC_MB21n, RPC_MB22n, RPC_MB23n, RPC_MB24n};

/// DT chamber types -- FIXME
enum DTType {DT_ALL = 0, 
	     DT_MB01,  DT_MB02,  DT_MB03,  DT_MB04, 
	     DT_MB11p, DT_MB12p, DT_MB13p, DT_MB14p, 
	     DT_MB21p, DT_MB22p, DT_MB23p, DT_MB24p, 
	     DT_MB11n, DT_MB12n, DT_MB13n, DT_MB14n, 
	     DT_MB21n, DT_MB22n, DT_MB23n, DT_MB24n};

const double ME11GEMdPhi[9][3] = {
  {-2 , 1.0, 1.0 },
  {3 , 0.03971647, 0.01710244 },
  {5 , 0.02123785, 0.00928431 },
  {7 , 0.01475524, 0.00650928 },
  {10, 0.01023299, 0.00458796 },
  {15, 0.00689220, 0.00331313 },
  {20, 0.00535176, 0.00276152 },
  {30, 0.00389050, 0.00224959 },
  {40, 0.00329539, 0.00204670 }
};
const double ME21GEMdPhi[9][3] = {
  {-2 , 1.0, 1.0 },
  {3 , 0.01832829, 0.01003643 },
  {5 , 0.01095490, 0.00631625 },
  {7 , 0.00786026, 0.00501017 },
  {10, 0.00596349, 0.00414560 },
  {15, 0.00462411, 0.00365550 },
  {20, 0.00435298, 0.00361550 },
  {30, 0.00465160, 0.00335700 },
  {40, 0.00372145, 0.00366262 }
};

namespace gemvalidation
{
inline bool 
is_dt(unsigned int detId) 
{
  return (DetId(detId)).subdetId() == MuonSubdetId::DT;
}

inline bool 
is_gem(unsigned int detId) 
{
  return (DetId(detId)).subdetId() == MuonSubdetId::GEM;
}

inline bool 
is_csc(unsigned int detId) 
{
  return (DetId(detId)).subdetId() == MuonSubdetId::CSC;
}

inline bool 
is_rpc(unsigned int detId) 
{
  return (DetId(detId)).subdetId() == MuonSubdetId::RPC;
}

inline bool 
is_me0(unsigned int detId) 
{
  return (DetId(detId)).subdetId() == MuonSubdetId::ME0;
}

int 
chamber(const DetId& id);

unsigned int
gemDetFromCSCDet(unsigned int id,int layer);

std::pair<unsigned int, unsigned int> 
gemDetsFromCSCDet(unsigned int id);

// return MuonType for a particular DetId
int toGEMType(int st, int ri);
int toRPCType(int re, int st, int ri);
int toDTType(int wh, int st);
int toCSCType(int st, int ri);

std::string toGEMTypeString(int st, int ri);
std::string toRPCTypeString(int re, int st, int ri);
std::string toDTTypeString(int wh, int st);
std::string toCSCTypeString(int st, int ri);

template<typename PROD>
bool
getByLabel(std::vector<edm::InputTag> const& tags, edm::Handle<PROD>& result, const edm::Event& iEvent)
{
  const bool verbose(false);
  bool inputTagIsNotValid(true);
  for (unsigned i=0; i<tags.size(); ++i){
    iEvent.getByLabel(tags[i], result);
    if (result.isValid()) {
      if (verbose) std::cout << tags[i] << " is a valid inputTag " << i << std::endl;
      inputTagIsNotValid = false;
      break;
    } else {
      if (verbose) std::cout << tags[i] << " is an invalid inputTag " << i << std::endl;
    }
  }
  return (!inputTagIsNotValid);
}
}

#endif
