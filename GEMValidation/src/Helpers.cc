#include "GEMCode/GEMValidation/interface/Helpers.h"

int 
gemvalidation::chamber(const DetId& id) 
{
  if (id.det() != DetId::Detector::Muon) return -99;
  int chamberN = 0;
  switch(id.subdetId()){
  case MuonSubdetId::GEM:
    chamberN = GEMDetId(id).chamber();
    break;
  case MuonSubdetId::RPC:
    // works only for endcap!!
    chamberN = RPCDetId(id).sector();
    break;
  case MuonSubdetId::CSC:
    chamberN = CSCDetId(id).chamber();
    break;
  case MuonSubdetId::ME0:
    chamberN = ME0DetId(id).chamber();
    break;
  };
  return chamberN;
}


unsigned int
gemvalidation::gemDetFromCSCDet(unsigned int id,int layer)
{
  CSCDetId cscId(id);
  // returns the gem superr chamber for a given ME1/1 chamber(ME1/1a + ME1/1b)
  GEMDetId gemId(cscId.zendcap(), 1, cscId.station(), layer, cscId.chamber(),0); 
  return gemId.rawId();
}

std::pair<unsigned int, unsigned int> 
gemvalidation::gemDetsFromCSCDet(unsigned int id)
{
  return std::make_pair(gemDetFromCSCDet(id,1),gemDetFromCSCDet(id,2));
}

// return MuonType for a particular DetId
int 
gemvalidation::toGEMType(GEMDetId id)
{
  const int re(id.region());
  const int st(id.station());
  const int ri(id.ring());
  // endcap
  if (abs(re)==1) {
    if (st ==1) {
      if (ri==1) return GEM_ME11;
    }
    else if (st ==2) {
      if (ri==1) return GEM_ME21;
    }
  }
  return GEM_ALL;
}

int 
gemvalidation::toRPCType(RPCDetId id)
{
  const int re(id.region());
  const int st(id.station());
  const int ri(id.ring());
  // endcap
  if (abs(re)==1) {
    if (ri ==1) {
      if (st==2) return RPC_ME12;
      if (st==3) return RPC_ME13;
    }
    else if (ri ==2) {
      if (st==2) return RPC_ME22;
      if (st==3) return RPC_ME23;
    }
    else if (ri ==3) {
      if (st==1) return RPC_ME31;
      if (st==2) return RPC_ME32;
      if (st==3) return RPC_ME33;
    }
    else if (ri ==4) {
      if (st==1) return RPC_ME41;
      if (st==2) return RPC_ME42;
      if (st==3) return RPC_ME43;
    }
  }
  // Barrel
  else {
    if (ri==-2) {
      if (st==1) return RPC_MB21n;
      if (st==2) return RPC_MB22n;
      if (st==3) return RPC_MB23n;
      if (st==4) return RPC_MB24n;
    }
    else if (ri==-1) {
      if (st==1) return RPC_MB11n;
      if (st==2) return RPC_MB12n;
      if (st==3) return RPC_MB13n;
      if (st==4) return RPC_MB14n;
    }
    else if (ri==0) {
      if (st==1) return RPC_MB01;
      if (st==2) return RPC_MB02;
      if (st==3) return RPC_MB03;
      if (st==4) return RPC_MB04;
    }
    else if (ri==1) {
      if (st==1) return RPC_MB11p;
      if (st==2) return RPC_MB12p;
      if (st==3) return RPC_MB13p;
      if (st==4) return RPC_MB14p;
    }
    else if (ri==2) {
      if (st==1) return RPC_MB21p;
      if (st==2) return RPC_MB22p;
      if (st==3) return RPC_MB23p;
      if (st==4) return RPC_MB24p;
    }
  }
  return RPC_ALL;
}

int 
gemvalidation::toDTType(DTChamberId id)
{
  const int wh(id.wheel());
  const int st(id.station());
  if (wh==-2) {
    if (st==1) return DT_MB21n;
    if (st==2) return DT_MB22n;
    if (st==3) return DT_MB23n;
    if (st==4) return DT_MB24n;
  }
  if (wh==-1) {
    if (st==1) return DT_MB11n;
    if (st==2) return DT_MB12n;
    if (st==3) return DT_MB13n;
    if (st==4) return DT_MB14n;
  }
  if (wh==0) {
    if (st==1) return DT_MB01;
    if (st==2) return DT_MB02;
    if (st==3) return DT_MB03;
    if (st==4) return DT_MB04;
  }
  if (wh==1) {
    if (st==1) return DT_MB11p;
    if (st==2) return DT_MB12p;
    if (st==3) return DT_MB13p;
    if (st==4) return DT_MB14p;
  }
  if (wh==2) {
    if (st==1) return DT_MB21p;
    if (st==2) return DT_MB22p;
    if (st==3) return DT_MB23p;
    if (st==4) return DT_MB24p;
  }
  return DT_ALL;
}

int 
gemvalidation::toDTType(DTWireId id)
{
  return toDTType(DTChamberId(id.rawId()));
  //  return toDTType(DTChamberId(id.wheel(),id.station(),id.sector()));
}

