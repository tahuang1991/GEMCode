#include "GEMCode/GEMValidation/interface/L1BaseMatcher.h"

L1BaseMatcher::L1BaseMatcher(SimHitMatcher& sh)
: BaseMatcher(csc.trk(), csc.vtx(), csc.conf(), csc.event(), csc.eventSetup())
{
  CSCTFSPset_ = conf().getParameter<edm::ParameterSet>("sectorProcessor");
  ptLUTset_ = CSCTFSPset_.getParameter<edm::ParameterSet>("PTLUT");
  srLUTset_ = CSCTFSPset_.getParameter<edm::ParameterSet>("SRLUT");
  
  hasMuScales_ = true;
  hasMuPtScale_ = true;

  ptLUT_ = nullptr;
  for(int e=0; e<2; e++) for (int s=0; s<6; s++) my_SPs_[e][s] = nullptr;
  //  CSCSectorReceiverLUT* srLUTs_[5][6][2];
  dtrc_ = nullptr;

  //std::cout<<" L1BaseMatcher constructor" <<std::endl;
  clear();
  init();
}

L1BaseMatcher::~L1BaseMatcher() 
{
 // if(ptLUT_) delete ptLUT_;
 // ptLUT_ = nullptr;
    //std::cout<<" L1BaseMatcher destructor" <<std::endl;
  delete ptLUT_;
  ptLUT_ = nullptr;

  for(int e=0; e<2; e++) for (int s=0; s<6; s++){
   // if  (my_SPs_[e][s])
       	delete my_SPs_[e][s];
        my_SPs_[e][s] = NULL;
    
  }
  for(int endcap = 1; endcap<=2; endcap++) {
    for(int sector=1; sector<=6; sector++) {
      for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
        if(station==1) for(int subSector=0; subSector<2; subSector++)
          delete srLUTs_[fpga++][sector-1][endcap-1];
        else
          delete srLUTs_[fpga++][sector-1][endcap-1];
      }
    }
  }


  if(dtrc_) 
    delete dtrc_;
  dtrc_ = nullptr;
  clear();
}

void 
L1BaseMatcher::clear()
{
}

void 
L1BaseMatcher::init()
{  
  try {
    eventSetup().get<L1MuTriggerScalesRcd>().get(muScalesHd_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerScalesRcd>& e) {
    hasMuScales_ = false;
    LogDebug("L1BaseMatcher") << "+++ Info: L1MuTriggerScalesRcd is unavailable. +++\n";
  }

  try {
    eventSetup().get<L1MuTriggerPtScaleRcd>().get(muPtScaleHd_);
  } catch (edm::eventsetup::NoProxyException<L1MuTriggerPtScaleRcd>& e) {
    hasMuPtScale_ = false;
    LogDebug("L1BaseMatcher") << "+++ Info: L1MuTriggerPtScaleRcd is unavailable. +++\n";
  }

  if (ptLUT_) delete ptLUT_;
  ptLUT_ = new CSCTFPtLUT(ptLUTset_, muScalesHd_.product(), muPtScaleHd_.product());

  /*// Sector Processor and Receiver
  for(int e=0; e<2; e++) {
    for (int s=0; s<6; s++){
      my_SPs_[e][s] = new CSCTFSectorProcessor(e+1, s+1, CSCTFSPset_, true, muScalesHd_.product(), muPtScaleHd_.product());
      my_SPs_[e][s]->initialize(eventSetup());
    }
  }*/

 const bool TMB07(true);  
  for(int endcap = 1; endcap<=2; endcap++) {
    for(int sector=1; sector<=6; sector++) {
      for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
        if(station==1) for(int subSector=0; subSector<2; subSector++)
          srLUTs_[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, subSector+1, station, srLUTset_, TMB07);
        else
          srLUTs_[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, 0, station, srLUTset_, TMB07);
      }
    }
  }
  dtrc_ = new CSCTFDTReceiver();
}
