#include "GEMCode/GEMValidation/interface/BarrelTriggerPtAssignmentHelper.h"

float getDirectionBasedPt2Stubs(float DPhi, int DT_type)
{
  std::string dt_type = DT_Types_2stubs_string[DT_type];
  return getDirectionBasedPt2Stubs(DPhi, dt_type);
}

float getDirectionBasedPt2Stubs(float DPhi, std::string DT_type)
{
  std::vector<float> pts = DirectionBasedPtLUT_2stubs[DT_type];
  std::vector<float> dphis = DirectionBasedDeltaPhiLUT_2stubs[DT_type];

  // very low pT
  if (DPhi > dphis[0]) return 3.5;

  // very high pT
  if (DPhi < dphis[-1]) return 55;

  // all other cases
  for (unsigned int i=0; i < dphis.size(); ++i){
    if (DPhi > dphis[i]){
      return pts[i-1];
    }
  }
  return 0;
}

float getDirectionBasedPt3Stubs(float DPhi1, float DPhi2, int DT_type)
{
  const std::string dt_type = DT_Types_3stubs_string[DT_type];
  return getDirectionBasedPt3Stubs(DPhi1, DPhi2, dt_type);
}

float getDirectionBasedPt3Stubs(float DPhi1, float DPhi2, std::string DT_type)
{
  return 0;
}

float getDirectionBasedPt4Stubs(float DPhi1, float DPhi2, int DT_type)
{
  const std::string dt_type = DT_Types_4stubs_string[DT_type];
  return getDirectionBasedPt4Stubs(DPhi1, DPhi2, dt_type);
}

float getDirectionBasedPt4Stubs(float DPhi1, float DPhi2, std::string DT_type)
{
  return 0;
}
