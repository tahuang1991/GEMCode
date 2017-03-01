#include "GEMCode/GEMValidation/interface/BarrelTriggerPtAssignmentHelper.h"

float BarrelTriggerPtAssignmentHelper::getEllipse(float x, float y, float a, float b, float alpha, float x0, float y0)
{
  // convert the angle to radian
  float alpha_rad = alpha / 180 * M_PI;
  float x1 = x*cos(alpha_rad)+y*sin(alpha_rad)-x0;
  float y1 = x*sin(alpha_rad)-y*cos(alpha_rad)-y0;
  return x1*x1/(a*a)+y1*y1/(b*b);
}

bool BarrelTriggerPtAssignmentHelper::passEllipse(float x, float y, float a, float b, float alpha, float x0, float y0)
{
  return getEllipse(x,y,a,b,alpha, x0, y0) <= 1.0;
}

bool BarrelTriggerPtAssignmentHelper::failEllipse(float x, float y, float a, float b, float alpha, float x0, float y0)
{
  return getEllipse(x,y,a,b,alpha,x0, y0) > 1.0;
}

float BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(float DPhi, int DT_type)
{
  std::string dt_type = DT_Types_2stubs_string[DT_type];
  return getDirectionBasedPt2Stubs(DPhi, dt_type);
}

float BarrelTriggerPtAssignmentHelper::getDirectionBasedPt2Stubs(float DPhi, std::string DT_type)
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

float BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(float DPhi1, float DPhi2, int DT_type)
{
  const std::string dt_type = DT_Types_3or4stubs_string[DT_type];
  return getDirectionBasedPt3or4Stubs(DPhi1, DPhi2, dt_type);
}

float BarrelTriggerPtAssignmentHelper::getDirectionBasedPt3or4Stubs(float DPhi1, float DPhi2, std::string DT_type)
{
  float default_pt = 2;

  std::map<int, std::vector<float> > lut = DirectionBasedDeltaPhiLUT_3or4stubs[DT_type];
  for (auto const& p : lut){
    // get the parameters for this ellipse
    float a_over_2 = p.second[0];
    float b_over_2 = p.second[1];
    float alpha_deg = p.second[2]; //degrees

    if (passEllipse(DPhi1, DPhi2, a_over_2, b_over_2, alpha_deg)){
      return p.first;
    }
  }

  return default_pt;
}
