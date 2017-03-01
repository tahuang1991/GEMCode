#include "GEMCode/GEMValidation/interface/BarrelTriggerPtAssignmentHelper.h"

// barrel trigger 2-station LUTs
std::map<std::string,std::vector<float> > DirectionBasedPtLUT_2stubs = {
  {"DT1_DT2", {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} },
  {"DT1_DT3", {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} },
  {"DT1_DT4", {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} },
  {"DT2_DT3", {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} },
  {"DT2_DT4", {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} },
  {"DT3_DT4", {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0, 28.0, 32.0, 36.0, 42.0} }
};

std::map<std::string,std::vector<float> > DirectionBasedDeltaPhiLUT_2stubs = {
  {"DT1_DT2", {0.4409166666666667, 0.1908051948051948, 0.11985024154589373, 0.09283333333333332, 0.07680000000000001, 0.0643731343283582, 0.05664723926380369, 0.0498625, 0.04719396551724138, 0.043472727272727274, 0.04341355932203391, 0.04370038910505837, 0.044390070921985825, 0.045082089552238816, 0.04535570469798658, 0.04629464285714286, 0.04681428571428572, 0.04832500000000001} },
  {"DT1_DT3", {0.7153999999999999, 0.52405, 0.28248750000000006, 0.19926244343891406, 0.16401477832512315, 0.13742857142857143, 0.11878896103896104, 0.10268253968253971, 0.08850628366247756, 0.07799284009546541, 0.07045631067961167, 0.06336458333333332, 0.05774454148471616, 0.05411052631578949, 0.0491657458563536, 0.047915094339622655, 0.04808536585365854, 0.048833333333333354} },
  {"DT1_DT4", {0.7152, 0.45605, 0.30988349514563107, 0.24174789915966385, 0.19359615384615386, 0.15707772020725388, 0.13356491228070178, 0.1064389534883721, 0.08990909090909091, 0.08115068493150687, 0.0722357142857143, 0.0646298076923077, 0.057612021857923516, 0.05361538461538461, 0.04822857142857144, 0.0498918918918919, 0.047575757575757584} },
  {"DT2_DT3", {0.5790000000000001, 0.38479411764705884, 0.19879850746268657, 0.1410168539325843, 0.11312716763005781, 0.09687843137254902, 0.08549166666666667, 0.07561097256857856, 0.06696994535519125, 0.05906714628297363, 0.054882575757575776, 0.049907161803713525, 0.05138912133891214, 0.0538978494623656, 0.05092452830188681, 0.04774766355140187, 0.05181159420289857, 0.053333333333333344} },
  {"DT2_DT4", {0.6513333333333334, 0.38794736842105265, 0.26101298701298703, 0.1996214285714286, 0.15861006289308177, 0.13258620689655173, 0.11008032128514056, 0.090284046692607, 0.0772006920415225, 0.06872500000000001, 0.06323648648648648, 0.05838075313807533, 0.05675675675675676, 0.05637974683544305, 0.05084782608695654, 0.0514047619047619, 0.052440000000000014} },
  {"DT3_DT4", {0.46840000000000004, 0.25570270270270273, 0.16560550458715598, 0.12402958579881657, 0.09708730158730158, 0.08288000000000001, 0.06968926553672317, 0.058116853932584284, 0.052528957528957546, 0.04980503144654089, 0.04863445378151261, 0.04883177570093458, 0.049464135021097054, 0.04774496644295303, 0.04876000000000001, 0.04692857142857144, 0.04775757575757576} }
};

// barrel trigger 3/4-station LUTs
std::map<std::string,std::map<int, std::vector<float> > > DirectionBasedDeltaPhiLUT_3or4stubs = {
  {"DT1_DT2__DT1_DT3", {
      {3 , { 0.15 , 0.49 , 32.0 } }, //Acceptance 0.900069294528  //Rejection 0.681036911806
      {5 , { 0.13 , 0.45 , 30.0 } }, //Acceptance 0.900258851035  //Rejection 0.76246830093
      {7 , { 0.11 , 0.15 , 33.0 } }, //Acceptance 0.90093437602  //Rejection 0.842209072978
      {10 , { 0.09 , 0.11 , 20.0 } }, //Acceptance 0.903721366586  //Rejection 0.898562975486
      {15 , { 0.11 , 0.07 , 28.0 } }, //Acceptance 0.900597585017  //Rejection 0.937447168216
      {20 , { 0.05 , 0.17 , 27.0 } }, //Acceptance 0.900077259851  //Rejection 0.942519019442
      {30 , { 0.25 , 0.05 , 22.0 } }, //Acceptance 0.901083537225  //Rejection 0.957171034094
      {40 , { 0.23 , 0.05 , 25.0 } }, //Acceptance 0.901442307692  //Rejection 0.964778810933
    }
  }
  /* {"DT1_DT2__DT1_DT4", {}}, */
  /* {"DT1_DT2__DT2_DT3", {}}, */
  /* {"DT1_DT2__DT2_DT4", {}}, */
  /* {"DT1_DT3__DT1_DT4", {}}, */
  /* {"DT1_DT3__DT2_DT3", {}}, */
  /* {"DT1_DT3__DT3_DT4", {}}, */
  /* {"DT1_DT4__DT2_DT4", {}}, */
  /* {"DT1_DT4__DT3_DT4", {}}, */
  /* {"DT2_DT3__DT2_DT4", {}}, */
  /* {"DT2_DT3__DT3_DT4", {}}, */
  /* {"DT3_DT4__DT3_DT4", {}} */
  /* {"DT1_DT2__DT3_DT4", {}}, */
  /* {"DT1_DT3__DT2_DT4", {}}, */
  /* {"DT1_DT4__DT2_DT3", {}} */
};


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
  float default_pt = 2.;

  std::map<int, std::vector<float> > lut = DirectionBasedDeltaPhiLUT_3or4stubs[DT_type];
  for (auto const& p : lut){
    // get the parameters for this ellipse
    float a_over_2 = p.second[0];
    float b_over_2 = p.second[1];
    float alpha_deg = p.second[2]; //degrees

    if (passEllipse(DPhi1, DPhi2, a_over_2, b_over_2, alpha_deg)){
      return float(p.first);
    }
  }

  return default_pt;
}

std::string BarrelTriggerPtAssignmentHelper::getBestDPhiPair(int station1, int station2, int station3)
{
  if (station1==1 and station2==2 and station3==3) return "DT1_DT2__DT1_DT3";
  if (station1==1 and station2==2 and station3==4) return "DT1_DT2__DT1_DT4";
  if (station1==1 and station2==3 and station3==4) return "DT1_DT3__DT1_DT4";
  if (station1==2 and station2==3 and station3==4) return "DT2_DT3__DT2_DT4";
  return "";
}

std::string BarrelTriggerPtAssignmentHelper::getBestDPhiPair(int station1, int station2, int station3, int station4)
{
  return "DT1_DT4__DT2_DT3";
}

