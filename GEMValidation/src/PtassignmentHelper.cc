#include "TF1.h"
#include "TGraphErrors.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "GEMCode/GEMValidation/interface/PtassignmentHelper.h"
//#include "GEMCode/GEMValidation/interface/DisplacedMuonTriggerPtassignment.h"
#include <iostream>
#include <math.h>       /* atan */

int PtassignmentHelper::GetEtaPartition(float eta)
{
  int neta=-1;
  if (fabs(eta)>=1.2 and fabs(eta)<1.4)
    neta=0;
  else if (fabs(eta)>=1.4 and fabs(eta)<1.6)
    neta=1;
  else if (fabs(eta)>=1.6 and fabs(eta)<1.8)
    neta=2;
  else if (fabs(eta)>=1.8 and fabs(eta)<2.0)
    neta=3;
  else if (fabs(eta)>=2.0 and fabs(eta)<2.2)
    neta=4;
  else if (fabs(eta)>=2.2 and fabs(eta)<2.4)
    neta=5;
  
  return neta;
}

float PtassignmentHelper::Ptassign_Position(float deltay12, float deltay23, float eta, int par){
    int neta = GetEtaPartition(eta);
    if (par<0 or par>3 or neta==-1) return -1;
    
    //std::cout <<"Pt position, npar "<< par <<" neta "<< neta <<" prop "<< PositionEpLUT[par][neta][0] <<" slope "<< PositionEpLUT[par][neta][1]<<" intercep "<< PositionEpLUT[par][neta][2] << " deltay12 " << deltay12 <<" deltay23 "<< deltay23 << std::endl;
    //if (fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))<0.005) return 100;
    //float pt=(1/fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))+PositionEpLUT[par][neta][2])/PositionEpLUT[par][neta][1]; 
    float pt=(1/fabs(deltay23-PositionEpLUT[par][neta][0]*deltay12)+PositionEpLUT[par][neta][2])/PositionEpLUT[par][neta][1]; 
    //std::cout <<" Ptassgin Position pt "<< pt << std::endl;
    return pt;
}

float PtassignmentHelper::deltaYcalculation(GlobalPoint gp1, GlobalPoint gp2){
   float anglea = gp2.phi();
   float newyst1 = -gp1.x()*sin(anglea) + gp1.y()*cos(anglea);
   float newyst2 = -gp2.x()*sin(anglea) + gp2.y()*cos(anglea);
   return (newyst2-newyst1);

}

float PtassignmentHelper::Ptassign_Position_gp(GlobalPoint gp1, GlobalPoint gp2, GlobalPoint gp3, float eta, int par){

   float anglea = gp2.phi();
   float newyst1 = -gp1.x()*sin(anglea) + gp1.y()*cos(anglea);
   float newyst2 = -gp2.x()*sin(anglea) + gp2.y()*cos(anglea);
	//float newxst3 = gp3.x()*cos(anglea) + gp3.y()*sin(anglea);
   float newyst3 = -gp3.x()*sin(anglea) + gp3.y()*cos(anglea);
   float deltay12 = newyst2-newyst1;
   float deltay23 = newyst3-newyst2;
   //std::cout <<" angle in st2 "<< anglea <<" newyst1 "<< newyst1 <<" newyst2 "<< newyst2 << " newyst3 "<< newyst3 << std::endl;
   
   return Ptassign_Position(deltay12,deltay23, eta, par);
}


float PtassignmentHelper::Ptassign_Direction(float bending_12, float eta, int par){
  int neta = GetEtaPartition(eta);
    if (par<0 or par>3 or neta==-1 or fabs(bending_12) > M_PI) return -1;

    float pt=(1/fabs(bending_12)+DirectionEpLUT[par][neta][1])/DirectionEpLUT[par][neta][0];
    //std::cout <<"Pt Direction, npar "<< par <<" neta "<<neta <<" slope "<< DirectionEpLUT[par][neta][0]<<" intercep "<< DirectionEpLUT[par][neta][1] << " bending_12 " << bending_12 <<" pt "<<pt <<std::endl;
    
    return pt;
}



float PtassignmentHelper::PhiMomentum(float dphi, float phi_position, int st, bool evenodd){
	
    //even->0, odd->1
    int cham = (evenodd ? 0:1);
    float slope = BendingAngleLUT[st-1][cham];
    float phi_m = dphi*slope+phi_position;
    //std::cout <<"st "<< st <<" cham "<< cham <<" gemcsc dphi "<< dphi <<" phi position "<< phi_position <<" slope "<< slope <<" phi momentum "<< phi_m << std::endl;
    return phi_m;


}


float PtassignmentHelper::PhiMomentum_Radius(float dphi, float phi_position, float radius_csc, float radius_gem){
    
     // usually radius_csc>radius_gem
     if (fabs(dphi) > M_PI or fabs(phi_position) > M_PI or radius_csc<radius_gem) return -9;
     float radius_diff = radius_gem-radius_csc*cos(dphi);
     float phi_diff = 0.0;
     if (fabs(radius_diff) > 0.0) phi_diff = atan(radius_csc*sin(dphi)/radius_diff);
     else phi_diff = M_PI/2.0;
     
     if (phi_diff <= -M_PI) phi_diff = phi_diff+2*M_PI;
     else if (phi_diff > M_PI) phi_diff = phi_diff-2*M_PI;

     float phiM = phi_position-phi_diff;
     if (phiM <= -M_PI) phiM = phiM+2*M_PI;
     else if (phiM > M_PI) phiM = phiM-2*M_PI;

     //std::cout <<" radius_csc "<< radius_csc <<" radius_gem "<< radius_gem <<" dphi "<< dphi << " phi_position "<< phi_position<<" radius_diff "<< radius_diff <<" phi_diff "<< phi_diff <<" phiM "<< phiM << std::endl;
     return phiM;

}


//dphi: local bending, phi_position: phi of GEM position, X: "x_factor"
//X: st1, X=D(GEM,CSC)*x_factor, st2: D(GEM,CSC)*x_factor/(D(ME11,ME21)*x+1)
float PtassignmentHelper::PhiMomentum_Xfactor(float dphi, float phi_position, float X){

      
     if (fabs(dphi) > M_PI or fabs(phi_position) > M_PI) return -9;
     float y = 1.0-cos(dphi)- X;
      
     float phi_diff = 0.0;
     if (fabs(y) > 0.0) phi_diff = atan(sin(dphi)/y);
     else phi_diff = M_PI/2.0;

     if (phi_diff <= -M_PI) phi_diff = phi_diff+2*M_PI;
     else if (phi_diff > M_PI) phi_diff = phi_diff-2*M_PI;

     float phiM = phi_position-phi_diff;
     if (phiM <= -M_PI) phiM = phiM+2*M_PI;
     else if (phiM > M_PI) phiM = phiM-2*M_PI;
     
     //std::cout <<"PhiMomentum_Xfactor: dphi "<< dphi <<" phi_position "<< phi_position <<" X "<<X <<" phi_diff "<< phi_diff <<" phiM "<< phiM << std::endl;
     return phiM;

}



float PtassignmentHelper::PhiMomentum_Xfactor_V2(float phi_CSC, float phi_GEM, float X){

     if (fabs(phi_CSC) > M_PI or fabs(phi_GEM) > M_PI) return -9;
     float dphi = deltaPhi(phi_CSC,phi_GEM);
     float y = 1.0-cos(dphi)- X;
      
     float phi_diff = 0.0;
     if (fabs(y) > 0.0) phi_diff = atan(sin(dphi)/y);
     else phi_diff = M_PI/2.0;

     if (phi_diff <= -M_PI) phi_diff = phi_diff+2*M_PI;
     else if (phi_diff > M_PI) phi_diff = phi_diff-2*M_PI;

     float phiM = phi_GEM-phi_diff;
     if (phiM <= -M_PI) phiM = phiM+2*M_PI;
     else if (phiM > M_PI) phiM = phiM-2*M_PI;
     
     //std::cout <<"PhiMomentum_Xfactor: dphi "<< dphi <<" phi_poistion1 "<< phi_GEM <<" phi_position2 "<< phi_CSC <<" Xfactor "<<X <<" phi_diff "<< phi_diff <<" phiM "<< phiM << std::endl;

     return phiM;

}



void PtassignmentHelper::calculateAlphaBeta(const std::vector<float>& v, 
                        const std::vector<float>& w, 
                        const std::vector<float>& ev, 
                        const std::vector<float>& ew, 
                        const std::vector<float>& status,
                        float& alpha, float& beta)
{
  //std::cout << "size of v: "<<v.size() << std::endl; 
  
  if (v.size()>=3) {
  
  float zmin;
  float zmax;
  if (v.front() < v.back()){
    zmin = v.front();
    zmax = v.back();
  }
  else{
    zmin = v.back();
    zmax = v.front();
  }

  TF1 *fit1 = new TF1("fit1","pol1",zmin,zmax); 
  //where 0 = x-axis_lowest and 48 = x_axis_highest 
  TGraphErrors* gr = new TGraphErrors(v.size(),&(v[0]),&(w[0]),&(ev[0]),&(ew[0]));
  gr->SetMinimum(w[2]-5*0.002);
  gr->SetMaximum(w[2]+5*0.002);
 
  gr->Fit(fit1,"RQ"); 
  
  alpha = fit1->GetParameter(0); //value of 0th parameter
  beta  = fit1->GetParameter(1); //value of 1st parameter

  delete fit1;
  delete gr;
  }
  else {alpha = -99; beta= 0.0;}
}



float PtassignmentHelper::normalizePhi(float phi) { 
    float result = phi;
    if(result > float(M_PI)) result -= float(2*M_PI);
    else if (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

