

#include "MuJetAnalysis/HLTBendingAngle/interface/Ptassignment.h"
#include <iostream>


int GetEtaPartition(float eta ){

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

float Ptassign_Position(float deltay12, float deltay23, float eta, int par){
    int neta = GetEtaPartition(eta);
    if (par<0 or par>3 or neta==-1) return -1;
    
    //std::cout <<"Pt position, npar "<< par <<" neta "<< neta <<" prop "<< PositionEpLUT[par][neta][0] <<" slope "<< PositionEpLUT[par][neta][1]<<" intercep "<< PositionEpLUT[par][neta][2] << " deltay12 " << deltay12 <<" deltay23 "<< deltay23 << std::endl;
    //if (fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))<0.005) return 100;
    //float pt=(1/fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))+PositionEpLUT[par][neta][2])/PositionEpLUT[par][neta][1]; 
    float pt=(1/fabs(deltay23-PositionEpLUT[par][neta][0]*deltay12)+PositionEpLUT[par][neta][2])/PositionEpLUT[par][neta][1]; 
    //std::cout <<" Ptassgin Position pt "<< pt << std::endl;
    return pt;
}


float Ptassign_Position_gp(GlobalPoint gp1, GlobalPoint gp2, GlobalPoint gp3, float eta, int par){

   float anglea = gp2.phi();
   float newyst1 = -gp1.x()*sin(anglea) + gp1.y()*cos(anglea);
   float newyst2 = -gp2.x()*sin(anglea) + gp2.y()*cos(anglea);
	//float newxst3 = gp3.x()*cos(anglea) + gp3.y()*sin(anglea);
   float newyst3 = -gp3.x()*sin(anglea) + gp3.y()*cos(anglea);
   float deltay12 = newyst2-newyst1;
   float deltay23 = newyst3-newyst2;
   
   return Ptassign_Position(deltay12,deltay23, eta, par);
}


float Ptassign_Direction(float bending_12, float eta, int par){
    int neta = GetEtaPartition(eta);
    if (par<0 or par>3 or neta==-1) return -1;

    //std::cout <<"Pt Direction, npar "<< par <<" neta "<<neta <<" slope "<< DirectionEpLUT[par][neta][0]<<" intercep "<< DirectionEpLUT[par][neta][1] << " bending_12 " << bending_12 << std::endl;
    float pt=(1/fabs(bending_12)+DirectionEpLUT[par][neta][1])/DirectionEpLUT[par][neta][0];
    
    return pt;
}
