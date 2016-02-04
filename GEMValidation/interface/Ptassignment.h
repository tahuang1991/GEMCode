

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"




//eta partitions:1.2-1.4,1.4-1.6,1.6-1.8, 1.8-2.0, 2.0-2.2, 2.2-2.4
enum {EtaPartitions=6, Parity=4};
const double PositionEpLUT[Parity][EtaPartitions][3] = {
    		//prop_factor, slope, intercept
		{{1.279, 0.04784, 0.1122}, //eta 1.2-1.4
		 {1.279, 0.65424, 0.09761},
    		 {0.648, 0.05527, 0.08944},
		 {0.648, 0.08295, 0.1294 },
		 {0.648, 0.1660, 0.2158},
		 {0.648, 0.4952, 0.7103},
		},
		{{0.6357, 0.0827, 0.2021},
		 {0.6357, 0.0906, 0.1773},
		 {0.3542, 0.1067, 0.1957},
		 {0.3542, 0.1561, 0.2645},
		 {0.3542, 0.3156, 0.4514},
		 {0.3542, 0.8242, 1.0712},
		},
		{{1.001, 0.038, 0.008345},
		 {1.001, 0.04157, 0.0617},
		 {0.5636, 0.0562, 0.08417},
		 {0.5636, 0.0870, 0.1426},
		 {0.5636, 0.1676, 0.2198},
		 {0.5636, 0.4953, 0.7272},
		},
		{{0.5252, 0.0739, 0.1714},
		 {0.5252, 0.07838, 0.1307},
		 {0.3217, 0.1026, 0.2026},
		 {0.3217, 0.1435, 0.2118},
		 {0.3217, 0.2874, 0.4055},
		 {0.3217, 0.7625, 1.468},
		},
	};


const double DirectionEpLUT[Parity][EtaPartitions][2]={
    	       {{2.907, 5.906},
		{2.600, 5.191},
		{4.530, 9.442},
		{5.788, 9.743},
		{8.367, 10.22},
		{11.02, 14.84},
	       },
	       {{2.409, 5.198},
		{2.467, 4.397},
		{4.779, 9.954},
		{6.273, 11.91},
		{9.315, 12.21},
		{10.34, 11.02},
	       },
	       {{2.301, 4.929},
		{2.230, 3.111},
		{7.677, 16.82},
		{7.726, 13.36},
		{9.621, 10.62},
		{11.23, 13.44},
	       },
	       {{2.401, 4.758},
		{2.383, 3.782},
		{7.720, 16.91},
		{8.643, 16.45},
		{10.02, 11.83},
		{11.83, 17.66},
	       },
	};

int GetEtaPartition_direction(float eta ){

    int neta=-1;
    if (fabs(eta)>=1.6 and fabs(eta)<1.8)
	neta=0;
    else if (fabs(eta)>=1.8 and fabs(eta)<2.0)
	neta=1;
    else if (fabs(eta)>=2.0 and fabs(eta)<2.2)
	neta=2;
    else if (fabs(eta)>2.2 and fabs(eta)<2.4)
	neta=3;

    return neta;

}

int GetEtaPartition_position(float eta ){

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
    else if (fabs(eta)>2.2 and fabs(eta)<2.4)
	neta=5;

    return neta;

}



float Ptassign_Position(float deltay12, float deltay23, float eta, int par){
    int neta = GetEtaPartition_position(eta);
    if (par<0 or par>3 or neta==-1) return -1;
    
    //std::cout <<" npar "<< par <<" neta "<< neta <<" prop "<< PositionEpLUT[par][neta][0] <<" slope "<< PositionEpLUT[par][neta][1]<<" intercep "<< PositionEpLUT[par][neta][2] << " ddY " <<fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12) << std::endl;
    if (fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))<0.005) return 100;
    float pt=(1/fabs(fabs(deltay23)-PositionEpLUT[par][neta][0]*fabs(deltay12))+PositionEpLUT[par][neta][2])/PositionEpLUT[par][neta][1]; 
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
    int neta = GetEtaPartition_direction(eta);
    if (par<0 or par>3 or neta==-1) return -1;

    //std::cout <<"Direction Based method  npar "<< par <<" neta "<<neta <<" slope "<< DirectionEpLUT[par][neta][1]<<" intercep "<< DirectionEpLUT[par][neta][1] << " bending_12 " << bending_12 << std::endl;
    float pt=(1/bending_12+DirectionEpLUT[par][neta][1])/DirectionEpLUT[par][neta][0];
    
    return pt;
}
