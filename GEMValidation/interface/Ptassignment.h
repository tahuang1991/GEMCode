

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"




//eta partitions: 1.6-1.8, 1.8-2.0, 2.0-2.2, 2.2-2.4
enum {EtaPartitions=4, Parity=4};
const double PositionEpLUT[Parity][EtaPartitions][3] = {
    		//prop_factor, slope, intercept
		{{0.649, 0.05517, 0.08284},
		 {0.649, 0.08129, 0.1122 },
		 {0.649, 0.1682, 0.2233},
		 {0.649, 0.5304, 1.061},
		},
		{{0.3533, 0.1121, 0.2312},
		 {0.3533, 0.1593, 0.2771},
		 {0.3533, 0.3293, 0.5923},
		 {0.3533, 0.8649, 1.492},
		},
		{{0.5724, 0.04756, 0.06255},
		 {0.5724, 0.08478, 0.1368},
		 {0.5724, 0.1608, 0.1612},
		 {0.5724, 0.4944, 1.043},
		},
		{{0.3172, 0.1026, 0.1795},
		 {0.3172, 0.1495, 0.2236},
		 {0.3172, 0.3166, 0.5733},
		 {0.3172, 0.8348, 1.468},
		},
	};

float Ptassign_Position(float deltay12, float deltay23, float eta, int par){
    int neta=-1;
    if (fabs(eta)>=1.6 and fabs(eta)<1.8)
	neta=0;
    else if (fabs(eta)>=1.8 and fabs(eta)<2.0)
	neta=1;
    else if (fabs(eta)>=2.0 and fabs(eta)<2.2)
	neta=2;
    else if (fabs(eta)>2.2 and fabs(eta)<2.4)
	neta=3;
    else return -1;
    if (par<0 or par>3) return -1;
    
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


