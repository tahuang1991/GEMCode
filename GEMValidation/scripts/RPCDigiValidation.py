from ROOT import *

from cuts import *
from drawPlots import *

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT 
ROOT.gROOT.SetBatch(1)


#_______________________________________________________________________________
def rpcDigiOccupancyXY(plotter):
    
    ##loop on station
    for i in range(1,5):
        ## loop on region
        for j in range(1,3):
            ## loop on even/odd
            for k in range(1,4):
                
                postfix = ['_even', '_odd', '']
                region = ['rp1', 'rm1']
                regions = ['region1','region-1']
                title = "Digi occupancy: %s, station%d; globalX [cm]; globalY [cm]"%(regions[j-1], i)
                
                draw_occ(plotter.targetDir, "rpc_strip_dg_xy_%s_st%d%s"%(region[j-1], i, postfix[k-1]), plotter.ext, plotter.treeRPCDigis, title,
                         "h_", "(700,-700,700,700,-700,700)", "g_x:g_y", AND(re(j), st(i), evenodd[k-1]), "COLZ")
