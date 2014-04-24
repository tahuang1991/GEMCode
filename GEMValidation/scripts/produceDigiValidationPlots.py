import sys

from ROOT import *

from cuts import *
from drawPlots import *

## run quiet mode
import sys
sys.argv.append( '-b' )

import ROOT 
ROOT.gROOT.SetBatch(1)

from BaseValidation import *
from GEMDigiValidation import *
from RPCDigiValidation import *
from CSCDigiValidation import *

if __name__ == "__main__":  

  ## Style
  gStyle.SetStatStyle(0);
  plotter = DigiPlotter()

  rpcDigiOccupancyXY(plotter)
  """
  gemDigiOccupancyXY(plotter)
  gemDigiOccupancyStripPhi(plotter)
  gemDigiOccupancyStrip(plotter)
  gemDigiBX(plotter)
  gemDigiOccupancyRZ(plotter)
  
  gemPadOccupancyXY(plotter)
  gemPadOccupancyPadPhi(plotter)
  gemPadOccupancyPad(plotter)
  gemPadBX(plotter)
  gemPadOccupancyRZ(plotter)

  gemCoPadOccupancyXY(plotter)
  gemCoPadOccupancyCoPadPhi(plotter)
  gemCoPadOccupancyCoPad(plotter)
  gemCoPadBX(plotter)
  gemCoPadOccupancyRZ(plotter)
  
  simTrackDigiMatchingEta(plotter)
  simTrackDigiMatchingPhi(plotter)
  simTrackDigiMatchingLX(plotter)
  simTrackDigiMatchingLY(plotter)
  simTrackPadMatchingEta(plotter)
  simTrackPadMatchingPhi(plotter)
  simTrackPadMatchingLX(plotter)
  simTrackPadMatchingLY(plotter)
  simTrackCoPadMatchingEta(plotter)
  simTrackCoPadMatchingPhi(plotter)
  simTrackCoPadMatchingLX(plotter)
  simTrackCoPadMatchingLY(plotter)
  """
