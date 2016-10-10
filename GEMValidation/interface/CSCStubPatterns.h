#ifndef MuJetAnalysis_DisplacedL1MuFilter_CSCStubPatterns_h
#define MuJetAnalysis_DisplacedL1MuFilter_CSCStubPatterns_h

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include <vector>

 std::vector<std::vector<int> > pat0delta {  
   { 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999}, 
     { 999, 999, 999, 999, 999}, 
       {999},            // pid=999: no pattern found 
         {999, 999, 999, 999, 999}, 
           {999, 999, 999, 999, 999, 999, 999, 999, 999}, 
             {999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999} 
 }; 

 std::vector<std::vector<int> > pat1delta {  
   { -5,   -4,   -3,   -2,   -1,   0,   1,   2,   3,   4,   5}, 
     {-2,   -1,   0,   1,   2}, 
       {0},             // pid=1: layer-OR trigger 
         {-2,   -1,   0,   1,   2}, 
           { -4,   -3,   -2,   -1,   0,   1,   2,   3,   4}, 
             { -5,   -4,   -3,   -2,   -1,   0,   1,   2,   3,   4,   5}, 
 }; 

 std::vector<std::vector<int> > pat2delta {  
   { 999, 999, 999, 999, 999, 999, 999, 999,   3,   4,   5}, 
     {999, 999, 999,   1,   2}, 
       {0},             // pid=2: right-bending (large) 
         {-2,   -1,   0, 999, 999}, 
           {-4,   -3,   -2, 999, 999, 999, 999, 999, 999}, 
             {-5,   -4,   -3, 999, 999, 999, 999, 999, 999, 999, 999} 
 }; 

 std::vector<std::vector<int> > pat3delta { 
   {   -5,   -4,   -3, 999, 999, 999, 999, 999, 999, 999, 999}, 
     {-2,   -1, 999, 999, 999}, 
       {0},             // pid=3: left-bending (large) 
         {999, 999,   0,   1,   2}, 
           {999, 999, 999, 999, 999, 999,   2,   3,   4}, 
             {999, 999, 999, 999, 999, 999, 999, 999,   3,   4,   5} 
 };  

 std::vector<std::vector<int> > pat4delta { 
   { 999, 999, 999, 999, 999, 999, 999,   2,   3,   4, 999}, 
     {999, 999, 999,   1,   2}, 
       {0},             // pid=4: right-bending (medium) 
         {-2,   -1, 999, 999, 999}, 
           {-4,   -3,   -2, 999, 999, 999, 999, 999, 999}, 
             {999,   -4,   -3,   -2, 999, 999, 999, 999, 999, 999, 999} 
            
 }; 

 std::vector<std::vector<int> > pat5delta { 
   { 999,   -4,   -3,   -2, 999, 999, 999, 999, 999, 999, 999}, 
     {-2,   -1, 999, 999, 999}, 
       {0},             // pid=5: left-bending (medium) 
         {999, 999, 999,   1,   2}, 
           {999, 999, 999, 999, 999, 999,   2,   3,   4}, 
             {999, 999, 999, 999, 999, 999, 999,   2,   3,   4, 999} 
            
 }; 

 std::vector<std::vector<int> > pat6delta { 
   { 999, 999, 999, 999, 999, 999,   1,   2,   3, 999, 999}, 
     {999, 999,   0,   1, 999}, 
       {0},             // pid=6: right-bending (medium) 
         {999,   -1,   0, 999, 999}, 
           {999, 999,   -2,   -1, 999, 999, 999, 999, 999}, 
             {999, 999,   -3,   -2,   -1, 999, 999, 999, 999, 999, 999} 
 }; 

 std::vector<std::vector<int> > pat7delta { 
   { 999, 999,   -3,   -2,   -1, 999, 999, 999, 999, 999, 999}, 
     {999,   -1,   0, 999, 999}, 
       {0},             // pid=7: left-bending (medium) 
         {999, 999,   0,   1, 999}, 
           {999, 999, 999, 999, 999,   1,   2, 999, 999}, 
             {999, 999, 999, 999, 999, 999,  1,  2,   3, 999, 999} 
 }; 

 std::vector<std::vector<int> > pat8delta { 
   { 999, 999, 999, 999, 999,   0,   1,   2, 999, 999, 999}, 
     {999, 999,   0,   1, 999}, 
       {0},             // pid=8: right-bending (small) 
         {999,   -1,   0, 999, 999}, 
           {999, 999,   -2,   -1,   0, 999, 999, 999, 999}, 
             {999, 999, 999,   -2,   -1,   0, 999, 999, 999, 999, 999} 
 }; 

 std::vector<std::vector<int> > pat9delta { 
   { 999, 999, 999,   -2,   -1,   0, 999, 999, 999, 999, 999}, 
     {999,   -1,   0, 999, 999}, 
       {0},             // pid=9: left-bending (small) 
         {999, 999,   0,   1, 999}, 
           {999, 999, 999, 999,   0,   1,   2, 999, 999}, 
             {999, 999, 999, 999, 999,   0,   1,   2, 999, 999, 999} 
 }; 

 std::vector<std::vector<int> > patAdelta { 
   { 999, 999, 999, 999,   -1,   0,   1, 999, 999, 999, 999}, 
     {999, 999,   0, 999, 999}, 
       {0},             // pid=A: straight-through 
         {999, 999,   0, 999, 999}, 
           {999, 999, 999,   -1,   0,   1, 999, 999, 999}, 
             {999, 999, 999, 999,   -1,   0,   1, 999, 999, 999, 999} 
 }; 

std::vector< std::vector<std::vector<int> > > patIndexToPatternDelta { 
  pat0delta, pat1delta, pat2delta, pat3delta, pat4delta, pat5delta, pat6delta, pat7delta, pat8delta, pat9delta, patAdelta 
 }; 

bool 
comparatorInLCTPattern(int keyStrip, int pattern, int layer, int halfStrip)
{
  bool debug = false;
  // first, get the (sub)pattern
  bool returnValue = false;
  std::vector<int> subpat = patIndexToPatternDelta[pattern].at(layer-1);
  if (debug) for (auto p: subpat) std::cout << "\t" << p << std::endl;
  if (debug) std::cout << "\tkeyStrip pattern layer halfstrip " << keyStrip << " " <<pattern << " " <<layer << " " <<halfStrip <<std::endl <<std::endl;
  // due to comparator digi time extension in the CLCT processor we need to  
  // search a bigger region around the key HS. +/-1, 0 should be sufficient
  int halfStripDelta = halfStrip - keyStrip;
  returnValue = ( std::find(subpat.begin(), subpat.end(), halfStripDelta+1) != subpat.end() or
                  std::find(subpat.begin(), subpat.end(), halfStripDelta)   != subpat.end() or
                  std::find(subpat.begin(), subpat.end(), halfStripDelta-1) != subpat.end() );
  return returnValue;
}

#endif
