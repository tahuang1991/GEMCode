from ROOT import TCut

#_______________________________________________________________________________
def ANDtwo(cut1,cut2):
    """AND of two TCuts in PyROOT"""
    return TCut("(%s) && (%s)"%(cut1.GetTitle(),cut2.GetTitle()))


#_______________________________________________________________________________
def ORtwo(cut1,cut2):
    """OR of two TCuts in PyROOT"""
    return TCut("(%s) || (%s)"%(cut1.GetTitle(),cut2.GetTitle()))


#_______________________________________________________________________________
def AND(*arg):
    """AND of any number of TCuts in PyROOT"""
    length = len(arg)
    if length == 0:
        print "ERROR: invalid number of arguments"
        return
    if length == 1:
        return arg[0] 
    if length==2:
        return ANDtwo(arg[0],arg[1])
    if length>2:
        result = arg[0]
        for i in range(1,len(arg)):
            result = ANDtwo(result,arg[i])
        return result


#_______________________________________________________________________________
def OR(*arg): 
    """OR of any number of TCuts in PyROOT"""
    length = len(arg)
    if length == 0:
        print "ERROR: invalid number of arguments"
        return
    if length == 1:
        return arg[0] 
    if length==2:
        return ORtwo(arg[0],arg[1])
    if length>2:
        result = arg[0]
        for i in range(1,len(arg)):
            result = ORtwo(result,arg[i])
        return result


#_______________________________________________________________________________
## ranges of eta partitions GE1/1 and GE2/1L

eta_partitions = {
    'GE11-6' : {
    },
    'GE11-8-8' : {
    },
    'GE11-9-10' : {
        1: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.5416, 'max' : 1.59633} },
        2: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.59652, 'max' : 1.65744} },
        3: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.65765, 'max' : 1.7236} },
        4: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.72382, 'max' : 1.78146} },
        5: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.7817, 'max' : 1.84353} },
        6: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.84379, 'max' : 1.89934} },
        7: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.89961, 'max' : 1.95892} },
        8: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 1.95922, 'max' : 2.02108} },
        9: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : 2.0214, 'max' : 2.08765} },
        10: {'even' : {'min' : -99, 'max' : -99}, 'odd' : {'min' : 2.08799, 'max' : 2.15167} },
    },
    'GE21' : {
        1: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : , 'max' : } },
        2: {'even' : {'min' : , 'max' : }, 'odd' : {'min' : , 'max' : } },
    },
}
    

#_______________________________________________________________________________
nocut = TCut("")

ok_eta = TCut("TMath::Abs(eta)>1.64 && TMath::Abs(eta)<2.14")
ok_pt = TCut("pt > 20.")

## CSC simhits & digis
ok_sh1 = TCut("(has_csc_sh&1) > 0")
ok_sh2 = TCut("(has_csc_sh&2) > 0")
ok_st1 = TCut("(has_csc_strips&1) > 0")
ok_st2 = TCut("(has_csc_strips&2) > 0")
ok_w1 = TCut("(has_csc_wires&1) > 0")
ok_w2 = TCut("(has_csc_wires&2) > 0")
ok_digi1 = AND(ok_st1,ok_w1)
ok_digi2 = AND(ok_st2,ok_w2)

## CSC stub
ok_lct1 = TCut("(has_lct&1) > 0")
ok_lct2 = TCut("(has_lct&2) > 0")
ok_alct1 = TCut("(has_alct&1) > 0")
ok_alct2 = TCut("(has_alct&2) > 0")
ok_clct1 = TCut("(has_clct&1) > 0")
ok_clct2 = TCut("(has_clct&2) > 0")
ok_lct_hs_min = TCut("hs_lct_odd > 4")
ok_lct_hs_max = TCut("hs_lct_odd < 125")
ok_lct_hs = AND(ok_lct_hs_min,ok_lct_hs_max)
ok_lcths1 = AND(ok_lct1,ok_lct_hs)
ok_lcths2 = AND(ok_lct2,ok_lct_hs)

## GEM simhit
ok_gsh1 = TCut("(has_gem_sh&1) > 0")
ok_gsh2 = TCut("(has_gem_sh&2) > 0")
ok_g2sh1 = TCut("(has_gem_sh2&1) > 0")
ok_g2sh2 = TCut("(has_gem_sh2&2) > 0")


## GEM digi
ok_gdg1 = TCut("(has_gem_dg&1) > 0")
ok_gdg2 = TCut("(has_gem_dg&2) > 0")
ok_pad1 = TCut("(has_gem_pad&1) > 0")
ok_pad2 = TCut("(has_gem_pad&2) > 0")
ok_rpcstrip1 = TCut("(has_rpc_dg&1) > 0")
ok_rpcstrip2 = TCut("(has_rpc_dg&2) > 0")

ok_dphi1 = TCut("dphi_pad_odd < 10.")
ok_dphi2 = TCut("dphi_pad_even < 10.")

ok_pad1_lct1 = AND(ok_pad1,ok_lct1)
ok_pad2_lct2 = AND(ok_pad2,ok_lct2)

ok_rpcstrip1_lct1 = AND(ok_rpcstrip1,ok_lct1)
ok_rpcstrip2_lct2 = AND(ok_rpcstrip2,ok_lct2)

ok_pad1_dphi1 = AND(ok_pad1,ok_dphi1)
ok_pad2_dphi2 = AND(ok_pad2,ok_dphi2)

ok_lct1_eta = AND(ok_eta,ok_lct1)
ok_lct2_eta = AND(ok_eta,ok_lct2)

ok_pad1_lct1_eta = AND(ok_pad1,ok_lct1,ok_eta)
ok_pad2_lct2_eta = AND(ok_pad2,ok_lct2,ok_eta)

ok_gsh1_lct1_eta = AND(ok_gsh1,ok_lct1,ok_eta)
ok_gsh2_lct2_eta = AND(ok_gsh2,ok_lct2,ok_eta)

ok_gsh1_eta = AND(ok_gsh1,ok_eta)
ok_gsh2_eta = AND(ok_gsh2,ok_eta)

ok_gdg1_eta = AND(ok_gdg1,ok_eta)
ok_gdg2_eta = AND(ok_gdg2,ok_eta)

ok_2pad1 = TCut("(has_gem_pad2&1) > 0")
ok_2pad2 = TCut("(has_gem_pad2&2) > 0")

ok_pad1_overlap = OR(ok_pad1,AND(ok_lct2,ok_pad2))
ok_pad2_overlap = OR(ok_pad2,AND(ok_lct1,ok_pad1))

ok_copad1 = TCut("(has_gem_copad&1) > 0")
ok_copad2 = TCut("(has_gem_copad&2) > 0")

ok_Qp = TCut("charge > 0")
ok_Qn = TCut("charge < 0")

ok_lct1_eta_Qn = AND(ok_lct1,ok_eta,ok_Qn)
ok_lct2_eta_Qn = AND(ok_lct2,ok_eta,ok_Qn)

ok_lct1_eta_Qp = AND(ok_lct1,ok_eta,ok_Qp)
ok_lct2_eta_Qp = AND(ok_lct2,ok_eta,ok_Qp)

Ep = TCut("endcap > 0")
En = TCut("endcap < 0")

