from __future__ import print_function
from __future__ import division

import sys
import ROOT as rt




###############
#   HELPERS   #
###############


##  GET_BOOSTED  ##

# Returns TLorentzVector "p4" boosted by TVector3 "beta"
def GetBoosted(p4_, beta):
    p4 = p4_.Clone()
    p4.Boost(-beta)
    return p4



##  GET_SUM  ##

# Returns (vector) sum of TLorentzVectors in list
def GetSum(p4_list):
    p4 = rt.TLorentzVector()
    for p4_ in p4_list:
        p4 += p4_
    return p4



##  SIGN  ##

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return x




####################
#  GET_MASS_PAIRS  #
####################

# Returns high- and low-mass pairs for same-flavor final state
# Takes TLorentzVectors (sorted by decreasing momentum) and charges (dict)
# Pairing is determined by maximizing mass difference


def GetMassPairs(p4_list, q_dict):

    ##  INITIALIZE  ##

    # Assign highest-p lepton as lep1, which is always in pair 1
    lep1p4, lep1q = p4_list[0], q_dict[p4_list[0]]
    pair1_list, pair2_list = [lep1p4], []

    config = 1



    ##  CONFIG 1  ##

    # OS leptons with largest momenta are Z1

    # Find subleading lepton of pair 1
    for lep in p4_list:
        if q_dict[lep] != lep1q:
            pair1_list.append(lep)
            break


    # Assemble pair 2
    for lep in p4_list:
        if lep not in pair1_list:
            pair2_list.append(lep)


    # Find difference between masses of pairs and assign them
    pair1p4, pair2p4 = GetSum(pair1_list), GetSum(pair2_list)
    mass_diff_1 = pair1p4.M() - pair2p4.M()

    if mass_diff_1 > 0:
        large_pair, small_pair = pair1p4, pair2p4
        large_list, small_list = pair1_list, pair2_list
    else:
        large_pair, small_pair = pair2p4, pair1p4
        large_list, small_list = pair2_list, pair1_list



    ##  CONFIG 2  ##

    # Swap leptons that are not SS as lep1
    # (i.e. look for pair 2 lepton with SS as subleading lep of pair 1)
    for i in range(len(pair2_list)):
        if q_dict[pair2_list[i]] == q_dict[pair1_list[1]]:
            pair2_list[i], pair1_list[1] = pair1_list[1], pair2_list[i]
            break


    # Reassign pair momenta
    pair1p4, pair2p4 = GetSum(pair1_list), GetSum(pair2_list)
    mass_diff_2 = pair1p4.M() - pair2p4.M()


    # If this config has larger mass difference, reassign returned values
    if abs(mass_diff_2) > abs(mass_diff_1):
        config = 2

        if mass_diff_2 > 0:
            large_pair, small_pair = pair1p4, pair2p4
            large_list, small_list = pair1_list, pair2_list
        else:
            large_pair, small_pair = pair2p4, pair1p4
            large_list, small_list = pair2_list, pair1_list



    ##  DEBUG  ##

#   if (False):
    if (small_pair.M() < 4):
        print("Configuration", config)
        print(mass_diff_1, "(1)", mass_diff_2, "(2)")
        print("\nPair 1: mass", large_pair.M())
        print(pair1_list[0], q_dict[pair1_list[0]])
        print(pair1_list[1], q_dict[pair1_list[1]])
        print("\nPair 2: mass", small_pair.M())
        print(pair2_list[0], q_dict[pair2_list[0]])
        print(pair2_list[1], q_dict[pair2_list[1]])
        print("\n")
        sys.stdout.flush()



    ##  RETURN ##
    return large_pair, small_pair, large_list, small_list




########################
#  GET_TRIPLE_PRODUCT  #
########################

# Returns the scalar triple product psi = p_c . (p_a * p_b)
#   and the quantity sin(phi) = (n_ab * n_cd) . z
#
# Takes TLorentzVectors (sorted by decreasing momentum) and charges (dict)
# Choose p_a and p_c to be antileptons with |p_a| > |p_c|
# Choose p_b as lepton with largest momentum


def GetTripleProduct(p4_list, q_dict):

    # Separate leptons by charge
    plus_leps, minus_leps = [], []

    for lep in p4_list:
        if q_dict[lep] > 0:
            plus_leps.append(lep)
        else:
            minus_leps.append(lep)


    # Debug
#   if (False):
    if (len(plus_leps) != 2) or (len(minus_leps) != 2):
        for lep in plus_leps:
            print(lep.E(), end=", ")
        print("\n", end="")
        for lep in minus_leps:
            print(lep.E(), end=", ")
#       print(plus_leps)
#       print(minus_leps)
        print("")
        sys.stdout.flush()
        return float('nan')


    # Assign as TVector3s
    lepAp3, lepCp3 = plus_leps[0].Vect(), plus_leps[1].Vect()
    lepBp3, lepDp3 = minus_leps[0].Vect(), minus_leps[1].Vect()


    # Calculate normals to planes (a, b) and (c, d)
    # and direction of p_a + p_b
    N_ab, N_cd = lepAp3.Cross(lepBp3), lepCp3.Cross(lepDp3)
    n_ab, n_cd = N_ab.Unit(), N_cd.Unit()
    z = (lepAp3 + lepBp3).Unit()


    # Calculate scalar triple product and sin(phi)
    psi = lepCp3.Dot(N_ab)
    phi = N_ab.Angle(N_cd)
    sin_phi = n_ab.Cross(n_cd).Dot(z)



    ##  RETURN  ##

    return psi, phi, sin_phi



##############################
#  GET_TRIPLE_PRODUCT_PAIRS  #
##############################

# Returns the scalar triple product psi = p_c . (p_a * p_b)
#   and the quantity sin(phi) = (n_ab * n_cd) . z
#
# Takes TLorentzVector pairs (two-item lists) and charges (dict)
# Choose p_a and p_b to be + and - q leptons of pair 1
# Choose p_c and p_d to be + and - q leptons of pair 2


def GetTripleProductPairs(p1_leps, p2_leps, q_dict):

    # Assign as TVector3s and check sign
    lepAp3, lepBp3 = p1_leps[0].Vect(), p1_leps[1].Vect()
    if q_dict[p1_leps[0]] < 0:
        lepAp3, lepBp3 = lepBp3, lepAp3

    lepCp3, lepDp3 = p2_leps[0].Vect(), p2_leps[1].Vect()
    if q_dict[p2_leps[0]] < 0:
        lepCp3, lepDp3 = lepDp3, lepCp3


    # Calculate normals to planes (a, b) and (c, d)
    # and direction of p_a + p_b
    N_ab, N_cd = lepAp3.Cross(lepBp3), lepCp3.Cross(lepDp3)
    n_ab, n_cd = N_ab.Unit(), N_cd.Unit()
    z = (lepAp3 + lepBp3).Unit()


    # Calculate scalar triple product and sin(phi)
    psi = lepCp3.Dot(N_ab)
    phi = N_ab.Angle(N_cd)
    sin_phi = n_ab.Cross(n_cd).Dot(z)



    ##  RETURN  ##

    return psi, phi, sin_phi
