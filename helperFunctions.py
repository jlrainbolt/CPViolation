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
# Takes TLorentzVectors (list) and charges (dict)
# Pairing is determined by maximizing mass difference


def GetMassPairs(p4_list, q_dict):

    ##  INITIALIZE  ##

    # Assign highest-p lepton as lep1, which is always in pair 1
    p4_list.sort(key=lambda fourvec: fourvec.P(), reverse=True)

    lep1p4, lep1q = p4_list[0], q_dict[p4_list[0]]
    p1_list, p2_list = [lep1p4], []

    config = 1



    ##  CONFIG 1  ##

    # OS leptons with largest momenta are Z1

    # Find subleading lepton of pair 1
    for lep in p4_list:
        if q_dict[lep] != lep1q:
            p1_list.append(lep)
            break


    # Assemble pair 2
    for lep in p4_list:
        if lep not in p1_list:
            p2_list.append(lep)


    # Find difference between masses of pairs and assign them
    pair1, pair2 = GetSum(p1_list), GetSum(p2_list)
    mass_diff_1 = pair1.M() - pair2.M()

    if mass_diff_1 > 0:
        p_pair, k_pair = pair1, pair2
        p_list, k_list = p1_list[:], p2_list[:]
    else:
        k_pair, p_pair = pair1, pair2
        k_list, p_list = p1_list[:], p2_list[:]



    ##  CONFIG 2  ##

    # Swap leptons that are not SS as lep1
    # (i.e. look for pair 2 lepton with SS as subleading lep of pair 1)
    for i in range(len(p2_list)):
        if q_dict[p2_list[i]] == q_dict[p1_list[1]]:
            p2_list[i], p1_list[1] = p1_list[1], p2_list[i]
            break


    # Reassign pair momenta
    pair1, pair2 = GetSum(p1_list), GetSum(p2_list)
    mass_diff_2 = pair1.M() - pair2.M()


    # If this config has larger mass difference, reassign returned values
    if abs(mass_diff_2) > abs(mass_diff_1):
        config = 2

        if mass_diff_2 > 0:
            p_pair, k_pair = pair1, pair2
            p_list, k_list = p1_list[:], p2_list[:]
        else:
            k_pair, p_pair = pair1, pair2
            k_list, p_list = p1_list[:], p2_list[:]



    ##  CHARGES  ##

    if q_dict[p_list[0]] > 0:
        p_plus, p_minus = p_list[0], p_list[1]
    else:
        p_minus, p_plus = p_list[0], p_list[1]

    if q_dict[k_list[0]] > 0:
        k_plus, k_minus = k_list[0], k_list[1]
    else:
        k_minus, k_plus = k_list[0], k_list[1]



    ##  DEBUG  ##

    if (False):
        print("Configuration", config)
        print(mass_diff_1, "(1)", mass_diff_2, "(2)")
        print("\npp mass:", p_pair.M(), (p_plus + p_minus).M())
        print(p_list[0].P(), q_dict[p_list[0]])
        print(p_list[1].P(), q_dict[p_list[1]])
        print("\nkk mass:", k_pair.M(), (k_plus + k_minus).M())
        print(k_list[0].P(), q_dict[k_list[0]])
        print(k_list[1].P(), q_dict[k_list[1]])
        print("\n")
        sys.stdout.flush()



    ##  RETURN ##
    return p_plus, p_minus, k_plus, k_minus
