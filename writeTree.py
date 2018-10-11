from __future__ import print_function
from __future__ import division

import sys
import ROOT as rt
from array import array
from LHEevent import *
from LHEfile import *
from helperFunctions import *
from scipy.io import loadmat
from scipy.interpolate import PPoly




if __name__ == '__main__':


    ################
    #  INITIALIZE  #
    ################

    debug = False

    ##  INPUT LHE FILE  ##

    lheName = "unweighted_events_fullSM.lhe"

    # Counters
    n_max = 100000
    if debug:
        n_max = 100
    n_evts = n_4m = n_2m2e = n_2e2m = n_4e = 0


    
    ##  FIT FUNCTION  ##

    # (for 1/alpha_QED; used to estimate theo. unc.)

    # Import "smoothingspline" fit struct from Matlab
    alpha_inv_p = loadmat("alpha_inv_p.mat")

    # Recreate as piecewise polynomial function
    # (takes mass; returns 1/alpha)
    get_alpha_inv = PPoly(alpha_inv_p['coefs'].transpose(), alpha_inv_p['breaks'][0])

    # Default MG5 value for 1/alpha_QED
    aEWM1 = 1.325070e+02



    ##  OUTPUT ROOT FILE  ##

    outName = "4l_fullSM.root"
    outFile = rt.TFile(outName, "RECREATE")
    tree = rt.TTree("tree", "")



    ##  TREE BRANCHES  ##

    # (pyROOT requires numbers to be stored as 1D arrays)

    # CPV observables
    triple_product, sin_phi = array("f", [0]), array("f", [0])

    tree.Branch("triple_product", triple_product, "triple_product/F")
    tree.Branch("sin_phi", sin_phi, "sin_phi/F")


    # Lepton pairs (sorted by mass)
    z1p4, z1pdg = rt.TLorentzVector(), array("i", [0])
    z2p4, z2pdg = rt.TLorentzVector(), array("i", [0])
    ttp4 = rt.TLorentzVector()

    tree.Branch("z1p4", z1p4)
    tree.Branch("z1pdg", z1pdg, "z1pdg/I")
    tree.Branch("z2p4", z2p4)
    tree.Branch("z2pdg", z2pdg, "z2pdg/I")
    tree.Branch("ttp4", z2p4)


    # Leptons (sorted by momentum in Z CM frame)
    l1p4, l1pdg = rt.TLorentzVector(), array("i", [0]) 
    l2p4, l2pdg = rt.TLorentzVector(), array("i", [0]) 
    l3p4, l3pdg = rt.TLorentzVector(), array("i", [0]) 
    l4p4, l4pdg = rt.TLorentzVector(), array("i", [0]) 

    tree.Branch("l1p4", l1p4)
    tree.Branch("l1pdg", l1pdg, "l1pdg/I")
    tree.Branch("l2p4", l2p4)
    tree.Branch("l2pdg", l2pdg, "l2pdg/I")
    tree.Branch("l3p4", l3p4)
    tree.Branch("l3pdg", l3pdg, "l3pdg/I")
    tree.Branch("l4p4", l4p4)
    tree.Branch("l4pdg", l4pdg, "l4pdg/I")


    # Angles
    theta, phi = array("f", [0]), array("f", [0])
    z1alpha, z2alpha = array("f", [0]), array("f", [0])
    z1theta, z2theta = array("f", [0]), array("f", [0])

    tree.Branch("theta", theta, "theta/F")
    tree.Branch("phi", phi, "phi/F")
    tree.Branch("z1alpha", z1alpha, "z1alpha/F")
    tree.Branch("z2alpha", z2alpha, "z2alpha/F")
    tree.Branch("z1theta", z1theta, "z1theta/F")
    tree.Branch("z2theta", z2theta, "z2theta/F")

    p1lpp4, p2lpp4 = rt.TLorentzVector(), rt.TLorentzVector()
    tree.Branch("p1lpp4", p1lpp4)
    tree.Branch("p2lpp4", p2lpp4)


    # "Weight" from alpha_QED 
    weight_alphaQED = array("f", [0])
    tree.Branch("weight_alphaQED", weight_alphaQED, "weight_alphaQED/F")




    ################
    #  EVENT LOOP  #
    ################


    myLHEfile = LHEfile(lheName)
    myLHEfile.setMax(n_max)
    
    print("Looping over", n_max, "events...", end='')
    sys.stdout.flush()

    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:

        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_evts += 1




        ################
        #  CATEGORIZE  #
        ################


        ##  GET LEPTONS  ##

        particles = myLHEevent.Particles
        elecs, muons = [], []
        lep_q, lep_id = {}, {}

        for i in range(0, len(particles)):
            p = particles[i]

            if abs(p['ID']) == 11:
                elecs.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                lep_q[elecs[-1]] = sign(p['ID'])
                lep_id[elecs[-1]] = abs(p['ID'])

            elif abs(p['ID']) == 13:
                muons.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                lep_q[muons[-1]] = sign(p['ID'])
                lep_id[muons[-1]] = abs(p['ID'])

 

        ##  PAIR LEPTONS  ##

        # Combine flavors and sort by momentum
        # (needed to choose pairs in same-flavor case)
        leptons = muons + elecs
        leptons.sort(key=lambda fourvec: fourvec.P(), reverse=True)


        # Determine channel and fill info
        if len(muons) == 4:
#           continue
            n_4m += 1
            z1pdg[0], z2pdg[0] = 13, 13
            pair1, pair2, p1_leps, p2_leps = GetMassPairs(leptons, lep_q)


        elif len(elecs) == 4:
#           continue
            n_4e += 1
            z1pdg[0], z2pdg[0] = 11, 11
            pair1, pair2, p1_leps, p2_leps = GetMassPairs(leptons, lep_q)


        elif len(muons) == 2 and len(elecs) == 2:
            muon_pair, elec_pair = GetSum(muons), GetSum(elecs)

            if muon_pair.M() > elec_pair.M():
#               continue
                n_2m2e += 1
                z1pdg[0], z2pdg[0] = 13, 11
                pair1, pair2 = muon_pair, elec_pair
                p1_leps, p2_leps = muons, elecs

            else:
#               continue
                n_2e2m += 1
                z1pdg[0], z2pdg[0] = 11, 13
                pair1, pair2 = elec_pair, muon_pair
                p1_leps, p2_leps = elecs, muons


        else:
            print("ERROR: wrong number of leptons")



        trio = leptons[1] + leptons[2] + leptons[3]


        p1_leps.sort(key=lambda fourvec: fourvec.P(), reverse=True)
        p2_leps.sort(key=lambda fourvec: fourvec.P(), reverse=True)


        p1_boost, p2_boost = pair1.BoostVector(), pair2.BoostVector()
        b_pair1, b_pair2 = GetBoosted(pair1, p2_boost), GetBoosted(pair2, p1_boost)

        if (lep_q[p1_leps[0]] > 0):
            b_p1plus = GetBoosted(p1_leps[0], p1_boost)
        else:
            b_p1plus = GetBoosted(p1_leps[1], p1_boost)

        if (lep_q[p2_leps[0]] > 0):
            b_p2plus = GetBoosted(p2_leps[0], p2_boost)
        else:
            b_p2plus = GetBoosted(p2_leps[1], p2_boost)




        ###############
        #  FILL TREE  #
        ###############


        ##  WEIGHT  ##

        # ("default"/"calculated" values of 1/alpha)
        weight_alphaQED[0] = aEWM1 / get_alpha_inv(pair2.M())



        ##  OBSERVABLES  ##

        triple_product[0], phi[0], sin_phi[0] = GetTripleProductPairs(p1_leps, p2_leps, lep_q)

        

        ##  MOMENTA  ##

        # (simply writing "p4 = _p4[0]" seems to produce a seg fault...)
        z1p4.SetPxPyPzE(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.E())
        z2p4.SetPxPyPzE(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.E())

        ttp4.SetPxPyPzE(trio.Px(), trio.Py(), trio.Pz(), trio.E())

        l1p4.SetPxPyPzE(leptons[0].Px(), leptons[0].Py(), leptons[0].Pz(), leptons[0].E())
        l2p4.SetPxPyPzE(leptons[1].Px(), leptons[1].Py(), leptons[1].Pz(), leptons[1].E())
        l3p4.SetPxPyPzE(leptons[2].Px(), leptons[2].Py(), leptons[2].Pz(), leptons[2].E())
        l4p4.SetPxPyPzE(leptons[3].Px(), leptons[3].Py(), leptons[3].Pz(), leptons[3].E())


        
        ##  CHARGE, FLAVOR  ##

        # (remember: charge and ID are stored as dictionaries)
        l1pdg[0] = lep_q[leptons[0]] * lep_id[leptons[0]]
        l2pdg[0] = lep_q[leptons[1]] * lep_id[leptons[1]]
        l3pdg[0] = lep_q[leptons[2]] * lep_id[leptons[2]]
        l4pdg[0] = lep_q[leptons[3]] * lep_id[leptons[3]]

        tree.Fill()



        ## ANGLES ##

        theta[0]= p1_leps[1].Angle(pair2.Vect())

        z1alpha[0] = p1_leps[0].Angle(p1_leps[1].Vect())
        z2alpha[0] = p2_leps[0].Angle(p2_leps[1].Vect())

        z1theta[0] = b_p1plus.Angle(b_pair2.Vect())
        z2theta[0] = b_p2plus.Angle(b_pair1.Vect())

        p1lpp4.SetPxPyPzE(b_p1plus.Px(), b_p1plus.Py(), b_p1plus.Pz(), b_p1plus.E())
        p2lpp4.SetPxPyPzE(b_p2plus.Px(), b_p2plus.Py(), b_p2plus.Pz(), b_p2plus.E())




        
        ###############
        #    DEBUG    #
        ###############


        if (debug):
            print(n_evts)
            print(p1_leps[0].P(), lep_q[p1_leps[0]] * lep_id[p1_leps[0]], "\t", p1_leps[1].P(), lep_q[p1_leps[1]] * lep_id[p1_leps[1]])
            print(p2_leps[0].P(), lep_q[p2_leps[0]] * lep_id[p2_leps[0]], "\t", p2_leps[1].P(), lep_q[p2_leps[1]] * lep_id[p2_leps[1]])
            print("")
            sys.stdout.flush()


        if (debug):
#           print(aEWM1 / get_alpha_inv(pair2.M()))
            print(l1p4.P(), l1pdg[0])
            print(l2p4.P(), l2pdg[0])
            print(l3p4.P(), l3pdg[0])
            print(l4p4.P(), l4pdg[0])
            print(z1pdg[0], z2pdg[0])
            print("")
            print("")
            sys.stdout.flush()



        del oneEvent, myLHEevent


    ##############
    #  END LOOP  #
    ##############




    ##  PRINTOUT  ##

    print("done!")
    print("\nFound:")
    print("\t", n_4m, "4mu decays")
    print("\t", n_2m2e, "2mu2e decays")
    print("\t", n_2e2m, "2e2mu decays")
    print("\t", n_4e, "4e decays")
    print("\nWrote events to", outName)




    ###############
    #  SAVE FILE  #
    ###############

    outFile.Write()
    outFile.Close()
