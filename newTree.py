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


    
#   ##  FIT FUNCTION  ##

#   # (for 1/alpha_QED; used to estimate theo. unc.)

#   # Import "smoothingspline" fit struct from Matlab
#   alpha_inv_p = loadmat("alpha_inv_p.mat")

#   # Recreate as piecewise polynomial function
#   # (takes mass; returns 1/alpha)
#   get_alpha_inv = PPoly(alpha_inv_p['coefs'].transpose(), alpha_inv_p['breaks'][0])

#   # Default MG5 value for 1/alpha_QED
#   aEWM1 = 1.325070e+02



    ##  OUTPUT ROOT FILE  ##

    outName = "4l_fullSM.root"
    outFile = rt.TFile(outName, "RECREATE")
    tree = rt.TTree("tree", "")




    ##############
    #  BRANCHES  #
    ##############

    # (pyROOT requires numbers to be stored as 1D arrays)



    ##  PAIRS  ##

    # (p is high-mass pair, k is low-mass pair)

    _pp_pair,   _p_pdg      = rt.TLorentzVector(),  array("i", [0])
    _kk_pair,   _k_pdg      = rt.TLorentzVector(),  array("i", [0])
    _trio                   = rt.TLorentzVector()

    tree.Branch("pp_pair",  _pp_pair)
    tree.Branch("p_pdg",    _p_pdg,     "p_pdg/I")
    tree.Branch("kk_pair",  _kk_pair)
    tree.Branch("k_pdg",    _k_pdg,     "k_pdg/I")
    tree.Branch("trio",     _trio)


    # Boosted pairs

    _pp_pair_kk,    _kk_pair_pp     = rt.TLorentzVector(),  rt.TLorentzVector()

    tree.Branch("pp_pair_kk",   _pp_pair_kk)
    tree.Branch("kk_pair_pp",   _kk_pair_pp)



    ##  LEPTONS  ##

    # "i" is P-sorted index

    _p_plus,    _p_plus_i   = rt.TLorentzVector(),  array("i", [0])
    _p_minus,   _p_minus_i  = rt.TLorentzVector(),  array("i", [0])
    _k_plus,    _k_plus_i   = rt.TLorentzVector(),  array("i", [0])
    _k_minus,   _k_minus_i  = rt.TLorentzVector(),  array("i", [0])

    tree.Branch("p_plus",       _p_plus)
    tree.Branch("p_plus_i",     _p_plus_i,  "p_plus_i/I")
    tree.Branch("p_minus",      _p_minus)
    tree.Branch("p_minus_i",    _p_minus_i, "p_minus_i/I")
    tree.Branch("k_plus",       _k_plus)
    tree.Branch("k_plus_i",     _k_plus_i,  "k_plus_i/I")
    tree.Branch("k_minus",      _k_minus)
    tree.Branch("k_minus_i",    _k_minus_i, "k_minus_i/I")


    # Boosted leptons

    _p_plus_pp, _k_plus_kk  = rt.TLorentzVector(),  rt.TLorentzVector()

    tree.Branch("p_plus_pp",    _p_plus_pp)
    tree.Branch("k_plus_kk",    _k_plus_kk)



    ##  ANGLES  ##

    _phi,       _p2_kk_angle    = array("f", [0]),  array("f", [0])
    _pp_angle,  _kk_angle       = array("f", [0]),  array("f", [0])
    _p_theta,   _k_theta        = array("f", [0]),  array("f", [0])

    tree.Branch("phi",          _phi,           "phi/F")
    tree.Branch("p2_kk_angle",  _p2_kk_angle,   "p2_kk_angle/F")
    tree.Branch("pp_angle",     _pp_angle,      "pp_angle/F")
    tree.Branch("kk_angle",     _kk_angle,      "kk_angle/F")
    tree.Branch("p_theta",      _p_theta,       "p_theta/F")
    tree.Branch("k_theta",      _k_theta,       "k_theta/F")



    ##  OBSERVABLES  ##

    _psi,       _sin_phi    = array("f", [0]),  array("f", [0])

    tree.Branch("psi",      _psi,       "psi/F")
    tree.Branch("sin_phi",  _sin_phi,   "sin_phi/F")



#   ##  WEIGHT  ##

#   # from alpha_QED 

#   _weight = array("f", [0])
#   tree.Branch("weight", _weight, "weight/F")




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
        charge = {}

        for i in range(0, len(particles)):
            p = particles[i]

            if abs(p['ID']) == 11:
                elecs.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                charge[elecs[-1]] = sign(p['ID'])

            elif abs(p['ID']) == 13:
                muons.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                charge[muons[-1]] = sign(p['ID'])

 

        ##  PAIR LEPTONS  ##

        # 4 muon
        if len(muons) == 4:
#           continue
            n_4m += 1
            p_pdg, k_pdg = 13, 13
            p_plus, p_minus, k_plus, k_minus = GetMassPairs(muons, charge)


        # 4 electron
        elif len(elecs) == 4:
#           continue
            n_4e += 1
            p_pdg, k_pdg = 11, 11
            p_plus, p_minus, k_plus, k_minus = GetMassPairs(elecs, charge)


        # Mixed cases
        elif len(muons) == 2 and len(elecs) == 2:
            muon_pair, elec_pair = GetSum(muons), GetSum(elecs)

            # 2 mu, 2 e
            if muon_pair.M() > elec_pair.M():
#               continue
                n_2m2e += 1
                p_pdg, k_pdg = 13, 11

                if charge[muons[0]] > 0:
                    p_plus, p_minus = muons[0], muons[1]
                else:
                    p_minus, p_plus = muons[0], muons[1]

                if charge[elecs[0]] > 0:
                    k_plus, k_minus = elecs[0], elecs[1]
                else:
                    k_minus, k_plus = elecs[0], elecs[1]


            # 2 e, 2 mu
            else:
#               continue
                n_2e2m += 1
                p_pdg, k_pdg = 11, 13

                if charge[elecs[0]] > 0:
                    p_plus, p_minus = elecs[0], elecs[1]
                else:
                    p_minus, p_plus = elecs[0], elecs[1]

                if charge[muons[0]] > 0:
                    k_plus, k_minus = muons[0], muons[1]
                else:
                    k_minus, k_plus = muons[0], muons[1]


        else:
            print("ERROR: wrong number of leptons")


        pp_pair = p_plus + p_minus
        kk_pair = k_plus + k_minus


        ##  SORT LEPTONS  ##

        leps = muons + elecs
        leps.sort(key=lambda fourvec: fourvec.P(), reverse=True)
        index = {leps[0]:0, leps[1]:1, leps[2]:2, leps[3]:3}

        lead = leps[0]
        trio = leps[1] + leps[2] + leps[3]




        ##################
        #  CALCULATIONS  #
        ##################



        ##  TRIPLE PRODUCT  ##
        
        p_plus_3, p_minus_3 = p_plus.Vect(), p_minus.Vect()
        k_plus_3, k_minus_3 = k_plus.Vect(), k_minus.Vect()

        # Normals to decay planes (TVector3s)
        N_pp, N_kk = p_plus_3.Cross(p_minus_3), k_plus_3.Cross(k_minus_3)
        n_pp, n_kk = N_pp.Unit(), N_kk.Unit()

        # Unit direction vectors
        z_pp, z_kk = pp_pair.Vect().Unit(), kk_pair.Vect().Unit()

        # Observables
        psi = k_plus_3.Dot(N_pp)
        phi = N_pp.Angle(N_kk)                  # 0 < phi < pi          (0 < sin phi < 1)
        sin_phi = n_pp.Cross(n_kk).Dot(z_pp)    # -1 < sin phi < 1      (-pi/2 < phi < pi/2?)



        ##  BOOST  ##

        pp_boost,   kk_boost    = pp_pair.BoostVector(),            kk_pair.BoostVector()

        pp_pair_kk, kk_pair_pp  = GetBoosted(pp_pair, kk_boost),    GetBoosted(kk_pair, pp_boost)
        p_plus_pp,  k_plus_kk   = GetBoosted(p_plus, pp_boost),     GetBoosted(k_plus, kk_boost)



        ##  OTHER ANGLES  ##

        # Find trailing p lepton
        if p_plus.P() > p_minus.P():
            p2 = p_minus
        else:
            p2 = p_plus

        p2_kk_angle = p2.Angle(kk_pair.Vect())


        # Angles between pairs
        pp_angle = p_plus.Angle(p_minus.Vect())
        kk_angle = k_plus.Angle(k_minus.Vect())


        # Pair frame angles
        p_theta = p_plus_pp.Angle(kk_pair_pp.Vect())
        k_theta = k_plus_kk.Angle(pp_pair_kk.Vect())



#       ##  WEIGHT  ##

#       # ("default"/"calculated" values of 1/alpha)

#       weight = aEWM1 / get_alpha_inv(kk_pair.M())




        ###############
        #  FILL TREE  #
        ###############
        

        ##  MOMENTA  ##

        # (simply writing "_p4 = p4" seems to produce a seg fault...)

        _pp_pair.SetPxPyPzE(pp_pair.Px(),   pp_pair.Py(),   pp_pair.Pz(),   pp_pair.E())
        _kk_pair.SetPxPyPzE(kk_pair.Px(),   kk_pair.Py(),   kk_pair.Pz(),   kk_pair.E())
        _trio.SetPxPyPzE(   trio.Px(),      trio.Py(),      trio.Pz(),      trio.E())

        _p_plus.SetPxPyPzE( p_plus.Px(),    p_plus.Py(),    p_plus.Pz(),    p_plus.E())
        _p_minus.SetPxPyPzE(p_minus.Px(),   p_minus.Py(),   p_minus.Pz(),   p_minus.E())
        _k_plus.SetPxPyPzE( k_plus.Px(),    k_plus.Py(),    k_plus.Pz(),    k_plus.E())
        _k_minus.SetPxPyPzE(k_minus.Px(),   k_minus.Py(),   k_minus.Pz(),   k_minus.E())

        _p_plus_pp.SetPxPyPzE(p_plus_pp.Px(), p_plus_pp.Py(), p_plus_pp.Pz(), p_plus_pp.E())
        _k_plus_kk.SetPxPyPzE(k_plus_kk.Px(), k_plus_kk.Py(), k_plus_kk.Pz(), k_plus_kk.E())

        _pp_pair_kk.SetPxPyPzE(pp_pair_kk.Px(), pp_pair_kk.Py(), pp_pair_kk.Pz(), pp_pair_kk.E())
        _kk_pair_pp.SetPxPyPzE(kk_pair_pp.Px(), kk_pair_pp.Py(), kk_pair_pp.Pz(), kk_pair_pp.E())

        
        ##  FLAVOR  ##

        _p_pdg[0],  _k_pdg[0]   = p_pdg,    k_pdg


        
        ##  INDEX  ##
        
        _p_plus_i[0],   _p_minus_i[0]    = index[p_plus],    index[p_minus]
        _k_plus_i[0],   _k_minus_i[0]    = index[k_plus],    index[k_minus]



        ## ANGLES ##

        _phi[0],        _p2_kk_angle[0] = phi,          p2_kk_angle
        _pp_angle[0],   _kk_angle[0]    = pp_angle,     kk_angle
        _p_theta[0],    _k_theta[0]     = p_theta,      k_theta



        ##  OBSERVABLES  ##

        _psi[0],    _sin_phi[0] = psi,  sin_phi



#       ##  WEIGHT  ##

#       _weight[0] = weight




        tree.Fill()

        


        ###############
        #    DEBUG    #
        ###############


        if (debug):
            print("Event", n_evts)
            print("p:",     pp_pair.M(),    "\tid",     p_pdg)
            print("p+:",    p_plus.P(),     "\tq",      charge[p_plus])
            print("p-:",    p_minus.P(),    "\tq",      charge[p_minus])
            print("")
            print("k:",     kk_pair.M(),    "\tid",     k_pdg)
            print("k+:",    k_plus.P(),     "\tq",      charge[k_plus])
            print("k-:",    k_minus.P(),    "\tq",      charge[k_minus])
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
    print("\t", n_4m,   "4mu decays")
    print("\t", n_2m2e, "2mu2e decays")
    print("\t", n_2e2m, "2e2mu decays")
    print("\t", n_4e,   "4e decays")
    print("\nWrote events to", outName)




    ###############
    #  SAVE FILE  #
    ###############

    outFile.Write()
    outFile.Close()
