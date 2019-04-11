// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

// Custom
#include "Cuts2017.hh"

using namespace std;


/*
**  DrawDists4l
**
**  Draws UNSCALED distributions for a "boosted_" sample
*/ 

void DrawCanvases4l()
{

    //
    //  OPTIONS
    //

    const TString suffix = "phase_space";


    //
    //  SAMPLE INFO
    //

    const unsigned N = 5;
    unsigned                   L4 = 0,  M4 = 1, ME = 2, EM = 3, E4 = 4;     // Indices
    TString selection[N]    = {"4l",    "4m",   "2m2e", "2e2m", "4e"};
    unsigned chanIdx[N]     = {5,       6,      7,      8,      9};
    TString lepChan[N]      = {_l,      _mu,    _l,     _l,     _e};



    //
    //  OUTPUT FILE
    //

    TString prefix  = "unscaled4l";
    TString outName = prefix + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");



    //
    //  HISTOGRAMS
    //

    vector<tuple<TString, TString, TString, Int_t, Double_t, Double_t>> v = {

        //          name            quantity         axis label             bins    xmin    xmax
//      make_tuple( "channel",      "channel",       "",                    4,      5.5,    9.5),
//      make_tuple( "nPV",          "nPV",           _nPV,                  20,     0,      60),

//      // Lab frame kinematics
        make_tuple( "zzm",          "zzp4.M()",      _m_(_4l),              20,     80,     100),
        make_tuple( "zzpt",         "zzp4.Pt()",     _pT_(_4l),             20,     0,      100),

        make_tuple( "z1m",          "z1p4.M()",      _m_(_Z1),              23,     0,      92),
        make_tuple( "z1pt",         "z1p4.Pt()",     _pT_(_Z1),             20,     0,      120),
//      make_tuple( "z1pdg",        "z1pdg",         "",                    3,      10.5,   13.5),

        make_tuple( "z2m",          "z2p4.M()",      _m_(_Z2),              22,     0,      44),
        make_tuple( "z2pt",         "z2p4.Pt()",     _pT_(_Z2),             20,     0,      60),
//      make_tuple( "z2pdg",        "z2pdg",         "",                    3,      10.5,   13.5),

//      make_tuple( "l1pt",         "l1p4.Pt()",     _pT_(_l_(1)),          24,     0,      120),
//      make_tuple( "l1eta",        "l1p4.Eta()",    _eta_(_l_(1)),         20,     -2.5,   2.5),
//      make_tuple( "l1pdg",        "l1pdg",         "",                    27,     -13.5,  13.5),
//      make_tuple( "l1z",          "l1z",           "",                    2,      0.5,    2.5),

//      make_tuple( "l2pt",         "l2p4.Pt()",     _pT_(_l_(2)),          24,     0,      60),
//      make_tuple( "l2eta",        "l2p4.Eta()",    _eta_(_l_(2)),         20,     -2.5,   2.5),
//      make_tuple( "l2pdg",        "l2pdg",         "",                    27,     -13.5,  13.5),
//      make_tuple( "l2z",          "l2z",           "",                    2,      0.5,    2.5),

//      make_tuple( "l3pt",         "l3p4.Pt()",     _pT_(_l_(3)),          20,     1,      41),
//      make_tuple( "l3eta",        "l3p4.Eta()",    _eta_(_l_(3)),         20,     -2.5,   2.5),
//      make_tuple( "l3pdg",        "l3pdg",         "",                    27,     -13.5,  13.5),
//      make_tuple( "l3z",          "l3z",           "",                    2,      0.5,    2.5),

//      make_tuple( "l4pt",         "l4p4.Pt()",     _pT_(_l_(4)),          25,     0,      25),
//      make_tuple( "l4eta",        "l4p4.Eta()",    _eta_(_l_(4)),         20,     -2.5,   2.5),
//      make_tuple( "l4pdg",        "l4pdg",         "",                    27,     -13.5,  13.5),
//      make_tuple( "l4z",          "l4z",           "",                    2,      0.5,    2.5),
 
        // Z rest frame kinematics                                                
        make_tuple( "b_z1p",        "b_z1p4.P()",    "",                    16,     0,      48),
        make_tuple( "b_ttm",        "b_ttp4.M()",    _m_(_l_("2,3,4")),     15,     0,      60),
                                    
        make_tuple( "b_l1p",        "b_l1v3.Mag()",  _p_(_l_(1)),           13,     24,     50),
        make_tuple( "b_l1pdg",      "b_l1pdg",       "",                    27,     -13.5,  13.5),
//      make_tuple( "b_l1z",        "b_l1z",         "",                    2,      0.5,    2.5),
                                    
        make_tuple( "b_l2p",        "b_l2v3.Mag()",  _p_(_l_(2)),           15,     15,     45),
        make_tuple( "b_l2pdg",      "b_l2pdg",       "",                    27,     -13.5,  13.5),
//      make_tuple( "b_l2z",        "b_l2z",         "",                    2,      0.5,    2.5),
                                    
        make_tuple( "b_l3p",        "b_l3v3.Mag()",  _p_(_l_(3)),           12,     2,      26),
        make_tuple( "b_l3pdg",      "b_l3pdg",       "",                    27,     -13.5,  13.5),
//      make_tuple( "b_l3z",        "b_l3z",         "",                    2,      0.5,    2.5),
                                    
        make_tuple( "b_l4p",        "b_l4v3.Mag()",  _p_(_l_(4)),           14,     0,      21),
        make_tuple( "b_l4pdg",      "b_l4pdg",       "",                    27,     -13.5,  13.5),
//      make_tuple( "b_l4z",        "b_l4z",         "",                    2,      0.5,    2.5),
                        
                        
        // Observables      
//      make_tuple( "psi",              "psi",              _psi,           20,     -5000,  5000),
        make_tuple( "sin_phi",          "sin_phi",          _sinphi,        20,     -1,     1),
//      make_tuple( "sin_phi_2",        "sin_phi",          _sinphi,        2,      -1,     1),
        make_tuple( "cos_theta_z1",     "cos_theta_z1",     _costheta_(_Z1),10,     -1,     1),
        make_tuple( "cos_theta_z2",     "cos_theta_z2",     _costheta_(_Z2),10,     -1,     1),
//      make_tuple( "cos_zeta_z1",      "cos_zeta_z1",      _coszeta_(_Z1), 10,     -1,     1),
//      make_tuple( "cos_zeta_z2",      "cos_zeta_z2",      _coszeta_(_Z2), 10,     -1,     1),
        make_tuple( "angle_z1leps",     "angle_z1leps/3.14",_alpha_(_Z1),   12,     0,      1),
        make_tuple( "angle_z2leps",     "angle_z2leps/3.14",_alpha_(_Z2),   10,     0,      1),
        make_tuple( "angle_z1l2_z2",    "angle_z1l2_z2/3.14",   _beta,      15,     0,      1)
    };



    //
    //  INPUT FILE
    //

    TString inName  = "boosted_" + suffix + ".root";
//  TString inPath  = HOME_PATH + "/Boosted/" + YEAR_STR + "_old/" + inName;
    TFile   *inFile = TFile::Open(inName);

    cout << endl << endl << "Opened " << inName<< endl << endl;


    // Get all "Events" histograms
    vector<TH1*> hEvents;
    TKey *key;
    TIter next(inFile->GetListOfKeys());
    while ((key = (TKey*) next()))
    {
        TString hname = key->GetName();
        if (hname.Contains("Events"))
        {
            TH1 *hist;
            inFile->GetObject(hname, hist);
            hist->SetDirectory(0);
            hEvents.push_back(hist);

            cout << "Got histogram " << hname << endl;
        }
    }






    ////
    ////
    ////    CHANNEL LOOP
    ////
    ////


    for (unsigned i = 0; i < N; i++)
    {
        outFile->mkdir(selection[i]);
        outFile->cd(selection[i]);

        TTree *tree;
        inFile->GetObject(selection[i] + "_" + suffix, tree);

        cout << selection[i] << " tree has " << tree->GetEntries() << " events." << flush;



        //
        //  DRAW
        //

        cout << "\t" << "Drawing histograms..." << flush;

        for (unsigned j = 0; j < v.size(); j++)
        {
            // Get parameters from tuple
            TString hname,  quantity,   xlabel;
            int     bins;
            float   xmin,   xmax;
            TString weight = "weight";
            tie(hname, quantity, xlabel, bins, xmin, xmax) = v[j];


            // Create and draw histogram
            TH1D *h = new TH1D(hname + "_" + suffix, "", 2*bins, xmin, xmax);
            tree->Draw(quantity + ">>+" + hname + "_" + suffix, weight);

            xlabel.ReplaceAll(_l, lepChan[i]);
            h->GetXaxis()->SetTitle(xlabel);
            h->SetStats(0);
            h->SetFillColor(lLightBlue);
            h->SetLineColor(lLightBlue);
            h->SetMinimum(0);
            h->Write();

            TCanvas *c = new TCanvas("c_" + hname + "_" + selection[i], "", 800, 800);
            c->cd();
            h->Draw("HIST");
            c->Write();
        }

        cout << "done!" << endl;

    }

    outFile->cd();
    for (unsigned k = 0; k < hEvents.size(); k++)
        hEvents[k]->Write();
    outFile->Close();
    inFile->Close();

    cout << endl << "Wrote histograms to " << outName << endl << endl << endl;
}
