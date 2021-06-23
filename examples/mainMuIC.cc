// main36.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; DIS;

// Basic setup for Deeply Inelastic Scattering at HERA.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "TF1.h"
#include "TRandom3.h"

using namespace Pythia8;

int main() {

  // Beam energies, minimal Q2, number of events to generate.
  double eProton   = 275.;
  double eMuon = 960.;
  double Q2min     = 1.0;

  double mProton = 0.938272;
  double pProton = sqrt(eProton*eProton - mProton*mProton);
  double mMuon = 0.105658;
  double pMuon = sqrt(eMuon*eMuon - mMuon*mMuon);
  double roots = sqrt((eMuon+eProton)*(eMuon+eProton)-(pMuon-pProton)*(pMuon-pProton));
//  int    nEvent    = 100000;
  int    nEvent    = 50000;

  double muAngleRes = 0.0002;
  TF1* muAngleResFunc = new TF1("muAngleResFunc","exp(-0.5*x*x/[0]/[0])",-10.*muAngleRes,10.*muAngleRes);  
  muAngleResFunc->SetParameter(0,muAngleRes);

  double muPRelRes = 0.01;
  TF1* muPRelResFunc = new TF1("muPRelResFunc","exp(-0.5*x*x/([0]*[0]+0.0001*0.0001*[1]*[1]))",-10,10);
  muPRelResFunc->SetParameter(0,muPRelRes);

  double chTrkPRelRes = 0.01;
  TF1* chTrkPRelResFunc = new TF1("chTrkPRelResFunc","exp(-0.5*x*x/([0]*[0]+0.001*0.001*[1]*[1]))",-10,10);
  chTrkPRelResFunc->SetParameter(0,chTrkPRelRes);

  double chTrkAngleRes = 0.0002;
  TF1* chTrkAngleResFunc = new TF1("chTrkAngleResFunc","exp(-0.5*x*x/([0]*[0]+0.002*0.002/[1]/[1]))",-0.01,0.01);
  chTrkAngleResFunc->SetParameter(0,chTrkAngleRes);
         
  double caloAngleRes = 0.087/sqrt(12);
  TF1* caloAngleResFunc = new TF1("caloAngleResFunc","exp(-0.5*x*x/[0]/[0])",-10.*caloAngleRes,10.*caloAngleRes);
  caloAngleResFunc->SetParameter(0,caloAngleRes);
 
  double emERelRes = 0.02;
  TF1* emERelResFunc = new TF1("emERelResFunc","exp(-0.5*x*x/([0]*[0]+0.01/[1]))",-10,10);
  emERelResFunc->SetParameter(0,emERelRes);

  double hadERelRes = 0.1;
  TF1* hadERelResFunc = new TF1("hadERelResFunc","exp(-0.5*x*x/([0]*[0]+0.25/[1]))",-10,10);
  hadERelResFunc->SetParameter(0,hadERelRes);

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = muon.
  pythia.readString("Beams:idB = 13");
  pythia.settings.parm("Beams:eB", eMuon);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Uncomment to allow charged current.
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Initialize.
  pythia.init();

  // Histograms.
  double Wmax = sqrt(4.* eProton * eMuon);
  TH1D* Nchhist = new TH1D("Nch","N_{ch}", 500, 0., 500.);
  TH1D* Qhist = new TH1D("Qhist","Q [GeV]", 2000, 0., 1000.);
  TH1D* Whist = new TH1D("Whist","W [GeV]", 200, 0., Wmax);
  TH1D* xhist = new TH1D("xhist","x", 1000000, 0., 1.);
  TH1D* yhist = new TH1D("yhist","y", 100, 0., 1.);
  TH1D* pTehist = new TH1D("pTehist","pT of scattered muon [GeV]", 200, 0., 500.);
  TH1D* pTrhist = new TH1D("pTrhist","pT of radiated parton [GeV]", 200, 0., 500.);
  TH1D* pTdhist = new TH1D("pTdhist","ratio pT_parton/pT_muon", 100, 0., 5.);
  TH2D* Q2xhist = new TH2D("Q2xhist",";x;Q^{2} [GeV]", 100000, 0.000005, 1., 200, 1., 200.);
  TH2D* Q2xhist_eta = new TH2D("Q2xhist_eta",";x;Q^{2} [GeV]", 100000, 0.000005, 1., 200, 1., 200.);

  TH2D* petahadhist = new TH2D("petahadhist",";#eta;p [GeV]", 80, -8., 8., 800, 0., 200.);
  TH2D* petamuhist = new TH2D("petamuhist",";#eta;p [GeV]", 80, -8., 8., 100, 0., 1.0);

  TH2D* hQ2resVsQ2 = new TH2D("hQ2resVsQ2",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsx = new TH2D("hQ2resVsx",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsy = new TH2D("hQ2resVsy",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2 = new TH2D("hxresVsQ2",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsx = new TH2D("hxresVsx",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsy = new TH2D("hxresVsy",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2 = new TH2D("hyresVsQ2",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsx = new TH2D("hyresVsx",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsy = new TH2D("hyresVsy",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2JB = new TH2D("hQ2resVsQ2JB",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxJB = new TH2D("hQ2resVsxJB",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyJB = new TH2D("hQ2resVsyJB",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2JB = new TH2D("hxresVsQ2JB",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxJB = new TH2D("hxresVsxJB",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyJB = new TH2D("hxresVsyJB",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2JB = new TH2D("hyresVsQ2JB",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxJB = new TH2D("hyresVsxJB",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyJB = new TH2D("hyresVsyJB",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2JB4 = new TH2D("hQ2resVsQ2JB4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxJB4 = new TH2D("hQ2resVsxJB4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyJB4 = new TH2D("hQ2resVsyJB4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2JB4 = new TH2D("hxresVsQ2JB4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxJB4 = new TH2D("hxresVsxJB4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyJB4 = new TH2D("hxresVsyJB4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2JB4 = new TH2D("hyresVsQ2JB4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxJB4 = new TH2D("hyresVsxJB4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyJB4 = new TH2D("hyresVsyJB4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2JB5 = new TH2D("hQ2resVsQ2JB5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxJB5 = new TH2D("hQ2resVsxJB5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyJB5 = new TH2D("hQ2resVsyJB5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2JB5 = new TH2D("hxresVsQ2JB5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxJB5 = new TH2D("hxresVsxJB5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyJB5 = new TH2D("hxresVsyJB5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2JB5 = new TH2D("hyresVsQ2JB5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxJB5 = new TH2D("hyresVsxJB5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyJB5 = new TH2D("hyresVsyJB5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2DA = new TH2D("hQ2resVsQ2DA",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxDA = new TH2D("hQ2resVsxDA",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyDA = new TH2D("hQ2resVsyDA",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2DA = new TH2D("hxresVsQ2DA",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxDA = new TH2D("hxresVsxDA",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyDA = new TH2D("hxresVsyDA",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2DA = new TH2D("hyresVsQ2DA",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxDA = new TH2D("hyresVsxDA",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyDA = new TH2D("hyresVsyDA",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2DA4 = new TH2D("hQ2resVsQ2DA4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxDA4 = new TH2D("hQ2resVsxDA4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyDA4 = new TH2D("hQ2resVsyDA4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2DA4 = new TH2D("hxresVsQ2DA4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxDA4 = new TH2D("hxresVsxDA4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyDA4 = new TH2D("hxresVsyDA4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2DA4 = new TH2D("hyresVsQ2DA4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxDA4 = new TH2D("hyresVsxDA4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyDA4 = new TH2D("hyresVsyDA4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2DA5 = new TH2D("hQ2resVsQ2DA5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxDA5 = new TH2D("hQ2resVsxDA5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyDA5 = new TH2D("hQ2resVsyDA5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2DA5 = new TH2D("hxresVsQ2DA5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxDA5 = new TH2D("hxresVsxDA5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyDA5 = new TH2D("hxresVsyDA5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2DA5 = new TH2D("hyresVsQ2DA5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxDA5 = new TH2D("hyresVsxDA5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyDA5 = new TH2D("hyresVsyDA5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A2 = new TH2D("hQ2resVsQ2A2",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA2 = new TH2D("hQ2resVsxA2",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA2 = new TH2D("hQ2resVsyA2",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A2 = new TH2D("hxresVsQ2A2",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA2 = new TH2D("hxresVsxA2",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA2 = new TH2D("hxresVsyA2",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A2 = new TH2D("hyresVsQ2A2",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA2 = new TH2D("hyresVsxA2",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA2 = new TH2D("hyresVsyA2",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A2Eta4 = new TH2D("hQ2resVsQ2A2Eta4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA2Eta4 = new TH2D("hQ2resVsxA2Eta4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA2Eta4 = new TH2D("hQ2resVsyA2Eta4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A2Eta4 = new TH2D("hxresVsQ2A2Eta4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA2Eta4 = new TH2D("hxresVsxA2Eta4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA2Eta4 = new TH2D("hxresVsyA2Eta4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A2Eta4 = new TH2D("hyresVsQ2A2Eta4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA2Eta4 = new TH2D("hyresVsxA2Eta4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA2Eta4 = new TH2D("hyresVsyA2Eta4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A2Eta5 = new TH2D("hQ2resVsQ2A2Eta5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA2Eta5 = new TH2D("hQ2resVsxA2Eta5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA2Eta5 = new TH2D("hQ2resVsyA2Eta5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A2Eta5 = new TH2D("hxresVsQ2A2Eta5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA2Eta5 = new TH2D("hxresVsxA2Eta5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA2Eta5 = new TH2D("hxresVsyA2Eta5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A2Eta5 = new TH2D("hyresVsQ2A2Eta5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA2Eta5 = new TH2D("hyresVsxA2Eta5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA2Eta5 = new TH2D("hyresVsyA2Eta5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A4 = new TH2D("hQ2resVsQ2A4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA4 = new TH2D("hQ2resVsxA4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA4 = new TH2D("hQ2resVsyA4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A4 = new TH2D("hxresVsQ2A4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA4 = new TH2D("hxresVsxA4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA4 = new TH2D("hxresVsyA4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A4 = new TH2D("hyresVsQ2A4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA4 = new TH2D("hyresVsxA4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA4 = new TH2D("hyresVsyA4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A4Eta4 = new TH2D("hQ2resVsQ2A4Eta4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA4Eta4 = new TH2D("hQ2resVsxA4Eta4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA4Eta4 = new TH2D("hQ2resVsyA4Eta4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A4Eta4 = new TH2D("hxresVsQ2A4Eta4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA4Eta4 = new TH2D("hxresVsxA4Eta4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA4Eta4 = new TH2D("hxresVsyA4Eta4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A4Eta4 = new TH2D("hyresVsQ2A4Eta4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA4Eta4 = new TH2D("hyresVsxA4Eta4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA4Eta4 = new TH2D("hyresVsyA4Eta4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A4Eta5 = new TH2D("hQ2resVsQ2A4Eta5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA4Eta5 = new TH2D("hQ2resVsxA4Eta5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA4Eta5 = new TH2D("hQ2resVsyA4Eta5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A4Eta5 = new TH2D("hxresVsQ2A4Eta5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA4Eta5 = new TH2D("hxresVsxA4Eta5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA4Eta5 = new TH2D("hxresVsyA4Eta5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A4Eta5 = new TH2D("hyresVsQ2A4Eta5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA4Eta5 = new TH2D("hyresVsxA4Eta5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA4Eta5 = new TH2D("hyresVsyA4Eta5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A5 = new TH2D("hQ2resVsQ2A5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA5 = new TH2D("hQ2resVsxA5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA5 = new TH2D("hQ2resVsyA5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A5 = new TH2D("hxresVsQ2A5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA5 = new TH2D("hxresVsxA5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA5 = new TH2D("hxresVsyA5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A5 = new TH2D("hyresVsQ2A5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA5 = new TH2D("hyresVsxA5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA5 = new TH2D("hyresVsyA5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A5Eta4 = new TH2D("hQ2resVsQ2A5Eta4",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA5Eta4 = new TH2D("hQ2resVsxA5Eta4",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA5Eta4 = new TH2D("hQ2resVsyA5Eta4",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A5Eta4 = new TH2D("hxresVsQ2A5Eta4",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA5Eta4 = new TH2D("hxresVsxA5Eta4",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA5Eta4 = new TH2D("hxresVsyA5Eta4",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A5Eta4 = new TH2D("hyresVsQ2A5Eta4",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA5Eta4 = new TH2D("hyresVsxA5Eta4",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA5Eta4 = new TH2D("hyresVsyA5Eta4",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  TH2D* hQ2resVsQ2A5Eta5 = new TH2D("hQ2resVsQ2A5Eta5",";Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hQ2resVsxA5Eta5 = new TH2D("hQ2resVsxA5Eta5",";x;#DeltaQ^{2}/Q^{2}", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hQ2resVsyA5Eta5 = new TH2D("hQ2resVsyA5Eta5",";y;#DeltaQ^{2}/Q^{2}", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hxresVsQ2A5Eta5 = new TH2D("hxresVsQ2A5Eta5",";Q^{2} [GeV];#Deltax/x", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hxresVsxA5Eta5 = new TH2D("hxresVsxA5Eta5",";x;#Deltax/x", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hxresVsyA5Eta5 = new TH2D("hxresVsyA5Eta5",";y;#Deltax/x", 100, 0., 1., 40, -1.5, 2.5);
  TH2D* hyresVsQ2A5Eta5 = new TH2D("hyresVsQ2A5Eta5",";Q^{2} [GeV];#Deltay/y", 100000, 1., 400000., 40, -1.5, 2.5);
  TH2D* hyresVsxA5Eta5 = new TH2D("hyresVsxA5Eta5",";x;#Deltay/y", 100000, 0.000005, 1., 40, -1.5, 2.5);
  TH2D* hyresVsyA5Eta5 = new TH2D("hyresVsyA5Eta5",";y;#Deltay/y", 100, 0., 1., 40, -1.5, 2.5);

  
  double xbins[] = {0.000005,0.00001,0.00002,0.00003,0.00004,0.00005,0.00007,0.0001,0.0002,0.0003,0.0004,0.0005,0.0007,0.001,0.002,0.003,0.004,0.005,0.007,0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  double q2bins[] = {1,2,3,4,5,8,10,20,30,40,50,80,100,200,300,400,500,800,1000,2000,3000,4000,5000,8000,10000,20000,30000,40000,50000,80000,100000,200000,300000,400000,500000};
  double resbins[] = {-1.5,-1,-0.7,-0.5,-0.3,-0.2,-0.15,-0.1,-0.07,-0.05,-0.03,-0.01,0,0.01,0.03,0.05,0.07,0.1,0.15,0.2,0.3,0.5,0.7,1.0,1.5,2.0,2.5,3.0};
  TH3D* hQ2resVsQ2x = new TH3D("hQ2resVsQ2x",";x;Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 34,xbins,34,q2bins,27,resbins);
  TH3D* hxresVsQ2x = new TH3D("hxresVsQ2x",";x;Q^{2} [GeV];#Deltax/x", 34,xbins,34,q2bins,27,resbins);
  TH3D* hyresVsQ2x = new TH3D("hyresVsQ2x",";x;Q^{2} [GeV];#Deltay/y", 34,xbins,34,q2bins,27,resbins);

  TH3D* hQ2resVsQ2xJB5 = new TH3D("hQ2resVsQ2xJB5",";x;Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 34,xbins,34,q2bins,27,resbins);
  TH3D* hxresVsQ2xJB5 = new TH3D("hxresVsQ2xJB5",";x;Q^{2} [GeV];#Deltax/x", 34,xbins,34,q2bins,27,resbins);
  TH3D* hyresVsQ2xJB5 = new TH3D("hyresVsQ2xJB5",";x;Q^{2} [GeV];#Deltay/y", 34,xbins,34,q2bins,27,resbins);

  TH3D* hQ2resVsQ2xDA5 = new TH3D("hQ2resVsQ2xDA5",";x;Q^{2} [GeV];#DeltaQ^{2}/Q^{2}", 34,xbins,34,q2bins,27,resbins);
  TH3D* hxresVsQ2xDA5 = new TH3D("hxresVsQ2xDA5",";x;Q^{2} [GeV];#Deltax/x", 34,xbins,34,q2bins,27,resbins);
  TH3D* hyresVsQ2xDA5 = new TH3D("hyresVsQ2xDA5",";x;Q^{2} [GeV];#Deltay/y", 34,xbins,34,q2bins,27,resbins);

  // Set up the ROOT TFile and TTree.
//  TFile *file = TFile::Open("pytree_mup_ycut001.root","recreate");
  TFile *file = TFile::Open("pytree_mup_muic_dp05_dang2_Q21.root","recreate");
//  TFile *file = TFile::Open("pytree_mup_muic.root","recreate");
//  Event *event = &pythia.event;
  TTree *T = new TTree("T","ev1 Tree");
  T->Branch("event",&event);

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Four-momenta of proton, muon, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn    = event[4].p();
    Vec4 peOut   = event[6].p();
    Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);

    double Q2rec = 4.*peIn.e()*peOut.e()*cos(peOut.theta()/2.)*cos(peOut.theta()/2.);
    double yrec = 1-peOut.e()*(1-cos(peOut.theta()))/2./peIn.e();
    double xrec = Q2rec/roots/roots/yrec;

    muPRelResFunc->SetParameter(1,peOut.pAbs());
    muAngleResFunc->SetParameter(1,peOut.pAbs());

//cout<<peOut.pAbs()<<" "<<muAngleResFunc->GetRandom()<<" "<<muPRelResFunc->GetRandom()<<endl;

    double muThetaSmear = peOut.theta() + muAngleResFunc->GetRandom();
    double muPSmear = peOut.pAbs() * (1 + muPRelResFunc->GetRandom());
    double muESmear = sqrt(muPSmear*muPSmear + mMuon*mMuon);
    double Q2recsm = 4*peIn.e()*muESmear*cos(muThetaSmear/2.)*cos(muThetaSmear/2.);
    double yrecsm = 1-muESmear*(1-cos(muThetaSmear))/2./peIn.e();
    double xrecsm = Q2recsm/roots/roots/yrecsm;

    // Fill kinematics histograms.
    Qhist->Fill( sqrt(Q2) );
    Whist->Fill( sqrt(W2) );
    xhist->Fill( x );
    yhist->Fill( y );
    pTehist->Fill( event[6].pT() );
    Q2xhist->Fill( x, Q2 );
    if(fabs(event[6].eta())<6) Q2xhist_eta->Fill( x, Q2 );

    if(yrec>0.999 || yrec<0.001) continue;

    hQ2resVsQ2->Fill(Q2rec,Q2recsm/Q2rec-1);
    hQ2resVsx->Fill(xrec,Q2recsm/Q2rec-1);
    hQ2resVsy->Fill(yrec,Q2recsm/Q2rec-1);
    hxresVsQ2->Fill(Q2rec,xrecsm/xrec-1);
    hxresVsx->Fill(xrec,xrecsm/xrec-1);
    hxresVsy->Fill(yrec,xrecsm/xrec-1);
    hyresVsQ2->Fill(Q2rec,yrecsm/yrec-1);
    hyresVsx->Fill(xrec,yrecsm/yrec-1);
    hyresVsy->Fill(yrec,yrecsm/yrec-1);

    hQ2resVsQ2x->Fill(xrec,Q2rec,Q2recsm/Q2rec-1);
    hxresVsQ2x->Fill(xrec,Q2rec,xrecsm/xrec-1);
    hyresVsQ2x->Fill(xrec,Q2rec,yrecsm/yrec-1);

/*
if(x<0.01) continue;
cout<<peOut.pAbs()<<" "<<muPSmear<<endl;
cout<<peOut.theta()<<" "<<muThetaSmear<<endl;
cout<<Q2<<" "<<Q2rec<<" "<<Q2recsm<<endl;
cout<<y<<" "<<yrec<<" "<<yrecsm<<endl;
cout<<x<<" "<<xrec<<" "<<xrecsm<<endl;
cout<<endl;
*/
    petamuhist->Fill( event[6].eta(), event[6].pT()*cosh(event[6].eta()) );

    // pT spectrum of partons being radiated in shower.
    int nMult = 0;
    double sumPx=0;
    double sumPy=0;
    double sumEPz=0;

    double sumPx4=0;
    double sumPy4=0;
    double sumEPz4=0;

    double sumPx5=0;
    double sumPy5=0;
    double sumEPz5=0;

    for (int i = 0; i < event.size(); ++i)
    {
      if (event[i].statusAbs() == 43) {
        pTrhist->Fill( event[i].pT() );
        pTdhist->Fill( event[i].pT() / event[6].pT() );
      }
//      if (event[i].status()>=81 && event[i].status()<=89) {
//        if(event[i].charge() == 0) continue;
//      if(fabs(event[i].id())==211 || fabs(event[i].id())==111 || fabs(event[i].id())==321 || fabs(event[i].id())==311 || fabs(event[i].id())==2112 || fabs(event[i].id())==2212 || fabs(event[i].id())==22 || fabs(event[i].id())==11) {
      if(event[i].isFinal() && event[i].id()!=13) {

        double pp = event[i].pAbs();
        double ppx = event[i].px();
        double ppy = event[i].py();
        double ppz = event[i].pz();
        double pe = event[i].e();
        double ppt = event[i].pT();
        double pphi = event[i].phi();
        double peta = event[i].eta();
        double pmass = event[i].m();

        // charged particle
        if(fabs(event[i].id())==211 || fabs(event[i].id())==321 || fabs(event[i].id())==2212 || fabs(event[i].id())==11)
        {
/*
          ppx *= (1 + hadPRelResFunc->GetRandom()); 
          ppy *= (1 + hadPRelResFunc->GetRandom());
          ppz *= (1 + hadPRelResFunc->GetRandom());
*/
          chTrkPRelResFunc->SetParameter(1,pp);

          pp *= (1 + chTrkPRelResFunc->GetRandom());
          pphi += chTrkAngleResFunc->GetRandom();
          peta += chTrkAngleResFunc->GetRandom();
          ppx = pp/cosh(peta)*cos(pphi);
          ppy = pp/cosh(peta)*sin(pphi);
          ppz = pp/cosh(peta)*sinh(peta);
          pe = sqrt(pmass*pmass+ppx*ppx+ppy*ppy+ppz*ppz);
        }
        // photon and electron
        if(fabs(event[i].id())==22)
        {
          emERelResFunc->SetParameter(1,pe);

          pe *= (1 + emERelResFunc->GetRandom());
          pphi += caloAngleResFunc->GetRandom();
          peta += caloAngleResFunc->GetRandom(); 
          ppx = pe/cosh(peta)*cos(pphi);
          ppy = pe/cosh(peta)*sin(pphi);
          ppz = pe/cosh(peta)*sinh(peta);
        }
        // neutrons
        if(fabs(event[i].id())==2112)
        {
          hadERelResFunc->SetParameter(1,pe);

          pe *= (1 + hadERelResFunc->GetRandom());
          pphi += caloAngleResFunc->GetRandom();
          peta += caloAngleResFunc->GetRandom();
          ppx = pe/cosh(peta)*cos(pphi);
          ppy = pe/cosh(peta)*sin(pphi);
          ppz = pe/cosh(peta)*sinh(peta);          
        }

        sumPx += ppx;
        sumPy += ppy;
        sumEPz += pe-ppz;

        if(fabs(event[i].eta())<4)
        {
          sumPx4 += ppx;
          sumPy4 += ppy; 
          sumEPz4 += pe-ppz;
        }

        if(fabs(event[i].eta())<5)
        {
          sumPx5 += ppx;  
          sumPy5 += ppy;
          sumEPz5 += pe-ppz;
        }

        if(fabs(event[i].eta())<5) nMult++;
        petahadhist->Fill( event[i].eta(), event[i].pT()*cosh(event[i].eta()) );
      }
    }
    double sumPt2 = sumPx*sumPx+sumPy*sumPy;
    double yjb = sumEPz/2.0/peIn.e();
    double Q2jb = sumPt2/(1-yjb);
    double xjb = Q2jb/yjb/roots/roots;
    double agljb = acos((sumPt2-sumEPz*sumEPz)/(sumPt2+sumEPz*sumEPz));
    double Q2DA = 4.*peIn.e()*peIn.e()*sin(agljb)*(cos(muThetaSmear)+1)/(sin(agljb)+sin(muThetaSmear)-sin(agljb+muThetaSmear));
    double yDA = sin(muThetaSmear)*(1-cos(agljb))/(sin(agljb)+sin(muThetaSmear)-sin(agljb+muThetaSmear));
    double xDA = Q2DA/yDA/roots/roots;
//    double yA2 = 0.5*(1-muESmear/2.0/peIn.e())*(1-cos(agljb))
//    double Q2A2 = 4.*peIn.e()*peIn.e()
    double yA4 = (muESmear-peIn.e())/(muESmear+sqrt(sumPt2)/sin(agljb)-2.0*peIn.e());
    double Q2A4 = 4.*peIn.e()*peIn.e()*(peIn.e()-sqrt(sumPt2)/sin(agljb))/(muESmear+sqrt(sumPt2)/sin(agljb)-2.0*peIn.e())+4.*peIn.e()*muESmear;
    double xA4 = Q2A4/yA4/roots/roots;

    double sumPt25 = sumPx5*sumPx5+sumPy5*sumPy5;
    double yjb5 = sumEPz5/2.0/peIn.e();
    double Q2jb5 = sumPt25/(1-yjb5);
    double xjb5 = Q2jb5/yjb5/roots/roots;
    double agljb5 = acos((sumPt25-sumEPz5*sumEPz5)/(sumPt25+sumEPz5*sumEPz5));
    double Q2DA5 = 4.*peIn.e()*peIn.e()*sin(agljb5)*(cos(muThetaSmear)+1)/(sin(agljb5)+sin(muThetaSmear)-sin(agljb5+muThetaSmear));
    double yDA5 = sin(muThetaSmear)*(1-cos(agljb5))/(sin(agljb5)+sin(muThetaSmear)-sin(agljb5+muThetaSmear));
    double xDA5 = Q2DA5/yDA5/roots/roots;
    double yA4Eta5 = (muESmear-peIn.e())/(muESmear+sqrt(sumPt25)/sin(agljb5)-2.0*peIn.e());
    double Q2A4Eta5 = 4.*peIn.e()*peIn.e()*(peIn.e()-sqrt(sumPt25)/sin(agljb5))/(muESmear+sqrt(sumPt25)/sin(agljb5)-2.0*peIn.e())+4.*peIn.e()*muESmear;
    double xA4Eta5 = Q2A4Eta5/yA4Eta5/roots/roots;

    double sumPt24 = sumPx4*sumPx4+sumPy4*sumPy4;
    double yjb4 = sumEPz4/2.0/peIn.e();
    double Q2jb4 = sumPt24/(1-yjb4);
    double xjb4 = Q2jb4/yjb4/roots/roots;
    double agljb4 = acos((sumPt24-sumEPz4*sumEPz4)/(sumPt24+sumEPz4*sumEPz4));
    double Q2DA4 = 4.*peIn.e()*peIn.e()*sin(agljb4)*(cos(muThetaSmear)+1)/(sin(agljb4)+sin(muThetaSmear)-sin(agljb4+muThetaSmear));
    double yDA4 = sin(muThetaSmear)*(1-cos(agljb4))/(sin(agljb4)+sin(muThetaSmear)-sin(agljb4+muThetaSmear));
    double xDA4 = Q2DA4/yDA4/roots/roots;
    double yA4Eta4 = (muESmear-peIn.e())/(muESmear+sqrt(sumPt24)/sin(agljb4)-2.0*peIn.e());
    double Q2A4Eta4 = 4.*peIn.e()*peIn.e()*(peIn.e()-sqrt(sumPt24)/sin(agljb4))/(muESmear+sqrt(sumPt24)/sin(agljb4)-2.0*peIn.e())+4.*peIn.e()*muESmear;
    double xA4Eta4 = Q2A4Eta4/yA4Eta4/roots/roots;
/*
cout<<endl;
cout<<sumPx<<" "<<sumPy<<" "<<sumEPz<<endl;
cout<<peOut.pAbs()<<" "<<muPSmear<<endl;
cout<<peOut.theta()<<" "<<muThetaSmear<<endl;
cout<<Q2<<" "<<Q2rec<<" "<<Q2recsm<<" "<<Q2jb<<" "<<Q2jb4<<" "<<Q2jb5<<endl;
cout<<y<<" "<<yrec<<" "<<yrecsm<<" "<<yjb<<" "<<yjb4<<" "<<yjb5<<endl;
cout<<x<<" "<<xrec<<" "<<xrecsm<<" "<<xjb<<" "<<xjb4<<" "<<xjb5<<endl;
cout<<endl;
*/
/*
    hQ2resVsQ2JB->Fill(Q2jb,Q2jb/Q2rec-1);
    hQ2resVsxJB->Fill(xjb,Q2jb/Q2rec-1);
    hQ2resVsyJB->Fill(yjb,Q2jb/Q2rec-1);
    hxresVsQ2JB->Fill(Q2jb,xjb/xrec-1);
    hxresVsxJB->Fill(xjb,xjb/xrec-1);
    hxresVsyJB->Fill(yjb,xjb/xrec-1);
    hyresVsQ2JB->Fill(Q2jb,yjb/yrec-1);
    hyresVsxJB->Fill(xjb,yjb/yrec-1);
    hyresVsyJB->Fill(yjb,yjb/yrec-1);

    hQ2resVsQ2JB4->Fill(Q2jb4,Q2jb4/Q2rec-1);
    hQ2resVsxJB4->Fill(xjb4,Q2jb4/Q2rec-1);
    hQ2resVsyJB4->Fill(yjb4,Q2jb4/Q2rec-1);
    hxresVsQ2JB4->Fill(Q2jb4,xjb4/xrec-1);
    hxresVsxJB4->Fill(xjb4,xjb4/xrec-1);
    hxresVsyJB4->Fill(yjb4,xjb4/xrec-1);
    hyresVsQ2JB4->Fill(Q2jb4,yjb4/yrec-1);
    hyresVsxJB4->Fill(xjb4,yjb4/yrec-1);
    hyresVsyJB4->Fill(yjb4,yjb4/yrec-1);

    hQ2resVsQ2JB5->Fill(Q2jb5,Q2jb5/Q2rec-1);
    hQ2resVsxJB5->Fill(xjb5,Q2jb5/Q2rec-1);
    hQ2resVsyJB5->Fill(yjb5,Q2jb5/Q2rec-1);
    hxresVsQ2JB5->Fill(Q2jb5,xjb5/xrec-1);
    hxresVsxJB5->Fill(xjb5,xjb5/xrec-1);
    hxresVsyJB5->Fill(yjb5,xjb5/xrec-1);
    hyresVsQ2JB5->Fill(Q2jb5,yjb5/yrec-1);
    hyresVsxJB5->Fill(xjb5,yjb5/yrec-1);
    hyresVsyJB5->Fill(yjb5,yjb5/yrec-1);

    hQ2resVsQ2DA->Fill(Q2DA,Q2DA/Q2rec-1);
    hQ2resVsxDA->Fill(xDA,Q2DA/Q2rec-1);
    hQ2resVsyDA->Fill(yDA,Q2DA/Q2rec-1);
    hxresVsQ2DA->Fill(Q2DA,xDA/xrec-1);
    hxresVsxDA->Fill(xDA,xDA/xrec-1);
    hxresVsyDA->Fill(yDA,xDA/xrec-1);
    hyresVsQ2DA->Fill(Q2DA,yDA/yrec-1);
    hyresVsxDA->Fill(xDA,yDA/yrec-1);
    hyresVsyDA->Fill(yDA,yDA/yrec-1);

    hQ2resVsQ2DA4->Fill(Q2DA4,Q2DA4/Q2rec-1);
    hQ2resVsxDA4->Fill(xDA4,Q2DA4/Q2rec-1);
    hQ2resVsyDA4->Fill(yDA4,Q2DA4/Q2rec-1);
    hxresVsQ2DA4->Fill(Q2DA4,xDA4/xrec-1);
    hxresVsxDA4->Fill(xDA4,xDA4/xrec-1);
    hxresVsyDA4->Fill(yDA4,xDA4/xrec-1);
    hyresVsQ2DA4->Fill(Q2DA4,yDA4/yrec-1);
    hyresVsxDA4->Fill(xDA4,yDA4/yrec-1);
    hyresVsyDA4->Fill(yDA4,yDA4/yrec-1);

    hQ2resVsQ2DA5->Fill(Q2DA5,Q2DA5/Q2rec-1);
    hQ2resVsxDA5->Fill(xDA5,Q2DA5/Q2rec-1);
    hQ2resVsyDA5->Fill(yDA5,Q2DA5/Q2rec-1);
    hxresVsQ2DA5->Fill(Q2DA5,xDA5/xrec-1);
    hxresVsxDA5->Fill(xDA5,xDA5/xrec-1);
    hxresVsyDA5->Fill(yDA5,xDA5/xrec-1);
    hyresVsQ2DA5->Fill(Q2DA5,yDA5/yrec-1);
    hyresVsxDA5->Fill(xDA5,yDA5/yrec-1);
    hyresVsyDA5->Fill(yDA5,yDA5/yrec-1);
*/


    hQ2resVsQ2JB->Fill(Q2rec,Q2jb/Q2rec-1);
    hQ2resVsxJB->Fill(xrec,Q2jb/Q2rec-1);
    hQ2resVsyJB->Fill(yrec,Q2jb/Q2rec-1);
    hxresVsQ2JB->Fill(Q2rec,xjb/xrec-1);
    hxresVsxJB->Fill(xrec,xjb/xrec-1);
    hxresVsyJB->Fill(yrec,xjb/xrec-1);
    hyresVsQ2JB->Fill(Q2rec,yjb/yrec-1);
    hyresVsxJB->Fill(xrec,yjb/yrec-1);
    hyresVsyJB->Fill(yrec,yjb/yrec-1);

    hQ2resVsQ2JB4->Fill(Q2rec,Q2jb4/Q2rec-1);
    hQ2resVsxJB4->Fill(xrec,Q2jb4/Q2rec-1);
    hQ2resVsyJB4->Fill(yrec,Q2jb4/Q2rec-1);
    hxresVsQ2JB4->Fill(Q2rec,xjb4/xrec-1);
    hxresVsxJB4->Fill(xrec,xjb4/xrec-1);
    hxresVsyJB4->Fill(yrec,xjb4/xrec-1);
    hyresVsQ2JB4->Fill(Q2rec,yjb4/yrec-1);
    hyresVsxJB4->Fill(xrec,yjb4/yrec-1);
    hyresVsyJB4->Fill(yrec,yjb4/yrec-1);

    hQ2resVsQ2JB5->Fill(Q2rec,Q2jb5/Q2rec-1);
    hQ2resVsxJB5->Fill(xrec,Q2jb5/Q2rec-1);
    hQ2resVsyJB5->Fill(yrec,Q2jb5/Q2rec-1);
    hxresVsQ2JB5->Fill(Q2rec,xjb5/xrec-1);
    hxresVsxJB5->Fill(xrec,xjb5/xrec-1);
    hxresVsyJB5->Fill(yrec,xjb5/xrec-1);
    hyresVsQ2JB5->Fill(Q2rec,yjb5/yrec-1);
    hyresVsxJB5->Fill(xrec,yjb5/yrec-1);
    hyresVsyJB5->Fill(yrec,yjb5/yrec-1);

    hQ2resVsQ2DA->Fill(Q2rec,Q2DA/Q2rec-1);
    hQ2resVsxDA->Fill(xrec,Q2DA/Q2rec-1);
    hQ2resVsyDA->Fill(yrec,Q2DA/Q2rec-1);
    hxresVsQ2DA->Fill(Q2rec,xDA/xrec-1);
    hxresVsxDA->Fill(xrec,xDA/xrec-1);
    hxresVsyDA->Fill(yrec,xDA/xrec-1);
    hyresVsQ2DA->Fill(Q2rec,yDA/yrec-1);
    hyresVsxDA->Fill(xrec,yDA/yrec-1);
    hyresVsyDA->Fill(yrec,yDA/yrec-1);

    hQ2resVsQ2DA4->Fill(Q2rec,Q2DA4/Q2rec-1);
    hQ2resVsxDA4->Fill(xrec,Q2DA4/Q2rec-1);
    hQ2resVsyDA4->Fill(yrec,Q2DA4/Q2rec-1);
    hxresVsQ2DA4->Fill(Q2rec,xDA4/xrec-1);
    hxresVsxDA4->Fill(xrec,xDA4/xrec-1);
    hxresVsyDA4->Fill(yrec,xDA4/xrec-1);
    hyresVsQ2DA4->Fill(Q2rec,yDA4/yrec-1);
    hyresVsxDA4->Fill(xrec,yDA4/yrec-1);
    hyresVsyDA4->Fill(yrec,yDA4/yrec-1);

    hQ2resVsQ2DA5->Fill(Q2rec,Q2DA5/Q2rec-1);
    hQ2resVsxDA5->Fill(xrec,Q2DA5/Q2rec-1);
    hQ2resVsyDA5->Fill(yrec,Q2DA5/Q2rec-1);
    hxresVsQ2DA5->Fill(Q2rec,xDA5/xrec-1);
    hxresVsxDA5->Fill(xrec,xDA5/xrec-1);
    hxresVsyDA5->Fill(yrec,xDA5/xrec-1);
    hyresVsQ2DA5->Fill(Q2rec,yDA5/yrec-1);
    hyresVsxDA5->Fill(xrec,yDA5/yrec-1);
    hyresVsyDA5->Fill(yrec,yDA5/yrec-1);

    hQ2resVsQ2A4->Fill(Q2rec,Q2A4/Q2rec-1);
    hQ2resVsxA4->Fill(xrec,Q2A4/Q2rec-1);
    hQ2resVsyA4->Fill(yrec,Q2A4/Q2rec-1);
    hxresVsQ2A4->Fill(Q2rec,xA4/xrec-1);
    hxresVsxA4->Fill(xrec,xA4/xrec-1);
    hxresVsyA4->Fill(yrec,xA4/xrec-1);
    hyresVsQ2A4->Fill(Q2rec,yA4/yrec-1);
    hyresVsxA4->Fill(xrec,yA4/yrec-1);
    hyresVsyA4->Fill(yrec,yA4/yrec-1);

    hQ2resVsQ2A4Eta4->Fill(Q2rec,Q2A4Eta4/Q2rec-1);
    hQ2resVsxA4Eta4->Fill(xrec,Q2A4Eta4/Q2rec-1);
    hQ2resVsyA4Eta4->Fill(yrec,Q2A4Eta4/Q2rec-1);
    hxresVsQ2A4Eta4->Fill(Q2rec,xA4Eta4/xrec-1);
    hxresVsxA4Eta4->Fill(xrec,xA4Eta4/xrec-1);
    hxresVsyA4Eta4->Fill(yrec,xA4Eta4/xrec-1);
    hyresVsQ2A4Eta4->Fill(Q2rec,yA4Eta4/yrec-1);
    hyresVsxA4Eta4->Fill(xrec,yA4Eta4/yrec-1);
    hyresVsyA4Eta4->Fill(yrec,yA4Eta4/yrec-1);

    hQ2resVsQ2A4Eta5->Fill(Q2rec,Q2A4Eta5/Q2rec-1);
    hQ2resVsxA4Eta5->Fill(xrec,Q2A4Eta5/Q2rec-1);
    hQ2resVsyA4Eta5->Fill(yrec,Q2A4Eta5/Q2rec-1);
    hxresVsQ2A4Eta5->Fill(Q2rec,xA4Eta5/xrec-1);
    hxresVsxA4Eta5->Fill(xrec,xA4Eta5/xrec-1);
    hxresVsyA4Eta5->Fill(yrec,xA4Eta5/xrec-1);
    hyresVsQ2A4Eta5->Fill(Q2rec,yA4Eta5/yrec-1);
    hyresVsxA4Eta5->Fill(xrec,yA4Eta5/yrec-1);
    hyresVsyA4Eta5->Fill(yrec,yA4Eta5/yrec-1);

    hQ2resVsQ2xJB5->Fill(xrec,Q2rec,Q2jb5/Q2rec-1);
    hxresVsQ2xJB5->Fill(xrec,Q2rec,xjb5/xrec-1);
    hyresVsQ2xJB5->Fill(xrec,Q2rec,yjb5/yrec-1);

    hQ2resVsQ2xDA5->Fill(xrec,Q2rec,Q2DA5/Q2rec-1);
    hxresVsQ2xDA5->Fill(xrec,Q2rec,xDA5/xrec-1);
    hyresVsQ2xDA5->Fill(xrec,Q2rec,yDA5/yrec-1);

    Nchhist->Fill(nMult);

    T->Fill();

  // End of event loop. Statistics and histograms.
  }
  pythia.stat();
  cout << Qhist << Whist << xhist << yhist << pTehist << pTrhist << pTdhist;

  //  Write tree.
  T->Print();
  T->Write();
  Nchhist->Write();
  Qhist->Write();
  Whist->Write();
  xhist->Write();
  yhist->Write();
  pTehist->Write();
  pTrhist->Write();
  pTdhist->Write();
  Q2xhist->Write();
  Q2xhist_eta->Write();
  petahadhist->Write();
  petamuhist->Write();

  hQ2resVsQ2->Write();
  hQ2resVsx->Write();
  hQ2resVsy->Write();
  hxresVsQ2->Write();
  hxresVsx->Write();
  hxresVsy->Write();
  hyresVsQ2->Write();
  hyresVsx->Write();
  hyresVsy->Write();

  hQ2resVsQ2JB->Write();
  hQ2resVsxJB->Write();
  hQ2resVsyJB->Write();
  hxresVsQ2JB->Write();
  hxresVsxJB->Write();
  hxresVsyJB->Write();
  hyresVsQ2JB->Write();
  hyresVsxJB->Write();
  hyresVsyJB->Write();

  hQ2resVsQ2JB4->Write();
  hQ2resVsxJB4->Write();
  hQ2resVsyJB4->Write();
  hxresVsQ2JB4->Write();
  hxresVsxJB4->Write();
  hxresVsyJB4->Write();
  hyresVsQ2JB4->Write();
  hyresVsxJB4->Write();
  hyresVsyJB4->Write();

  hQ2resVsQ2JB5->Write();
  hQ2resVsxJB5->Write();
  hQ2resVsyJB5->Write();
  hxresVsQ2JB5->Write();
  hxresVsxJB5->Write();
  hxresVsyJB5->Write();
  hyresVsQ2JB5->Write();
  hyresVsxJB5->Write();
  hyresVsyJB5->Write();

  hQ2resVsQ2DA->Write();
  hQ2resVsxDA->Write();
  hQ2resVsyDA->Write();
  hxresVsQ2DA->Write();
  hxresVsxDA->Write();
  hxresVsyDA->Write();
  hyresVsQ2DA->Write();
  hyresVsxDA->Write();
  hyresVsyDA->Write();

  hQ2resVsQ2DA4->Write();
  hQ2resVsxDA4->Write();
  hQ2resVsyDA4->Write();
  hxresVsQ2DA4->Write();
  hxresVsxDA4->Write();
  hxresVsyDA4->Write();
  hyresVsQ2DA4->Write();
  hyresVsxDA4->Write();
  hyresVsyDA4->Write();

  hQ2resVsQ2DA5->Write();
  hQ2resVsxDA5->Write();
  hQ2resVsyDA5->Write();
  hxresVsQ2DA5->Write();
  hxresVsxDA5->Write();
  hxresVsyDA5->Write();
  hyresVsQ2DA5->Write();
  hyresVsxDA5->Write();
  hyresVsyDA5->Write();

  hQ2resVsQ2A4->Write();
  hQ2resVsxA4->Write();
  hQ2resVsyA4->Write();
  hxresVsQ2A4->Write();
  hxresVsxA4->Write();
  hxresVsyA4->Write();
  hyresVsQ2A4->Write();
  hyresVsxA4->Write();
  hyresVsyA4->Write();

  hQ2resVsQ2A4Eta4->Write();
  hQ2resVsxA4Eta4->Write();
  hQ2resVsyA4Eta4->Write();
  hxresVsQ2A4Eta4->Write();
  hxresVsxA4Eta4->Write();
  hxresVsyA4Eta4->Write();
  hyresVsQ2A4Eta4->Write();
  hyresVsxA4Eta4->Write();
  hyresVsyA4Eta4->Write();

  hQ2resVsQ2A4Eta5->Write();
  hQ2resVsxA4Eta5->Write();
  hQ2resVsyA4Eta5->Write();
  hxresVsQ2A4Eta5->Write();
  hxresVsxA4Eta5->Write();
  hxresVsyA4Eta5->Write();
  hyresVsQ2A4Eta5->Write();
  hyresVsxA4Eta5->Write();
  hyresVsyA4Eta5->Write();

  hQ2resVsQ2x->Write();
  hxresVsQ2x->Write();
  hyresVsQ2x->Write();

  hQ2resVsQ2xJB5->Write();
  hxresVsQ2xJB5->Write();
  hyresVsQ2xJB5->Write();

  hQ2resVsQ2xDA5->Write();
  hxresVsQ2xDA5->Write();
  hyresVsQ2xDA5->Write();

  delete file;

  // Done.
  return 0;
}
