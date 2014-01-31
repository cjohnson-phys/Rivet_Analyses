#! /usr/bin/env python

# James Henderson 01/10/12
# Updated to Rivet 2.0.0 on 22/01/14

# This little script will write out MC_VBF.cc to the custom_rivet folder with specific flags for cuts etc.
# This takes a while to compile, so use sparingly if possible

## ATLAS Evgen seems to be confused with units so all "TeV"s in this code really mean "GeV"!!!

import sys, os, copy, subprocess
from optparse import OptionParser
import time

parser = OptionParser()
# default values are already = 0 but, this makes thing clearer for me
parser.add_option("-m", action="store_true", dest="muon",     help="Do analysis for muons rather than default electrons",       default = 0)
parser.add_option("-t", action="store_true", dest="tau",      help="Do analysis for taus rather than default electrons",        default = 0)
parser.add_option("-c", action="store_true", dest="cuts",     help="Apply the standard VBF type cuts",                          default = 0)
parser.add_option("-w", action="store_true", dest="forceW",   help="Force the analysis to only consider VBF W",                 default = 0)
parser.add_option("-z", action="store_true", dest="forceZ",   help="Force the analysis to only consider VBF Z",                 default = 0)
parser.add_option("-s", action="store_true", dest="shape",    help="Don't consider the cross-section, just do shape analysis",  default = 0)
parser.add_option("-e", action="store_true", dest="ByEta",    help="Order the Jets by eta rather than by pT",                   default = 0) 
parser.add_option("-u", action="store_true", dest="units",    help="Change the units to be in GeV rather than TeV default",     default = 0)
parser.add_option("--diM", dest="dijetMassCut", type=int,     help="Put a veto on events that have a dijet mass of < X GeV",    default = 0)

(options, args) = parser.parse_args()

# This is the location of Jim's custom rivet analysis - please alter as approriate
f = open('/Users/Alex/Desktop/HEP/Rivet_Analyses/MC_VBF/MC_VBF.cc', 'w')

unit = "TeV"
if options.units: unit = "GeV"

section ="""// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Analysis.hh"

double PassedWeight = 0;
double TotWeight = 0;

namespace Rivet {
  class MC_VBF : public MC_JetAnalysis {
  public:
  /// Default constructor
    MC_VBF()
      :  MC_JetAnalysis("MC_VBF", 4, "Jets")
    {    }

    void init() {

      FinalState fs;
      WFinder wfinder(fs, -2.5, 2.5, 25.0*"""+unit+""", """
f.write(section)

if options.muon:
    f.write( "PID::MUON, ")
elif options.tau:
    f.write( "PID:TAU, ")
else:
    f.write( "PID::ELECTRON, ")
    pass

section = """60.0*"""+unit+""", 100.0*"""+unit+""", 25.0*"""+unit+""", 0.2);
  addProjection(wfinder, "WFinder");

  ZFinder zfinder(fs, -2.5, 2.5, 25.0*"""+unit+""", """

f.write(section)

if options.muon:
    f.write( "PID::MUON, ")
elif options.tau:
    f.write( "PID:TAU, ")
else:
    f.write( "PID::ELECTRON, ")
    pass

section = """65.0*"""+unit+""", 115.0*"""+unit+""", 0.2, true, true);
  addProjection(zfinder, "ZFinder");

  // The jets are build from 'everything but the boson' in the final state so we need different jet collections depending on what boson we consider
  FastJets jets_z( zfinder.remainingFinalState() , FastJets::ANTIKT, 0.4 );
  jets_z.useInvisibles();
  addProjection(jets_z, "Jets_z");

  FastJets jets_w( wfinder.remainingFinalState() , FastJets::ANTIKT, 0.4);
  jets_w.useInvisibles();
  addProjection(jets_w, "Jets_w");

  //lepton histos
  _h_lepton_phi        = bookHisto1D("Lepton_phi",50,0,4);
  _h_lepton_eta        = bookHisto1D("Lepton_eta",50,-3,3);
  _h_lepton_pT         = bookHisto1D("Lepton_pT",50,0,1000);
  _h_neutrino_phi      = bookHisto1D("Neutrino_phi",50,0,4);
  _h_neutrino_eta      = bookHisto1D("Neutrino_eta",50,-3,3);
  _h_neutrino_pT       = bookHisto1D("Neutrino_pT",50,0,1000);
  
  _h_met               = bookHisto1D("Missing_eT", 30, 0., 1000.);
  _h_met_phi           = bookHisto1D("Missing_eT_phi", 30, 0, 4);
  _h_MT_dilep          = bookHisto1D("Dilepton_mT", 30, 0., 400.);
  _h_mll               = bookHisto1D("Dilepton_mass", 30, 0, 500);
  _h_boson_pT          = bookHisto1D("boson_pT", 30, 0, 2000);
  _h_RATIOboson_pT     = bookHisto1D("RATIOboson_pT", 10, 0, 1000);
  _h_ALTboson_pT       = bookHisto1D("boson_pT_50", 50, 0, 1000);      

  //jet histos
  _h_jet1_phi          = bookHisto1D("Jet1_phi",30,0,4);
  _h_jet1_eta          = bookHisto1D("Jet1_eta",30,-7,7);
  _h_jet1_pT           = bookHisto1D("Jet1_pT",30,0,1000);
  _h_jet2_phi          = bookHisto1D("Jet2_phi",30,0,4);
  _h_jet2_eta          = bookHisto1D("Jet2_eta",30,-7,7);
  _h_jet2_pT           = bookHisto1D("Jet2_pT",30,0,1000);
  _h_jet3_phi          = bookHisto1D("Jet3_phi",30,0,4);
  _h_jet3_eta          = bookHisto1D("Jet3_eta",30,-7,7);
  _h_jet3_pT           = bookHisto1D("Jet3_pT",30,0,1000);
  _h_mjj               = bookHisto1D("Dijet_Mass",30, 0, 6000);
  _h_ALTmjj            = bookHisto1D("hmj1j2_w",200, 0, 5000);
  _h_zoomed_mjj        = bookHisto1D("Dijet_Mass_zoomed",30,0,200);
  _h_RATIOmjj          = bookHisto1D("RATIODijet_Mass",10, 0, 6000);
  _h_scalarSumPt       = bookHisto1D("Scalar_sum_pT",30,0,3000);
  _h_vectorSumPt       = bookHisto1D("Vector_sum_pT",30,0,600);
  _h_deltaEta          = bookHisto1D("Jet_delta_rapidity",30,-10,10);
  _h_delatPhi          = bookHisto1D("Jet_delta_Phi",30,0,4);
  _h_deltaR            = bookHisto1D("Jet_delta_R",30,0,10);
  _h_etaProduct        = bookHisto1D("Eta_product",30,-20,20);
  _h_dijetPt           = bookHisto1D("Dijet_pT",30,0,2000); 
  _h_ALTdijet_pT       = bookHisto1D("Dijet_pT_50",50,0,1000);
  _h_RATIOdijetPt      = bookHisto1D("RATIODijet_pT",10,0,1000); 
  _h_gapref            = bookHisto1D("GapRef_PtAvg",30,0,1000); 
  _h_OppHemi           = bookHisto1D("Opposite_Hemisphere",5,-1,3);
  _h_njet_inclusive    = bookHisto1D("njet_inc",11,-0.5,10.5);

  //Other histos
  _h_MT                = bookHisto1D("Total_mT", 30, 0., 2500.);
  _h_deltaPhiJetLep    = bookHisto1D("Delta_Phi_dijet_dilep",30,-7,7); 
  _h_totVisibleM       = bookHisto1D("Total_visible_Mass",30,0,500); 
  _h_totVBFM           = bookHisto1D("Total_VBF_Mass",30,0,2000); 
  _h_delabosonJetM     = bookHisto1D("Delta_boson_Jet_Mass",30,0,4); 
  _h_ClosestBosM       = bookHisto1D("Closest_to_boson_Mass",30,75,95); 
  _h_dijetEtaStar      = bookHisto1D("Dijet_Eta*",30,0,10);
  _h_WorZ              = bookHisto1D("W_or_Z", 6, -0.5, 5.5);
  _h_dR_lepJet         = bookHisto1D("dRLepJet",30,0,10);
  _h_PID               = bookHisto1D("PID", 26, -0.5, 25.5);
  _h_MomFac            = bookHisto1D("MomFacTot", 31, 0, 10);
  _h_MomFac350         = bookHisto1D("MomFacTot350", 31, 0, 10);
  _h_MomFac1000        = bookHisto1D("MomFacTot1000", 31, 0, 10);
  _h_MomFacX           = bookHisto1D("MomFacX",41,-0.05,4.05);
  _h_MomFacY           = bookHisto1D("MomFacY",41,-0.05,4.05);
  _h_Weight            = bookHisto1D("Weights",30,-0.5,100.5);
  _h_PtBal             = bookHisto1D("PtJ_mB_Balance", 30, 0, 100); 
  _h_Zeppenfeld        = bookHisto1D("Zeppenfeld",30,-7,7);  
  _h_BDPtSum           = bookHisto1D("BosonDijetPtSum",31,-100,100);
  _h_BDPtSum350        = bookHisto1D("BosonDijetPtSum350",31,-100,100);
  _h_BDPtSum1000       = bookHisto1D("BosonDijetPtSum1000",31,-100,100); 
  _h_BDPtSum2Jets      = bookHisto1D("BosonDijetPtSum2Jets",31,-100,100);
  _h_DPtInW            = bookHisto1D("DPtInW",31,-100,100);
  _h_DPtInW350         = bookHisto1D("DPtInW350",31,-100,100);
  _h_DPtInW1000        = bookHisto1D("DPtInW1000",31,-100,100);

  _h_VBFAngle          = bookHisto1D("VBFAngle",50, -2, 2 );
  _h_VBFAngle_cuts     = bookHisto1D("VBFAngle_cuts",50, -2, 2 );
  _h_VBFAngle_B        = bookHisto1D("VBFAngle_B",50, -2, 2 );
  _h_VBFAngle_cuts_B   = bookHisto1D("VBFAngle_cuts_B",50, -2, 2 );
  }

   
  void analyze(const Event& event) {

   const double weight = event.weight();
  
   TotWeight+=weight;
   Jets all_jets;
   FourMomentum boson, lepton, neutrino, second_lepton;

   const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
   const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");
      
   bool isZ;
   if ( zfinder.bosons().size() == 1 && wfinder.bosons().size() == 0 ) {
      isZ = true;
      all_jets      = applyProjection<FastJets>(event, "Jets_z").jetsByPt(30.0*"""+unit+""");
      boson         = zfinder.bosons().front().momentum();
      lepton        = zfinder.constituents()[0].momentum();
      second_lepton = zfinder.constituents()[1].momentum();
   }
   else if ( zfinder.bosons().size() == 0 && wfinder.bosons().size() == 1 ){
      isZ = false;
      all_jets = applyProjection<FastJets>(event, "Jets_w").jetsByPt(30.0*"""+unit+""");
      boson    = wfinder.bosons().front().momentum();
      lepton   = wfinder.constituentLeptons().front().momentum();
      neutrino = wfinder.constituentNeutrinos().front().momentum();
   }
   else 
      vetoEvent;

   """
f.write(section)

if options.forceW:
    section = """
   if ( isZ ) vetoEvent;
   """
    f.write(section)
elif options.forceZ:
    section = """
    if ( !isZ ) vetoEvent;
    """
    f.write(section)
    pass

section = """
   FourMomentum jet1, jet2, jet3;
   Jets my_jets;
   foreach (const Jet& jet, all_jets ){
     if ( fabs( jet.eta() ) > 4.4 ) continue;
     if ( fabs( deltaR( jet, lepton ) ) < 0.3 ) continue;
     if ( isZ ) if ( fabs( deltaR( jet, second_lepton ) ) < 0.3 ) continue;
     my_jets.push_back( jet );
   }

   if ( my_jets.size() < 2 ) vetoEvent;

   jet1 = my_jets[0].momentum();
   jet2 = my_jets[1].momentum();
   if ( my_jets.size() > 2 ) jet3 = my_jets[2].momentum();
   FourMomentum dijet = jet1 + jet2;
   """
f.write(section)

if options.dijetMassCut:
    section="""
    if ( sqrt(fabs(dijet.mass2()))/"""+unit+""" < """+str(options.dijetMassCut)+""" ) vetoEvent;
    """
    f.write(section)
    pass

if options.cuts:
    cuts_section="""
    
    if ( jet1.pT()/"""+unit+""" < 80 )                                                 vetoEvent;
    if ( jet2.pT()/"""+unit+""" < 60 )                                                 vetoEvent;
    if ( lepton.pT()/"""+unit+""" < 25  )                                              vetoEvent;
    // TRANSVERSE MASS CUT!!
    if ( sqrt(fabs(dijet.mass2()))/"""+unit+""" < 500 )                                vetoEvent;
    """

    f.write(cuts_section)
    pass
    
section="""
      // From this point we assume we have a VBF type event - Please be careful with cuts above to ensure this!
            
      foreach (const GenParticle* gp, particles(event.genEvent()))
	_h_PID->fill(abs(gp->pdg_id()), weight);
     
      _h_Weight->fill( weight, 1.0);
      
      //Fill Jets Histograms
      _h_jet1_phi->fill( acos(cos(jet1.phi())), weight);
      _h_jet1_eta->fill( jet1.eta(), weight);
      _h_jet1_pT->fill(  jet1.pT()/"""+unit+""", weight);
      _h_jet2_phi->fill( acos(cos(jet2.phi())), weight);
      _h_jet2_eta->fill( jet2.eta(), weight);
      _h_jet2_pT->fill(  jet2.pT() /"""+unit+""", weight);
 
      if ( my_jets.size() > 2 ){
	_h_jet3_phi->fill( acos(cos(jet3.phi())), weight);
	_h_jet3_eta->fill( jet3.eta(), weight);
	_h_jet3_pT->fill(  jet3.pT() /"""+unit+""", weight);
        _h_Zeppenfeld->fill( (jet3.eta() - ( jet1.eta() + jet2.eta() ) / 2.0)   , weight);
      }

      _h_mjj->fill( sqrt(fabs(dijet.mass2())) /"""+unit+""", weight);
      _h_ALTmjj->fill( sqrt(fabs(dijet.mass2())) /"""+unit+""", weight);
      _h_zoomed_mjj->fill( sqrt(fabs(dijet.mass2())) /"""+unit+""", weight);
      _h_RATIOmjj->fill( sqrt(fabs(dijet.mass2())) /"""+unit+""", weight);
      _h_scalarSumPt->fill( ( jet1.pT() + jet2.pT() ) /"""+unit+""", weight);
      _h_vectorSumPt->fill( dijet.pT() /"""+unit+""", weight);
      _h_etaProduct->fill( (jet1.eta() * jet2.eta()), weight);
      _h_dijetPt->fill( dijet.pT() /"""+unit+""", weight);
      _h_RATIOdijetPt->fill( dijet.pT() /"""+unit+""", weight);
      _h_ALTdijet_pT->fill( dijet.pT() /"""+unit+""", weight);
      _h_gapref->fill( ((jet1.pT() + jet2.pT())/2)  /"""+unit+""", weight);

      if(jet1.eta()*jet2.eta()<0) _h_OppHemi->fill(2, weight);
      else _h_OppHemi->fill(0, weight);
    
      double deltaEta = jet1.eta() - jet2.eta();
      double deltaPhi = acos(cos(jet1.phi() - jet2.phi() ));
      _h_deltaEta->fill( deltaEta, weight);
      _h_delatPhi->fill( fabs(deltaPhi) , weight);
      _h_deltaR->fill( deltaR( jet1, jet2) , weight);
      _h_dR_lepJet->fill( deltaR( lepton, jet1 ) , weight);
      double mT_dilep = sqrt( boson.Et2() - boson.pT2() ) /"""+unit+""";	
      FourMomentum finalstate = (boson + jet1 + jet2);

      for ( uint n_jet = 0; n_jet <= my_jets.size() ; n_jet++ )
	_h_njet_inclusive->fill(n_jet, weight);
            
      double MomX = fabs( (boson.px() - dijet.px() ) / ( dijet.px() ) ); 
      double MomY = fabs( (boson.py() - dijet.py() ) / ( dijet.py() ) );
      
      _h_BDPtSum->fill( (boson.px() + dijet.px()) /"""+unit+""" + (boson.py() + dijet.py()) /"""+unit+""", weight);
      _h_MomFac->fill( MomX * MomY, weight);
      _h_DPtInW->fill(boson.pT()/"""+unit+""" + dijet.pT()/"""+unit+""" *cos( boson.phi() - dijet.phi() ), weight);
      if( sqrt(dijet.mass2()) /"""+unit+""" > 350){
        _h_BDPtSum350->fill( (boson.px() + dijet.px())/"""+unit+""" + (boson.py() + dijet.py())/"""+unit+""", weight);
        _h_MomFac350->fill( MomX * MomY, weight);
        _h_DPtInW350->fill(boson.pT()/"""+unit+""" + dijet.pT()/"""+unit+""" *cos( boson.phi() - dijet.phi() ), weight);
        if ( my_jets.size() == 2) _h_BDPtSum2Jets->fill( (boson.px()/"""+unit+""" + dijet.px()/"""+unit+""") + (boson.py()/"""+unit+""" + dijet.py()/"""+unit+"""),  weight);  
        if( sqrt(dijet.mass2()) /"""+unit+""" > 1000){
          _h_BDPtSum1000->fill( (boson.px()/"""+unit+""" + dijet.px()/"""+unit+""") + (boson.py()/"""+unit+""" + dijet.py()/"""+unit+"""), weight);
          _h_MomFac1000->fill( MomX * MomY, weight);
          _h_DPtInW1000->fill(boson.pT()/"""+unit+""" + dijet.pT()/"""+unit+"""*cos( boson.phi() - dijet.phi() ) , weight);
        }      
      }
      
      _h_MomFacX->fill( MomX, weight);
      _h_MomFacY->fill( MomY, weight);

      double vbf_angle = (lepton.eta()-((my_jets.at(0).momentum().eta()+my_jets.at(1).momentum().eta())/2.0))/(fabs(my_jets.at(0).momentum().eta()-my_jets.at(1).momentum().eta()));
      double vbf_angle_B = (boson.eta()-((my_jets.at(0).momentum().eta()+my_jets.at(1).momentum().eta())/2.0))/(fabs(my_jets.at(0).momentum().eta()-my_jets.at(1).momentum().eta()));

      _h_VBFAngle->fill( vbf_angle , weight );
      _h_VBFAngle_B->fill( vbf_angle_B , weight );

      if ( sqrt(fabs(dijet.mass2())) /"""+unit+""" > 500 && my_jets.at(0).momentum().pT() /"""+unit+""" > 80 && my_jets.at(1).momentum().pT() /"""+unit+""" > 60 ){
          _h_VBFAngle_cuts->fill( vbf_angle , weight );
          _h_VBFAngle_cuts_B->fill( vbf_angle_B , weight );
      }
      
      //Fill Lepton Histograms
      if ( isZ ) _h_WorZ->fill( 2, weight);
      else if ( wfinder.constituentLeptons().front().charge() > 0 )	_h_WorZ->fill( 4, weight);
      else _h_WorZ->fill( 0, weight);
      
      _h_lepton_phi->fill( acos(cos(lepton.phi())), weight);        
      _h_lepton_eta->fill(lepton.eta(), weight);        
      _h_lepton_pT->fill(lepton.pT()/"""+unit+""", weight);
      
      if ( isZ ) {
      _h_neutrino_phi->fill( acos(cos(second_lepton.phi())), weight);
      _h_neutrino_eta->fill(second_lepton.eta(), weight);
      _h_neutrino_pT->fill(second_lepton.pT() /"""+unit+""", weight);
      }
      else {
        _h_neutrino_phi->fill( acos(cos(neutrino.phi())), weight); 
        _h_neutrino_eta->fill(neutrino.eta(), weight);
        _h_neutrino_pT->fill(neutrino.pT() /"""+unit+""", weight);
      }
      
      _h_MT->fill( sqrt( finalstate.Et2() - finalstate.pT2() ) /"""+unit+""", weight); 	
      _h_MT_dilep->fill(mT_dilep, weight);
     
      if ( !isZ ){
        _h_met->fill( neutrino.Et() /"""+unit+""", weight);
        _h_met_phi->fill( acos(cos(neutrino.phi())), weight);
      }

      _h_PtBal->fill( fabs(dijet.pT() - boson.pT()) /"""+unit+""", weight);
      _h_mll->fill( sqrt(fabs(boson.mass2())) /"""+unit+""", weight);
      _h_boson_pT->fill( boson.pT() /"""+unit+""", weight);
      _h_RATIOboson_pT->fill( boson.pT() /"""+unit+""", weight);
      _h_ALTboson_pT->fill( boson.pT() /"""+unit+""", weight);
	
      _h_deltaPhiJetLep->fill( dijet.phi() - boson.phi() , weight);
      
      double tot_jet_mass = 0;
      for( unsigned int i = 0 ; i<my_jets.size()-1 ; i++)
        tot_jet_mass += sqrt(fabs(my_jets[i].momentum().mass2())) /"""+unit+""";
      _h_totVisibleM->fill( tot_jet_mass, weight);
        
      _h_totVBFM->fill( sqrt(fabs((dijet+boson).mass2())) /"""+unit+""", weight);
      _h_dijetEtaStar->fill( fabs( (jet1.eta() + jet2.eta())/2 - boson.eta() ), weight);
	
      const double Wmass = 80.385;
      const double Zmass = 91.1876;
      double massdiff = 1000000;
      double minMassDiff = 10000000.0;
      double currentMass = 0.;
      double minMass = 0.;
      double boson_mass= Wmass;

      if (isZ) boson_mass = Zmass;
     
      for (uint i=0;i<my_jets.size();i++) {
	FourMomentum jeti = my_jets[i].momentum();

	for (uint j=i+1;j<my_jets.size();j++) {
	  currentMass = 0;
	  massdiff = 10000;      
	  FourMomentum jetj = my_jets[j].momentum();
	  currentMass = sqrt(fabs((jeti + jetj).mass2())) /"""+unit+""";
	  massdiff = fabs(boson_mass - currentMass);

	  if ( massdiff < minMassDiff) {
	    //this is the minimum pairing found so far
	    minMassDiff = massdiff;
	    minMass = currentMass;
	  }
	}//end second loop
      }//end first loop
     
      _h_delabosonJetM->fill( minMassDiff, weight);
      _h_ClosestBosM->fill( minMass, weight);
           
      PassedWeight += weight;
    }

    void finalize() {

    """
f.write(section)

if options.shape:
    # This has to be "PassedWeight" as we want the distribution curves to be the same height
    norm_section="""double normfac = 1.0/PassedWeight;"""
else:
    # This has to be total weight as this is the true input weight to the curves
    norm_section="""double normfac = crossSection()/sumOfWeights();"""

f.write(norm_section)

section="""
      cout << TotWeight << "\t" << PassedWeight << "\t" << sumOfWeights() << "\t" << crossSection() << endl;
      cout << endl << "Cross-section used for normalisation: " << crossSection() << endl << endl;
      cout << "Normalisation factor used here:       " << normfac << endl;

      //leptons
      scale( _h_lepton_phi, normfac);
      scale( _h_lepton_eta, normfac);
      scale( _h_lepton_pT, normfac);
      
      scale( _h_neutrino_phi, normfac);
      scale( _h_neutrino_eta, normfac);
      scale( _h_neutrino_pT, normfac);
      
      scale( _h_met, normfac);
      scale( _h_met_phi, normfac);
      scale( _h_MT_dilep, normfac);
      scale( _h_mll, normfac);
      scale( _h_boson_pT, normfac);
      scale( _h_RATIOboson_pT, normfac);
      scale( _h_ALTboson_pT, normfac);
      scale( _h_njet_inclusive, normfac);
      //my_jets
      scale( _h_jet1_phi, normfac);
      scale( _h_jet1_eta, normfac);
      scale( _h_jet1_pT, normfac);
      scale( _h_jet2_phi, normfac);
      scale( _h_jet2_eta, normfac);
      scale( _h_jet2_pT, normfac);
      scale( _h_jet3_phi, normfac);
      scale( _h_jet3_eta, normfac);
      scale( _h_jet3_pT, normfac);
      scale( _h_mjj, normfac);
      scale( _h_ALTmjj, normfac);
      scale( _h_zoomed_mjj, normfac);
      scale( _h_RATIOmjj, normfac);
      scale( _h_scalarSumPt, normfac);
      scale( _h_vectorSumPt, normfac);
      scale( _h_deltaEta, normfac);
      scale( _h_delatPhi, normfac);
      scale( _h_deltaR, normfac);
      scale( _h_etaProduct, normfac);    
      scale( _h_dijetPt, normfac);  
      scale( _h_ALTdijet_pT, normfac);
      scale( _h_RATIOdijetPt, normfac); 
      scale( _h_gapref, normfac);           
      scale( _h_OppHemi, normfac);

      //other
      scale( _h_MT, normfac);
      //scale( _h_mllmjj, normfac);
      scale( _h_deltaPhiJetLep, normfac);
      scale( _h_totVisibleM, normfac);
      scale( _h_totVBFM, normfac);
      scale( _h_delabosonJetM, normfac);
      scale( _h_ClosestBosM, normfac);
      scale( _h_dijetEtaStar, normfac);
      scale( _h_WorZ, normfac);
      scale( _h_dR_lepJet, normfac);
      scale( _h_PID, normfac);
      scale( _h_MomFac, normfac);
      scale( _h_MomFac350, normfac);
      scale( _h_MomFac1000, normfac);
      scale( _h_MomFacX, normfac);
      scale( _h_MomFacY, normfac);
      scale( _h_PtBal, normfac);
      scale( _h_Weight, normfac);
      scale( _h_Zeppenfeld, normfac);
      scale( _h_BDPtSum, normfac);
      scale( _h_BDPtSum350, normfac);
      scale( _h_BDPtSum1000, normfac);
      scale( _h_BDPtSum2Jets, normfac);
      scale( _h_DPtInW, normfac);
      scale( _h_DPtInW350, normfac);
      scale( _h_DPtInW1000, normfac);

      scale( _h_VBFAngle, normfac );
      scale( _h_VBFAngle_cuts, normfac );
      scale( _h_VBFAngle_B, normfac );
      scale( _h_VBFAngle_cuts_B, normfac );
    }


  private:

    //leptons
    Histo1DPtr _h_lepton_phi;
    Histo1DPtr _h_lepton_eta;
    Histo1DPtr _h_lepton_pT;
    
    Histo1DPtr _h_neutrino_phi;
    Histo1DPtr _h_neutrino_eta;
    Histo1DPtr _h_neutrino_pT;
   
    Histo1DPtr _h_met;
    Histo1DPtr _h_met_phi;
    Histo1DPtr _h_MT_dilep;
    Histo1DPtr _h_mll;
    Histo1DPtr _h_boson_pT;
    Histo1DPtr _h_RATIOboson_pT;
    Histo1DPtr _h_ALTboson_pT;
    Histo1DPtr _h_njet_inclusive;
    //jets
    Histo1DPtr _h_jet1_phi;
    Histo1DPtr _h_jet1_eta;
    Histo1DPtr _h_jet1_pT;
    Histo1DPtr _h_jet2_phi;
    Histo1DPtr _h_jet2_eta;
    Histo1DPtr _h_jet2_pT;
    Histo1DPtr _h_jet3_phi;
    Histo1DPtr _h_jet3_eta;
    Histo1DPtr _h_jet3_pT;
    Histo1DPtr _h_mjj;
    Histo1DPtr _h_ALTmjj;
    Histo1DPtr _h_zoomed_mjj;
    Histo1DPtr _h_RATIOmjj;
    Histo1DPtr _h_scalarSumPt;
    Histo1DPtr _h_vectorSumPt;
    Histo1DPtr _h_deltaEta;
    Histo1DPtr _h_delatPhi;
    Histo1DPtr _h_deltaR;
    Histo1DPtr _h_etaProduct;   
    Histo1DPtr _h_dijetPt;
    Histo1DPtr _h_ALTdijet_pT;
    Histo1DPtr _h_RATIOdijetPt;
    Histo1DPtr _h_gapref;        
    Histo1DPtr _h_OppHemi;
    //other
    Histo1DPtr _h_MT;
    Histo1DPtr _h_deltaPhiJetLep;
    Histo1DPtr _h_totVisibleM;
    Histo1DPtr _h_totVBFM;
    Histo1DPtr _h_delabosonJetM;
    Histo1DPtr _h_ClosestBosM;
    Histo1DPtr _h_dijetEtaStar;
    Histo1DPtr _h_WorZ;
    Histo1DPtr _h_dR_lepJet;
    Histo1DPtr _h_PID;
    Histo1DPtr _h_MomFac;
    Histo1DPtr _h_MomFac350;
    Histo1DPtr _h_MomFac1000; 
    Histo1DPtr _h_MomFacX;
    Histo1DPtr _h_MomFacY;
    Histo1DPtr _h_PtBal;
    Histo1DPtr _h_Weight;
    Histo1DPtr _h_Zeppenfeld; 
    Histo1DPtr _h_BDPtSum;
    Histo1DPtr _h_BDPtSum350;
    Histo1DPtr _h_BDPtSum1000;
    Histo1DPtr _h_BDPtSum2Jets;
    Histo1DPtr _h_DPtInW;
    Histo1DPtr _h_DPtInW350;
    Histo1DPtr _h_DPtInW1000;

    Histo1DPtr _h_VBFAngle;
    Histo1DPtr _h_VBFAngle_cuts;
    Histo1DPtr _h_VBFAngle_B;
    Histo1DPtr _h_VBFAngle_cuts_B;
    
    
    };
  DECLARE_RIVET_PLUGIN(MC_VBF);
}
"""

f.write(section)
f.close()

# This line is used to compile the analysis just after it is written - the standard rivet-build was seen to not find a library so this custom script is used
os.system( "cd /Users/Alex/Desktop/HEP/Rivet_Analyses/MC_VBF; ./compile.sh")
