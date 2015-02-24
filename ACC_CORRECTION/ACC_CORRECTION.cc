// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/RivetAIDA.hh"
//#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {

  using namespace Cuts;

  class ACC_CORRECTION : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ACC_CORRECTION()
      : Analysis("ACC_CORRECTION")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      WFinder wfinder(fs, etaIn(-2.5, 2.5) & (pT >= 25.0*GeV), PID::MUON, 60.0*GeV, 100.0*GeV, 0*GeV, 0.2);
      addProjection(wfinder, "WFinder");
      
      FastJets jets( wfinder.remainingFinalState() , FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      addProjection(jets, "Jets_w");
	  
      MissingMomentum missmom(FinalState(-5.0,5.0,0*GeV));
      addProjection(missmom, "mm");

      // BOOK HISTOGRAMS HERE!!!!!
      _h_CutFlow = bookHisto1D("CutFlow",12,0.0,12.0);
      _h_WeightCutFlow = bookHisto1D("WeightCutFlow",12,0.0,12.0);
	  
    }
	
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      _h_CutFlow->fill(1.0);
      _h_WeightCutFlow->fill(1.0,weight);

      FourMomentum boson, lepton, neutrino, second_lepton;
      const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets_w");
	  const MissingMomentum& missmom = applyProjection<MissingMomentum>(event, "mm");
      if ( wfinder.bosons().size() == 1 ){
            boson    = wfinder.bosons().front().momentum();
            lepton   = wfinder.constituentLeptons().front().momentum();
            neutrino = wfinder.constituentNeutrinos().front().momentum();
         }
      else vetoEvent;
	  
      //if ( neutrino.pT() < 25.0*GeV ) vetoEvent;
      if ( missmom.vectorEt().mod() < 25.0*GeV ) vetoEvent;
      //if ( lepton.pT() < 25.0*GeV ) vetoEvent;
      //if ( fabs(lepton.rapidity()) > 2.5 ) vetoEvent;
      _h_CutFlow->fill(2.0);
      _h_WeightCutFlow->fill(2.0,weight);

      // Jet Projection (only cares about jets with pT > 30 GeV)
      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(30.0*GeV)) {
          if ( fabs(jet.momentum().rapidity()) > 4.4 ) continue;
          //if ( fabs(deltaR(jet, lepton)) < 0.3 ) continue;
          jets.push_back(jet.momentum());
      }
	  
      
      double mT=sqrt(2.0*lepton.pT()*neutrino.Et()*(1.0-cos(lepton.phi()-neutrino.phi())));
      if (mT<40.0*GeV) vetoEvent;									// mT(W) must be greater than 40.0*GeV
      _h_CutFlow->fill(3.0);
      _h_WeightCutFlow->fill(3.0,weight);
	  
      if (jets.size() < 2) vetoEvent;								// Keep dijet events only
      _h_CutFlow->fill(4.0);
      _h_WeightCutFlow->fill(4.0,weight);
	  
      if (jets[0].pT() < 80*GeV) vetoEvent;							// pT_1 must be greater than 80*GeV
      _h_CutFlow->fill(5.0);
      _h_WeightCutFlow->fill(5.0,weight);
      
      if (jets[1].pT() < 60*GeV) vetoEvent;							// pT_2 must be greater than 60*GeV
      _h_CutFlow->fill(6.0);
      _h_WeightCutFlow->fill(6.0,weight);
	  
      double dijet_mass2 = FourMomentum(jets[0]+jets[1]).mass2();
      if (dijet_mass2 < 0.0) vetoEvent;								// Veto events with negative m_jj^2
      double dijet_mass = sqrt(dijet_mass2);

      if (dijet_mass < 500.0*GeV) vetoEvent;						// Veto event with m_jj < 500*GeV
      _h_CutFlow->fill(7.0);
      _h_WeightCutFlow->fill(7.0,weight);
		
      // Third Jet Centrality Cut (Vetos any event that has a jet between Etas of the two leading jets
      size_t nojets = jets.size();
      bool passJet3Centrality = false;
      double jet_1_eta = jets[0].rapidity();
      double jet_2_eta = jets[1].rapidity();
      if ( nojets < 3) passJet3Centrality = true;
      if ( nojets > 2 ) {
          double jet_3_eta = jets[2].rapidity();
          double jet3gap = fabs(( jet_3_eta - ((jet_1_eta + jet_2_eta)/2.0))/(jet_1_eta - jet_2_eta));
          if ( jet3gap > 0.4 ) passJet3Centrality = true;
      }
    
      // Lepton Centrality Cut (Vetos any event with leptons outside of the Etas of the two leading jets)
      double lepton_eta = lepton.rapidity();
      double lep_cent = fabs((lepton_eta - ((jet_1_eta + jet_2_eta)/2.0) )/(jet_1_eta - jet_2_eta));
      if ( lep_cent > 0.4 ) {
          vetoEvent;
      }
      _h_CutFlow->fill(8.0);
      _h_WeightCutFlow->fill(8.0,weight);
    
      if ( passJet3Centrality ) {
          _h_CutFlow->fill(9.0);
          _h_WeightCutFlow->fill(9.0,weight);
      }
      else {
          vetoEvent;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //double factor = crossSection()/sumOfWeights();
      
      //scale(_h_WeightCutFlow, factor);
      
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_CutFlow;
    Histo1DPtr _h_WeightCutFlow;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ACC_CORRECTION);

}
