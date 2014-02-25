// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"

double PassedWeight = 0;
double TotWeight = 0;

namespace Rivet {

  class WJETS_SYST_NEWANALYSIS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    WJETS_SYST_NEWANALYSIS()
      : Analysis("WJETS_SYST_NEWANALYSIS")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      WFinder wfinder(fs, -2.5, 2.5, 25.0*GeV, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfinder, "WFinder");
      FastJets jets( wfinder.remainingFinalState() , FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      addProjection(jets, "Jets_w");

      for (size_t i=0; i<2; ++i) {
        // New Histograms defined here (must book them here and define them in .plots file.)
        string String = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
		
        _h_NJetExcl[i] = bookHistogram1D("NJetExcl_"+String,4,0,4);
        _h_FirstJetPt_2jet[i] = bookHistogram1D("FirstJetPt_2jet_"+String,50,0.0,500.0);
        _h_FirstJetPt_3jet[i] = bookHistogram1D("FirstJetPt_3jet_"+String,50,0.0,500.0);
        _h_FirstJetPt_4jet[i] = bookHistogram1D("FirstJetPt_4jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_2jet[i] = bookHistogram1D("SecondJetPt_2jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_3jet[i] = bookHistogram1D("SecondJetPt_3jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_4jet[i] = bookHistogram1D("SecondJetPt_4jet_"+String,50,0.0,500.0);
        _h_ThirdJetPt_3jet[i] = bookHistogram1D("ThirdJetPt_3jet_"+String,50,0.0,500.0);
        _h_ThirdJetPt_4jet[i] = bookHistogram1D("ThirdJetPt_4jet_"+String,50,0.0,500.0);
        _h_FourthJetPt_4jet[i] = bookHistogram1D("FourthJetPt_4jet_"+String,50,0.0,500.0);
        _h_Ht_2jet[i] = bookHistogram1D("Ht_2jet_"+String,200,0.0,5000.0);
        _h_Ht_3jet[i] = bookHistogram1D("Ht_3jet_"+String,200,0.0,5000.0);
        _h_Ht_4jet[i] = bookHistogram1D("Ht_4jet_"+String,200,0.0,5000.0);
        _h_Minv_2jet[i] = bookHistogram1D("Minv_2jet_"+String,200,0.0,5000.0);
        _h_Minv_3jet[i] = bookHistogram1D("Minv_3jet_"+String,200,0.0,5000.0);
        _h_Minv_4jet[i] = bookHistogram1D("Minv_4jet_"+String,200,0.0,5000.0);
        _h_JetRapidity[i] = bookHistogram1D("JetRapidity_"+String,100,-5.0,5.0);
        _h_DeltaYElecJet[i] = bookHistogram1D("DeltaYElecJet_"+String,140,-7.0,7.0);
        _h_SumYElecJet[i] = bookHistogram1D("SumYElecJet_"+String,140,-7.0,7.0);
        _h_DeltaR_2jet[i] = bookHistogram1D("DeltaR_2jet_"+String,40,0.0,10.0);
        _h_DeltaR13_3jet[i] = bookHistogram1D("DeltaR13_3jet_"+String,40,0.0,10.0);
        _h_DeltaR23_3jet[i] = bookHistogram1D("DeltaR23_3jet_"+String,40,0.0,10.0);
        _h_DeltaY_2jet[i] = bookHistogram1D("DeltaY_2jet_"+String,140,-7.0,7.0);
        _h_DeltaPhi_2jet[i] = bookHistogram1D("DeltaPhi_2jet_"+String,35,0,3.5);
        _h_DijetMass_2jet[i] = bookHistogram1D("DijetMass_2jet_"+String,200,0.0,5000.0);
        _h_DijetMass_3jet[i] = bookHistogram1D("DijetMass_3jet_"+String,200,0.0,5000.0);
        _h_DijetMass_4jet[i] = bookHistogram1D("DijetMass_4jet_"+String,200,0.0,5000.0);
        _h_AntiDijetMass_2jet[i] = bookHistogram1D("AntiDijetMass_2jet_"+String,200,0.0,5000.0);
        _h_AntiDijetMass_3jet[i] = bookHistogram1D("AntiDijetMass_3jet_"+String,200,0.0,5000.0);
        _h_AntiDijetMass_4jet[i] = bookHistogram1D("AntiDijetMass_4jet_"+String,200,0.0,5000.0);
        _h_ThirdZep_3jet[i] = bookHistogram1D("ThirdZep_3jet_"+String,25,-5.0,5.0);
        _h_ThirdZep_4jet[i] = bookHistogram1D("ThirdZep_4jet_"+String,25,-5.0,5.0);
        _h_FourthZep_4jet[i] = bookHistogram1D("FourthZep_4jet_"+String,25,-5.0,5.0);
        _h_AntiDijetEtaDiff_2jet[i] = bookHistogram1D("AntiDijetEtaDiff_2jet_"+String,50,0.0,10.0);
        _h_AntiDijetEtaDiff_3jet[i] = bookHistogram1D("AntiDijetEtaDiff_3jet_"+String,50,0.0,10.0);
        _h_AntiDijetEtaDiff_4jet[i] = bookHistogram1D("AntiDijetEtaDiff_4jet_"+String,50,0.0,10.0);
        _h_AntiDijetPhiDiff_2jet[i] = bookHistogram1D("AntiDijetPhiDiff_2jet_"+String,20,-1.0,1.0);
        _h_AntiDijetPhiDiff_3jet[i] = bookHistogram1D("AntiDijetPhiDiff_3jet_"+String,20,-1.0,1.0);
        _h_AntiDijetPhiDiff_4jet[i] = bookHistogram1D("AntiDijetPhiDiff_4jet_"+String,20,-1.0,1.0);
        _h_CutFlow[i] = bookHistogram1D("CutFlow_"+String,12,0.0,12.0);
        _h_WeightCutFlow[i] = bookHistogram1D("WeightCutFlow_"+String,10,0.0,10.0);
		
        _h_Mjj_0ex[i] = bookHistogram1D("Mjj_Excl_00_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_1ex[i] = bookHistogram1D("Mjj_Excl_01_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_2ex[i] = bookHistogram1D("Mjj_Excl_02_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_3ex[i] = bookHistogram1D("Mjj_Excl_03_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_4ex[i] = bookHistogram1D("Mjj_Excl_04_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_5ex[i] = bookHistogram1D("Mjj_Excl_05_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_6ex[i] = bookHistogram1D("Mjj_Excl_06_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_7ex[i] = bookHistogram1D("Mjj_Excl_07_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_8ex[i] = bookHistogram1D("Mjj_Excl_08_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_9ex[i] = bookHistogram1D("Mjj_Excl_09_Jet_"+String,200,0.0,5000.0);
        _h_Mjj_10ex[i] = bookHistogram1D("Mjj_Excl_10_Jet_"+String,200,0.0,5000.0);

      }
	  
	  _h_NJetsNoCuts = bookHistogram1D("NJetsNoCuts",20.0,0.0,20.0);
    }
	
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      TotWeight+=weight;
      _h_CutFlow[0]->fill(1.0);
      _h_CutFlow[1]->fill(1.0);
      _h_WeightCutFlow[0]->fill(1.0,weight);
      _h_WeightCutFlow[1]->fill(1.0,weight);

      FourMomentum boson, lepton, neutrino, second_lepton;
      const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets_w");
      if ( wfinder.bosons().size() == 1 ){
            boson    = wfinder.bosons().front().momentum();
            lepton   = wfinder.constituentLeptons().front().momentum();
            neutrino = wfinder.constituentNeutrinos().front().momentum();
         }
      else 
            vetoEvent;

      double mT=sqrt(2.0*lepton.pT()*neutrino.Et()*(1.0-cos(lepton.phi()-neutrino.phi())));
      if (mT<40.0*GeV) {											// mT(W) must be greater than 40.0*GeV
        vetoEvent;
      }
      _h_CutFlow[0]->fill(2.0);
	  _h_CutFlow[1]->fill(2.0);
      _h_WeightCutFlow[0]->fill(2.0,weight);
      _h_WeightCutFlow[1]->fill(2.0,weight);

      // Jet Projection (only cares about jets with pT > 20 GeV)
      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(20.0*GeV)) {
          if ( fabs(jet.momentum().rapidity()) > 4.4 ) continue;
          if ( fabs(deltaR(jet, lepton)) < 0.3 ) continue;
          jets.push_back(jet.momentum());
      }
      _h_NJetsNoCuts->fill(jets.size());
	
      // Jet Selection	
      if (jets.size() < 2) vetoEvent;								// Keep dijet events only
      _h_CutFlow[0]->fill(3.0);
	  _h_CutFlow[1]->fill(3.0);
      _h_WeightCutFlow[0]->fill(3.0,weight);
      _h_WeightCutFlow[1]->fill(3.0,weight);
      if (jets[0].pT() < 80*GeV) vetoEvent;							// pT_1 must be greater than 80*GeV
      _h_CutFlow[0]->fill(4.0);
	  _h_CutFlow[1]->fill(4.0);
      _h_WeightCutFlow[0]->fill(4.0,weight);
      _h_WeightCutFlow[1]->fill(4.0,weight);
      if (jets[1].pT() < 60*GeV) vetoEvent;							// pT_2 must be greater than 60*GeV
      _h_CutFlow[0]->fill(5.0);
	  _h_CutFlow[1]->fill(5.0);
      _h_WeightCutFlow[0]->fill(5.0,weight);
      _h_WeightCutFlow[1]->fill(5.0,weight);

      double dijet_mass2 = FourMomentum(jets[0]+jets[1]).mass2();
      if (dijet_mass2 < 0.0) vetoEvent;								// Veto events with negative m_jj^2
      double dijet_mass = sqrt(dijet_mass2);
      if (dijet_mass < 500.0*GeV) vetoEvent;						// Veto event with m_jj < 500*GeV
      _h_CutFlow[0]->fill(6.0);
	  _h_CutFlow[1]->fill(6.0);
      _h_WeightCutFlow[0]->fill(6.0,weight);
      _h_WeightCutFlow[1]->fill(6.0,weight);

      double jetcuts[] = {30.0*GeV, 20.0*GeV};						// All jets should be greater than 30*GeV (20*GeV)
      for (size_t i=0; i<2; ++i) {
        vector<FourMomentum> jets;
        vector<FourMomentum> jetsByEta;
        double HT=lepton.pT()+neutrino.pT();
        foreach (const Jet& jet, jetpro.jetsByPt(jetcuts[i])) {
          if (fabs(jet.momentum().rapidity())<4.4) {
            jets.push_back(jet.momentum());
            jetsByEta.push_back(jet.momentum());
            HT += jet.momentum().pT();
          }
        }
		
        // Central Jet Veto (Vetos any event that has a jet between Etas of the two leading jets
        size_t nojets = jets.size();
        for (size_t j=2; j<nojets; ++j) {
          if ( nojets>=3 && jets[j].pT()>jetcuts[i] &&
             (
             ( jets[j].rapidity()<jets[0].rapidity() && jets[j].rapidity()>jets[1].rapidity() ) ||
             ( jets[j].rapidity()>jets[0].rapidity() && jets[j].rapidity()<jets[1].rapidity() )
             )
             ) {
                  vetoEvent;
          }
        }
        _h_CutFlow[i]->fill(7.0);
        _h_WeightCutFlow[i]->fill(7.0,weight);
	
        // Outside Lepton Veto (Vetos any event with leptons outside of the Etas of the two leading jets)
        bool outsideLeptonVeto = true;
  	    double lepton_eta = lepton.rapidity();
        if ( (lepton_eta>jets[0].rapidity() && lepton_eta<jets[1].rapidity()) ||
  	       (lepton_eta<jets[0].rapidity() && lepton_eta>jets[1].rapidity()) ) {
               outsideLeptonVeto = true; // in the rapidity interval defined by lead two jets
        } 
  	    else {
               outsideLeptonVeto = false;
  			 vetoEvent;
        }
        _h_CutFlow[i]->fill(8.0);
        _h_WeightCutFlow[i]->fill(8.0,weight);
          
        sort(jetsByEta.begin(),jetsByEta.end(),cmpMomByDescPseudorapidity);		// Sorts Jet list by Pseudorapidity
        double antidijet_mass2 = FourMomentum(jetsByEta.front()+jetsByEta.back()).mass2();
        double antidijet_eta_diff = jetsByEta.front().eta() - jetsByEta.back().eta();
        double antidijet_phi_diff = cos(fabs(jetsByEta.front().phi()-jetsByEta.back().phi()));
		
		if (jets.size()==0) _h_Mjj_0ex[i]->fill(dijet_mass, weight);
		
		if (jets.size()==1) _h_Mjj_1ex[i]->fill(dijet_mass, weight);
        if (jets.size()<1) {
          cout << "Event had less than 2 jets...";
          continue;
        }
		
        // Njet>=2 observables
		if (jets.size()==2) _h_Mjj_2ex[i]->fill(dijet_mass, weight);
        if (jets.size()<2) continue;
        _h_NJetExcl[i]->fill(2.0, weight);
        _h_FirstJetPt_2jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_2jet[i]->fill(jets[1].pT(), weight);
        _h_Ht_2jet[i]->fill(HT, weight);
        _h_JetRapidity[i]->fill(jets[0].rapidity(), weight);
        _h_DeltaYElecJet[i]->fill(lepton.rapidity()-jets[0].rapidity(), weight);
        _h_SumYElecJet[i]->fill(lepton.rapidity()+jets[0].rapidity(), weight);
        double m2_2jet = FourMomentum(jets[0]+jets[1]).mass2();
        _h_Minv_2jet[i]->fill(m2_2jet>0.0 ? sqrt(m2_2jet) : 0.0, weight);
        _h_DeltaR_2jet[i]->fill(deltaR(jets[0], jets[1]), weight);
        _h_DeltaY_2jet[i]->fill(jets[0].rapidity()-jets[1].rapidity(), weight);
        _h_DeltaPhi_2jet[i]->fill(deltaPhi(jets[0], jets[1]), weight);

        _h_DijetMass_2jet[i]->fill(dijet_mass, weight);								// Dijet mass for Njet>=2, same as _h_Minv_2jet
        _h_AntiDijetMass_2jet[i]->fill(antidijet_mass2>0.0 ? sqrt(antidijet_mass2) : 0.0, weight);		// Anti-Dijet mass for Njet>=2
        _h_AntiDijetEtaDiff_2jet[i]->fill(antidijet_eta_diff, weight);
        _h_AntiDijetPhiDiff_2jet[i]->fill(antidijet_phi_diff, weight);

        // Njet>=3 observables
		if (jets.size()==3) _h_Mjj_3ex[i]->fill(dijet_mass, weight);
        if (jets.size()<3) continue;
        _h_NJetExcl[i]->fill(3.0, weight);
        _h_FirstJetPt_3jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_3jet[i]->fill(jets[1].pT(), weight);
        _h_ThirdJetPt_3jet[i]->fill(jets[2].pT(), weight);
        _h_Ht_3jet[i]->fill(HT, weight);
        double m2_3jet = FourMomentum(jets[0]+jets[1]+jets[2]).mass2();
        _h_Minv_3jet[i]->fill(m2_3jet>0.0 ? sqrt(m2_3jet) : 0.0, weight);

        _h_DijetMass_3jet[i]->fill(dijet_mass, weight);													// Dijet mass for Njet>=3
        _h_AntiDijetMass_3jet[i]->fill(antidijet_mass2>0.0 ? sqrt(antidijet_mass2) : 0.0, weight);		// Anti-Dijet mass for Njet>=3
        _h_AntiDijetEtaDiff_3jet[i]->fill(antidijet_eta_diff, weight);
        _h_AntiDijetPhiDiff_3jet[i]->fill(antidijet_phi_diff, weight);
        double zep3 = jets[2].eta()-0.5*(jets[0].eta()+jets[1].eta());
        _h_ThirdZep_3jet[i]->fill(zep3, weight);														// Third Zeppenfeld for Njet>=3
        _h_DeltaR13_3jet[i]->fill(deltaR(jets[0], jets[2]), weight);
        _h_DeltaR23_3jet[i]->fill(deltaR(jets[1], jets[2]), weight);

        // Njet>=4 observables
		if (jets.size()==4) _h_Mjj_4ex[i]->fill(dijet_mass, weight);
        if (jets.size()<4) continue;
        _h_NJetExcl[i]->fill(4.0, weight);
        _h_FirstJetPt_4jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_4jet[i]->fill(jets[1].pT(), weight);
        _h_ThirdJetPt_4jet[i]->fill(jets[2].pT(), weight);
        _h_FourthJetPt_4jet[i]->fill(jets[3].pT(), weight);
        _h_Ht_4jet[i]->fill(HT, weight);
        double m2_4jet = FourMomentum(jets[0]+jets[1]+jets[2]+jets[3]).mass2();
        _h_Minv_4jet[i]->fill(m2_4jet>0.0 ? sqrt(m2_4jet) : 0.0, weight);

        _h_DijetMass_4jet[i]->fill(dijet_mass, weight);													// Dijet mass for Njet>=4
        _h_AntiDijetMass_4jet[i]->fill(antidijet_mass2>0.0 ? sqrt(antidijet_mass2) : 0.0, weight);		// Anti-Dijet mass for Njet>=4
        _h_AntiDijetEtaDiff_4jet[i]->fill(antidijet_eta_diff, weight);
        _h_AntiDijetPhiDiff_4jet[i]->fill(antidijet_phi_diff, weight);
        double zep4 = jets[3].eta()-0.5*(jets[0].eta()+jets[1].eta());
        _h_ThirdZep_4jet[i]->fill(zep3, weight);														// Third Zeppenfeld for Njet>=4
        _h_FourthZep_4jet[i]->fill(zep4, weight);														// Fourth Zeppenfeld for Njet>=4

        // Njet>=5 observables
		if (jets.size()==5) _h_Mjj_5ex[i]->fill(dijet_mass, weight);
        if (jets.size()<5) continue;
        _h_NJetExcl[i]->fill(5.0, weight);
		
		if (jets.size()==6) _h_Mjj_6ex[i]->fill(dijet_mass, weight);
		if (jets.size()==7) _h_Mjj_7ex[i]->fill(dijet_mass, weight);
		if (jets.size()==8) _h_Mjj_8ex[i]->fill(dijet_mass, weight);
		if (jets.size()==9) _h_Mjj_9ex[i]->fill(dijet_mass, weight);
		if (jets.size()==10) _h_Mjj_10ex[i]->fill(dijet_mass, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i=0; i<2; ++i) {

        // scale all histos to the cross section
		cout << crossSection() << endl << sumOfWeights();
        double factor = crossSection()/sumOfWeights();
        cout << TotWeight << "	" << PassedWeight << "	" << sumOfWeights() << "	" << crossSection() << endl;
        cout << endl << "Cross-section used for normalisation: " << crossSection() << endl << endl;
        cout << "Normalisation factor used here:       " << factor << endl;
		
        scale(_h_DeltaPhi_2jet[i], factor);
        scale(_h_DeltaR_2jet[i], factor);
		scale(_h_DeltaR13_3jet[i], factor);
		scale(_h_DeltaR23_3jet[i], factor);
        scale(_h_DeltaY_2jet[i], factor);
        scale(_h_DeltaYElecJet[i], factor);
        scale(_h_FirstJetPt_2jet[i], factor);
        scale(_h_FirstJetPt_3jet[i], factor);
        scale(_h_FirstJetPt_4jet[i], factor);
        scale(_h_FourthJetPt_4jet[i], factor);
        scale(_h_Ht_2jet[i], factor);
        scale(_h_Ht_3jet[i], factor);
        scale(_h_Ht_4jet[i], factor);
        scale(_h_JetRapidity[i], factor);
        scale(_h_Minv_2jet[i], factor);
        scale(_h_Minv_3jet[i], factor);
        scale(_h_Minv_4jet[i], factor);
        scale(_h_NJetExcl[i], factor);
        scale(_h_SecondJetPt_2jet[i], factor);
        scale(_h_SecondJetPt_3jet[i], factor);
        scale(_h_SecondJetPt_4jet[i], factor);
        scale(_h_SumYElecJet[i], factor);
        scale(_h_ThirdJetPt_3jet[i], factor);
        scale(_h_ThirdJetPt_4jet[i], factor);
        scale(_h_DijetMass_2jet[i], factor);
        scale(_h_DijetMass_3jet[i], factor);
        scale(_h_DijetMass_4jet[i], factor);
        scale(_h_AntiDijetMass_2jet[i], factor);
        scale(_h_AntiDijetMass_3jet[i], factor);
        scale(_h_AntiDijetMass_4jet[i], factor);
        scale(_h_ThirdZep_3jet[i], factor);
        scale(_h_ThirdZep_4jet[i], factor);
        scale(_h_FourthZep_4jet[i], factor);
        scale(_h_AntiDijetEtaDiff_2jet[i], factor);
        scale(_h_AntiDijetEtaDiff_3jet[i], factor);
        scale(_h_AntiDijetEtaDiff_4jet[i], factor);
        scale(_h_AntiDijetPhiDiff_2jet[i], factor);
        scale(_h_AntiDijetPhiDiff_3jet[i], factor);
        scale(_h_AntiDijetPhiDiff_4jet[i], factor);
		scale(_h_Mjj_0ex[i], factor);
		scale(_h_Mjj_1ex[i], factor);
		scale(_h_Mjj_2ex[i], factor);
		scale(_h_Mjj_3ex[i], factor);
		scale(_h_Mjj_4ex[i], factor);
		scale(_h_Mjj_5ex[i], factor);
		scale(_h_Mjj_6ex[i], factor);
		scale(_h_Mjj_7ex[i], factor);
		scale(_h_Mjj_8ex[i], factor);
		scale(_h_Mjj_9ex[i], factor);
		scale(_h_Mjj_10ex[i], factor);
      }
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
	AIDA::IHistogram1D *_h_DeltaPhi_2jet[2];
	AIDA::IHistogram1D *_h_DeltaR_2jet[2];
	AIDA::IHistogram1D *_h_DeltaR13_3jet[2];
	AIDA::IHistogram1D *_h_DeltaR23_3jet[2];
	AIDA::IHistogram1D *_h_DeltaY_2jet[2];
	AIDA::IHistogram1D *_h_DeltaYElecJet[2];
	AIDA::IHistogram1D *_h_FirstJetPt_2jet[2];
	AIDA::IHistogram1D *_h_FirstJetPt_3jet[2];
	AIDA::IHistogram1D *_h_FirstJetPt_4jet[2];
	AIDA::IHistogram1D *_h_FourthJetPt_4jet[2];
	AIDA::IHistogram1D *_h_Ht_2jet[2];
	AIDA::IHistogram1D *_h_Ht_3jet[2];
	AIDA::IHistogram1D *_h_Ht_4jet[2];
	AIDA::IHistogram1D *_h_JetRapidity[2];
	AIDA::IHistogram1D *_h_Minv_2jet[2];
	AIDA::IHistogram1D *_h_Minv_3jet[2];
	AIDA::IHistogram1D *_h_Minv_4jet[2];
	AIDA::IHistogram1D *_h_NJetExcl[2];
	AIDA::IHistogram1D *_h_SecondJetPt_2jet[2];
	AIDA::IHistogram1D *_h_SecondJetPt_3jet[2];
	AIDA::IHistogram1D *_h_SecondJetPt_4jet[2];
	AIDA::IHistogram1D *_h_SumYElecJet[2];
	AIDA::IHistogram1D *_h_ThirdJetPt_3jet[2];
	AIDA::IHistogram1D *_h_ThirdJetPt_4jet[2];
	AIDA::IHistogram1D *_h_DijetMass_2jet[2];
	AIDA::IHistogram1D *_h_DijetMass_3jet[2];
	AIDA::IHistogram1D *_h_DijetMass_4jet[2];
	AIDA::IHistogram1D *_h_AntiDijetMass_2jet[2];
	AIDA::IHistogram1D *_h_AntiDijetMass_3jet[2];
	AIDA::IHistogram1D *_h_AntiDijetMass_4jet[2];
	AIDA::IHistogram1D *_h_ThirdZep_3jet[2];
	AIDA::IHistogram1D *_h_ThirdZep_4jet[2];
	AIDA::IHistogram1D *_h_FourthZep_4jet[2];
	AIDA::IHistogram1D *_h_AntiDijetEtaDiff_2jet[2];
	AIDA::IHistogram1D *_h_AntiDijetEtaDiff_3jet[2];
	AIDA::IHistogram1D *_h_AntiDijetEtaDiff_4jet[2];
	AIDA::IHistogram1D *_h_AntiDijetPhiDiff_2jet[2];
	AIDA::IHistogram1D *_h_AntiDijetPhiDiff_3jet[2];
	AIDA::IHistogram1D *_h_AntiDijetPhiDiff_4jet[2];
	AIDA::IHistogram1D *_h_CutFlow[2];
    AIDA::IHistogram1D *_h_WeightCutFlow[2];
	AIDA::IHistogram1D *_h_NJetsNoCuts;
	AIDA::IHistogram1D *_h_Mjj_0ex[2];
	AIDA::IHistogram1D *_h_Mjj_1ex[2];
	AIDA::IHistogram1D *_h_Mjj_2ex[2];
	AIDA::IHistogram1D *_h_Mjj_3ex[2];
	AIDA::IHistogram1D *_h_Mjj_4ex[2];
	AIDA::IHistogram1D *_h_Mjj_5ex[2];
	AIDA::IHistogram1D *_h_Mjj_6ex[2];
	AIDA::IHistogram1D *_h_Mjj_7ex[2];
	AIDA::IHistogram1D *_h_Mjj_8ex[2];
	AIDA::IHistogram1D *_h_Mjj_9ex[2];
	AIDA::IHistogram1D *_h_Mjj_10ex[2];
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WJETS_SYST_NEWANALYSIS);

}
