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
      IdentifiedFinalState allleptons;
      allleptons.acceptIdPair(ELECTRON);
      allleptons.acceptIdPair(MUON);
      std::vector<std::pair<double, double> > etaRanges;
      etaRanges.push_back(make_pair(-2.5, 2.5));
	  LeptonClusters leptons(fs, allleptons, 0.1, true, etaRanges, 15.0*GeV);		// pT_min of lepton is 15.0*GeV
      addProjection(leptons, "leptons");

      // Leading neutrinos for Etmiss
      LeadingParticlesFinalState neutrinos(fs);
      neutrinos.addParticleIdPair(NU_E);
      neutrinos.addParticleIdPair(NU_MU);
      neutrinos.setLeadingOnly(true);
      addProjection(neutrinos, "neutrinos");

      // Input for the jets: "Neutrinos, electrons, and muons from decays of the
      // massive W boson were not used"
      VetoedFinalState veto;
      veto.addVetoOnThisFinalState(leptons);
      veto.addVetoOnThisFinalState(neutrinos);
      FastJets jets(veto, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      addProjection(jets, "jets");

      for (size_t i=0; i<2; ++i) {
		// New Histograms defined here (must book them here and define them in .plots file.)
		string String = static_cast<ostringstream*>( &(ostringstream() << i+1) )->str();
		
        _h_NJetExcl[i] = bookHistogram1D("NJetExcl_"+String,4,0,4);
        _h_RatioNJetExcl[i] = bookDataPointSet("RatioNJetExcl_"+String,5,1,5);
        _h_FirstJetPt_2jet[i] = bookHistogram1D("FirstJetPt_2jet_"+String,50,0.0,500.0);
        _h_FirstJetPt_3jet[i] = bookHistogram1D("FirstJetPt_3jet_"+String,50,0.0,500.0);
        _h_FirstJetPt_4jet[i] = bookHistogram1D("FirstJetPt_4jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_2jet[i] = bookHistogram1D("SecondJetPt_2jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_3jet[i] = bookHistogram1D("SecondJetPt_3jet_"+String,50,0.0,500.0);
        _h_SecondJetPt_4jet[i] = bookHistogram1D("SecondJetPt_4jet_"+String,50,0.0,500.0);
        _h_ThirdJetPt_3jet[i] = bookHistogram1D("ThirdJetPt_3jet_"+String,50,0.0,500.0);
        _h_ThirdJetPt_4jet[i] = bookHistogram1D("ThirdJetPt_4jet_"+String,50,0.0,500.0);
        _h_FourthJetPt_4jet[i] = bookHistogram1D("FourthJetPt_4jet_"+String,50,0.0,500.0);
        _h_Ht_2jet[i] = bookHistogram1D("Ht_2jet_"+String,150,0.0,1500.0);
        _h_Ht_3jet[i] = bookHistogram1D("Ht_3jet_"+String,150,0.0,1500.0);
        _h_Ht_4jet[i] = bookHistogram1D("Ht_4jet_"+String,150,0.0,1500.0);
        _h_Minv_2jet[i] = bookHistogram1D("Minv_2jet_"+String,150,0.0,1500.0);
        _h_Minv_3jet[i] = bookHistogram1D("Minv_3jet_"+String,150,0.0,1500.0);
        _h_Minv_4jet[i] = bookHistogram1D("Minv_4jet_"+String,150,0.0,1500.0);
        _h_JetRapidity[i] = bookHistogram1D("JetRapidity_"+String,100,-5.0,5.0);
        _h_DeltaYElecJet[i] = bookHistogram1D("DeltaYElecJet_"+String,140,-7.0,7.0);
        _h_SumYElecJet[i] = bookHistogram1D("SumYElecJet_"+String,140,-7.0,7.0);
        _h_DeltaR_2jet[i] = bookHistogram1D("DeltaR_2jet_"+String,40,0.0,10.0);
        _h_DeltaY_2jet[i] = bookHistogram1D("DeltaY_2jet_"+String,140,-7.0,7.0);
        _h_DeltaPhi_2jet[i] = bookHistogram1D("DeltaPhi_2jet_"+String,35,0,3.5);
		_h_DijetMass_2jet[i] = bookHistogram1D("DijetMass_2jet_"+String,60,300.0,1500.0);
		_h_DijetMass_3jet[i] = bookHistogram1D("DijetMass_3jet_"+String,60,300.0,1500.0);
		_h_DijetMass_4jet[i] = bookHistogram1D("DijetMass_4jet_"+String,60,300.0,1500.0);
		_h_AntiDijetMass_2jet[i] = bookHistogram1D("AntiDijetMass_2jet_"+String,75,0.0,1500.0);
		_h_AntiDijetMass_3jet[i] = bookHistogram1D("AntiDijetMass_3jet_"+String,75,0.0,1500.0);
		_h_AntiDijetMass_4jet[i] = bookHistogram1D("AntiDijetMass_4jet_"+String,75,0.0,1500.0);
		_h_ThirdZep_3jet[i] = bookHistogram1D("ThirdZep_3jet_"+String,25,-5.0,5.0);
		_h_ThirdZep_4jet[i] = bookHistogram1D("ThirdZep_4jet_"+String,25,-5.0,5.0);
		_h_FourthZep_4jet[i] = bookHistogram1D("FourthZep_4jet_"+String,25,-5.0,5.0);
		_h_AntiDijetEtaDiff_2jet[i] = bookHistogram1D("AntiDijetEtaDiff_2jet_"+String,50,0.0,10.0);
		_h_AntiDijetEtaDiff_3jet[i] = bookHistogram1D("AntiDijetEtaDiff_3jet_"+String,50,0.0,10.0);
		_h_AntiDijetEtaDiff_4jet[i] = bookHistogram1D("AntiDijetEtaDiff_4jet_"+String,50,0.0,10.0);
		_h_AntiDijetPhiDiff_2jet[i] = bookHistogram1D("AntiDijetPhiDiff_2jet_"+String,20,-1.0,1.0);
		_h_AntiDijetPhiDiff_3jet[i] = bookHistogram1D("AntiDijetPhiDiff_3jet_"+String,20,-1.0,1.0);
		_h_AntiDijetPhiDiff_4jet[i] = bookHistogram1D("AntiDijetPhiDiff_4jet_"+String,20,-1.0,1.0);
		
	  _h_CutFlow = bookHistogram1D("CutFlow",10,0.0,10.0);
	  _h_Weights = bookHistogram1D("Weights",25,0.0,1.0);

      }
    }
	
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
	  _h_Weights->fill(weight);
	  //size_t event_number = numEvents();	//Could be used to skip broken events in MC containers

      const vector<ClusteredLepton>& leptons = applyProjection<LeptonClusters>(event, "leptons").clusteredLeptons();
      ParticleVector neutrinos = applyProjection<FinalState>(event, "neutrinos").particlesByPt();
	  _h_CutFlow->fill(1.0,weight);									// All Events (1.0) -> No Cuts
	
      if (leptons.size()!=1 || (neutrinos.size()==0)) {				// Keep events with exactly one lepton and at least one paired neutrino
        vetoEvent;
      }
	  _h_CutFlow->fill(2.0,weight);									// First Cut (2.0) -> Events with one lepton and at least one neutrino
	
      FourMomentum lepton = leptons[0].momentum();
      FourMomentum p_miss = neutrinos[0].momentum();
      if (p_miss.Et()<25.0*GeV) {									// MET must be greater than 25.0*GeV
        vetoEvent;
      }
	  _h_CutFlow->fill(3.0,weight);									// Second Cut (3.0) -> Events with MET greater than 25.0*GeV

      double mT=sqrt(2.0*lepton.pT()*p_miss.Et()*(1.0-cos(lepton.phi()-p_miss.phi())));
      if (mT<40.0*GeV) {											// mT(W) must be greater than 40.0*GeV
        vetoEvent;
      }
	  _h_CutFlow->fill(4.0,weight);									// Third Cut (4.0) -> Events with mT(W) greater than 40.0*GeV

	  double W_phi = FourMomentum(lepton + p_miss).phi();

	  const FastJets& jetpro = applyProjection<FastJets>(event, "jets");
	  vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(20.0*GeV)) {
		if (fabs(jet.momentum().rapidity())<4.4) {
		  jets.push_back(jet.momentum());
		}
      }
		
	  if (jets.size() < 2) vetoEvent;								// Keep dijet events only
	  _h_CutFlow->fill(5.0,weight);									// Fourth Cut (5.0) -> Events that are dijets
	  if (jets[0].pT() < 45*GeV) vetoEvent;							// pT_1 must be greater than 45*GeV
	  _h_CutFlow->fill(6.0,weight);									// Fifth Cut (6.0) -> Events with pT_1 greater than 45*GeV
	  if (jets[1].pT() < 35*GeV) vetoEvent;							// pT_2 must be greater than 35*GeV
	  _h_CutFlow->fill(7.0,weight);									// Sixth Cut (7.0) -> Events with pT_2 greater than 35*GeV
	  if (fabs(W_phi - jets[0].phi()) < 2.5) vetoEvent;				// DPhi(W,jet1) must be greater than 2.5
	  _h_CutFlow->fill(8.0,weight);									// Seventh Cut (8.0) -> Events with DeltaPhi(W,Jet_1) > 2.5

	  double dijet_mass2 = FourMomentum(jets[0]+jets[1]).mass2();
	  if (dijet_mass2 < 0.0) vetoEvent;								// Veto events with negative m_jj^2
	  double dijet_mass = sqrt(dijet_mass2);
	  if (dijet_mass < 350.0*GeV) vetoEvent;						// Veto event with m_jj < 350*GeV
	  _h_CutFlow->fill(9.0,weight);									// Eigth Cut (9.0) -> Events with Dijet Mass > 350.0*GeV

	  double jetcuts[] = {30.0*GeV, 25.0*GeV};						// All jets should be greater than 30*GeV or 25*GeV
      for (size_t i=0; i<2; ++i) {									
        vector<FourMomentum> jets;
		vector<FourMomentum> jetsByEta;
        double HT=lepton.pT()+p_miss.pT();
        foreach (const Jet& jet, jetpro.jetsByPt(jetcuts[i])) {
          if (fabs(jet.momentum().rapidity())<4.4) {
            jets.push_back(jet.momentum());
			jetsByEta.push_back(jet.momentum());
            HT += jet.momentum().pT();
          }
        }
		sort(jetsByEta.begin(),jetsByEta.end(),cmpMomByDescPseudorapidity);
		
		double antidijet_mass2 = FourMomentum(jetsByEta.front()+jetsByEta.back()).mass2();
		double antidijet_eta_diff = jetsByEta.front().eta() - jetsByEta.back().eta();
		double antidijet_phi_diff = cos(fabs(jetsByEta.front().phi()-jetsByEta.back().phi()));
		
        if (jets.size()<1) {
			cout << "Event had less than 2 jets...";
			continue;
		}
		
        // Njet>=2 observables
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

		_h_DijetMass_2jet[i]->fill(dijet_mass, weight);													// Dijet mass for Njet>=2, same as _h_Minv_2jet
		_h_AntiDijetMass_2jet[i]->fill(antidijet_mass2>0.0 ? sqrt(antidijet_mass2) : 0.0, weight);		// Anti-Dijet mass for Njet>=2
		_h_AntiDijetEtaDiff_2jet[i]->fill(antidijet_eta_diff, weight);
		_h_AntiDijetPhiDiff_2jet[i]->fill(antidijet_phi_diff, weight);

        // Njet>=3 observables
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

        // Njet>=4 observables
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
        if (jets.size()<5) continue;
        _h_NJetExcl[i]->fill(5.0, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i=0; i<2; ++i) {
        // first construct jet multi ratio
        int Nbins = _h_NJetExcl[i]->axis().bins();
        std::vector<double> ratio(Nbins-1, 0.0);
        std::vector<double> err(Nbins-1, 0.0);
        for (int n = 0; n < Nbins-1; ++n) {
          if (_h_NJetExcl[i]->binHeight(n) > 0.0 && _h_NJetExcl[i]->binHeight(n+1) > 0.0) {
            ratio[n] = _h_NJetExcl[i]->binHeight(n+1)/_h_NJetExcl[i]->binHeight(n);
            double relerr_n = _h_NJetExcl[i]->binError(n)/_h_NJetExcl[i]->binHeight(n);
            double relerr_m = _h_NJetExcl[i]->binError(n+1)/_h_NJetExcl[i]->binHeight(n+1);
            err[n] = ratio[n] * (relerr_n + relerr_m);
          }
        }
        _h_RatioNJetExcl[i]->setCoordinate(1, ratio, err);

        // scale all histos to the cross section
        double factor = crossSection()/sumOfWeights();
        scale(_h_DeltaPhi_2jet[i], factor);
        scale(_h_DeltaR_2jet[i], factor);
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
    AIDA::IDataPointSet *_h_RatioNJetExcl[2];
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
	AIDA::IHistogram1D *_h_CutFlow;
	AIDA::IHistogram1D *_h_Weights;
	
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WJETS_SYST_NEWANALYSIS);

}

// Just find max/min value in jetsEta, then keep the index for jets[]
// This finds the most forward and rearward jets in eta
// cout << "Event: " << event_number;
// cout << jet_string;
// cout << jetsEta;
// double eta_max;
// double eta_min;
// FourMomentum jet_eta_max;
// FourMomentum jet_eta_min;
// if (jetsEta[0] > jetsEta[1]) {
// 	eta_max = jetsEta[0];
// 	jet_eta_max = jets[0];
// 	eta_min = jetsEta[1];
// 	jet_eta_min = jets[1];
// }
// else{
// 	eta_max = jetsEta[1];
// 	jet_eta_max = jets[1];
// 	eta_min = jetsEta[0];
// 	jet_eta_min = jets[0];
// }
// 
// if (jetsEta.size() > 2) {
// 	for (size_t index=2; index<jetsEta.size(); index++) {
// 		if (eta_max < jetsEta[index]) {
// 			eta_max = jetsEta[index];
// 			jet_eta_max = jets[index];
// 		}
// 		else if (eta_min > jetsEta[index]) {
// 			eta_min = jetsEta[index];
// 			jet_eta_min = jets[index];
// 		}
// 	}
// }

// cout << "Event: " << event_number << endl;
// cout << "Jet 1 p_T: " << jets[0].pT() << "Jet 2 p_T: " << jets[1].pT() << "W Phi: " << W_phi << endl;
// cout << "Dijet Mass2: " << dijet_mass2 << ", Dijet Mass: " << dijet_mass << endl;
// cout << "Jets: " << jets << endl;
// cout << "Jets By Eta: " << jetsByEta << endl;
// cout << "Front Jet Eta: " << jetsByEta.front().eta() << ", Rear Jet Eta: " << jetsByEta.back().eta() << endl;
// cout << "FR Dijet Mass2: " << antidijet_mass2 << ", FR Eta Diff: " << antidijet_eta_diff << ", FR Phi Diff: " << antidijet_phi_diff << endl;
