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


  class WJETS_SYST_ANALYSIS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    WJETS_SYST_ANALYSIS()
      : Analysis("WJETS_SYST_ANALYSIS")
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
      LeptonClusters leptons(fs, allleptons,
                     0.1, true,
                     etaRanges, 20.0*GeV);
					// LeptonClusters leptons(fs, allleptons,
				    //               0.1, true,
				    //               etaRanges, 15.0*GeV)
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
        _h_NjetIncl[i] = bookHistogram1D(1, 1, i+1);
        _h_RatioNjetIncl[i] = bookDataPointSet(2, 1, i+1);
        _h_FirstJetPt_1jet[i] = bookHistogram1D(3, 1, i+1);
        _h_FirstJetPt_2jet[i] = bookHistogram1D(4, 1, i+1);
        _h_FirstJetPt_3jet[i] = bookHistogram1D(5, 1, i+1);
        _h_FirstJetPt_4jet[i] = bookHistogram1D(6, 1, i+1);
        _h_SecondJetPt_2jet[i] = bookHistogram1D(7, 1, i+1);
        _h_SecondJetPt_3jet[i] = bookHistogram1D(8, 1, i+1);
        _h_SecondJetPt_4jet[i] = bookHistogram1D(9, 1, i+1);
        _h_ThirdJetPt_3jet[i] = bookHistogram1D(10, 1, i+1);
        _h_ThirdJetPt_4jet[i] = bookHistogram1D(11, 1, i+1);
        _h_FourthJetPt_4jet[i] = bookHistogram1D(12, 1, i+1);
        _h_Ht_1jet[i] = bookHistogram1D(13, 1, i+1);
        _h_Ht_2jet[i] = bookHistogram1D(14, 1, i+1);
        _h_Ht_3jet[i] = bookHistogram1D(15, 1, i+1);
        _h_Ht_4jet[i] = bookHistogram1D(16, 1, i+1);
        _h_Minv_2jet[i] = bookHistogram1D(17, 1, i+1);
        _h_Minv_3jet[i] = bookHistogram1D(18, 1, i+1);
        _h_Minv_4jet[i] = bookHistogram1D(19, 1, i+1);
        _h_JetRapidity[i] = bookHistogram1D(20, 1, i+1);
        _h_DeltaYElecJet[i] = bookHistogram1D(21, 1, i+1);
        _h_SumYElecJet[i] = bookHistogram1D(22, 1, i+1);
        _h_DeltaR_2jet[i] = bookHistogram1D(23, 1, i+1);
        _h_DeltaY_2jet[i] = bookHistogram1D(24, 1, i+1);
        _h_DeltaPhi_2jet[i] = bookHistogram1D(25, 1, i+1);
      }
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const vector<ClusteredLepton>& leptons = applyProjection<LeptonClusters>(event, "leptons").clusteredLeptons();
      ParticleVector neutrinos = applyProjection<FinalState>(event, "neutrinos").particlesByPt();

      if (leptons.size()!=1 || (neutrinos.size()==0)) {
        vetoEvent;
      }

      FourMomentum lepton=leptons[0].momentum();
      FourMomentum p_miss = neutrinos[0].momentum();
      if (p_miss.Et()<25.0*GeV) {
        vetoEvent;
      }
	  //if (p_miss.Et()<20.0*GeV) {
	  //	vetoEvent;
	  //}

      double mT=sqrt(2.0*lepton.pT()*p_miss.Et()*(1.0-cos(lepton.phi()-p_miss.phi())));
      if (mT<40.0*GeV) {
        vetoEvent;
      }
	  //Need to Define an if statement that vetos the event if the dijet mass is less/greater than 350*GeV

      double jetcuts[] = {30.0*GeV, 20.0*GeV};
	  //double jetcuts[] = {45.0*GeV, 35.0*GeV} //We want pT1 > 45.0*GeV and pT2 > 35.0*GeV
      const FastJets& jetpro = applyProjection<FastJets>(event, "jets");

      for (size_t i=0; i<2; ++i) {
        vector<FourMomentum> jets;
        double HT=lepton.pT()+p_miss.pT();
        foreach (const Jet& jet, jetpro.jetsByPt(jetcuts[i])) {
          if (fabs(jet.momentum().rapidity())<4.4 && //Remove deltaR cut and use deltaPhi cut instead
              deltaR(lepton, jet.momentum())>0.5) {
            jets.push_back(jet.momentum());
            HT += jet.momentum().pT();
          }
        }

        _h_NjetIncl[i]->fill(0.0, weight);

        // Njet>=1 observables
        if (jets.size()<1) continue;
        _h_NjetIncl[i]->fill(1.0, weight);
        _h_FirstJetPt_1jet[i]->fill(jets[0].pT(), weight);
        _h_JetRapidity[i]->fill(jets[0].rapidity(), weight);
        _h_Ht_1jet[i]->fill(HT, weight);
        _h_DeltaYElecJet[i]->fill(lepton.rapidity()-jets[0].rapidity(), weight);
        _h_SumYElecJet[i]->fill(lepton.rapidity()+jets[0].rapidity(), weight);

        // Njet>=2 observables
        if (jets.size()<2) continue;
        _h_NjetIncl[i]->fill(2.0, weight);
        _h_FirstJetPt_2jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_2jet[i]->fill(jets[1].pT(), weight);
        _h_Ht_2jet[i]->fill(HT, weight);
        double m2_2jet = FourMomentum(jets[0]+jets[1]).mass2();
        _h_Minv_2jet[i]->fill(m2_2jet>0.0 ? sqrt(m2_2jet) : 0.0, weight);
        _h_DeltaR_2jet[i]->fill(deltaR(jets[0], jets[1]), weight);
        _h_DeltaY_2jet[i]->fill(jets[0].rapidity()-jets[1].rapidity(), weight);
        _h_DeltaPhi_2jet[i]->fill(deltaPhi(jets[0], jets[1]), weight);

        // Njet>=3 observables
        if (jets.size()<3) continue;
        _h_NjetIncl[i]->fill(3.0, weight);
        _h_FirstJetPt_3jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_3jet[i]->fill(jets[1].pT(), weight);
        _h_ThirdJetPt_3jet[i]->fill(jets[2].pT(), weight);
        _h_Ht_3jet[i]->fill(HT, weight);
        double m2_3jet = FourMomentum(jets[0]+jets[1]+jets[2]).mass2();
        _h_Minv_3jet[i]->fill(m2_3jet>0.0 ? sqrt(m2_3jet) : 0.0, weight);

        // Njet>=4 observables
        if (jets.size()<4) continue;
        _h_NjetIncl[i]->fill(4.0, weight);
        _h_FirstJetPt_4jet[i]->fill(jets[0].pT(), weight);
        _h_SecondJetPt_4jet[i]->fill(jets[1].pT(), weight);
        _h_ThirdJetPt_4jet[i]->fill(jets[2].pT(), weight);
        _h_FourthJetPt_4jet[i]->fill(jets[3].pT(), weight);
        _h_Ht_4jet[i]->fill(HT, weight);
        double m2_4jet = FourMomentum(jets[0]+jets[1]+jets[2]+jets[3]).mass2();
        _h_Minv_4jet[i]->fill(m2_4jet>0.0 ? sqrt(m2_4jet) : 0.0, weight);

        // Njet>=5 observables
        if (jets.size()<5) continue;
        _h_NjetIncl[i]->fill(5.0, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i=0; i<2; ++i) {
        // first construct jet multi ratio
        int Nbins = _h_NjetIncl[i]->axis().bins();
        std::vector<double> ratio(Nbins-1, 0.0);
        std::vector<double> err(Nbins-1, 0.0);
        for (int n = 0; n < Nbins-1; ++n) {
          if (_h_NjetIncl[i]->binHeight(n) > 0.0 && _h_NjetIncl[i]->binHeight(n+1) > 0.0) {
            ratio[n] = _h_NjetIncl[i]->binHeight(n+1)/_h_NjetIncl[i]->binHeight(n);
            double relerr_n = _h_NjetIncl[i]->binError(n)/_h_NjetIncl[i]->binHeight(n);
            double relerr_m = _h_NjetIncl[i]->binError(n+1)/_h_NjetIncl[i]->binHeight(n+1);
            err[n] = ratio[n] * (relerr_n + relerr_m);
          }
        }
        _h_RatioNjetIncl[i]->setCoordinate(1, ratio, err);

        // scale all histos to the cross section
        double factor = crossSection()/sumOfWeights();
        scale(_h_DeltaPhi_2jet[i], factor);
        scale(_h_DeltaR_2jet[i], factor);
        scale(_h_DeltaY_2jet[i], factor);
        scale(_h_DeltaYElecJet[i], factor);
        scale(_h_FirstJetPt_1jet[i], factor);
        scale(_h_FirstJetPt_2jet[i], factor);
        scale(_h_FirstJetPt_3jet[i], factor);
        scale(_h_FirstJetPt_4jet[i], factor);
        scale(_h_FourthJetPt_4jet[i], factor);
        scale(_h_Ht_1jet[i], factor);
        scale(_h_Ht_2jet[i], factor);
        scale(_h_Ht_3jet[i], factor);
        scale(_h_Ht_4jet[i], factor);
        scale(_h_JetRapidity[i], factor);
        scale(_h_Minv_2jet[i], factor);
        scale(_h_Minv_3jet[i], factor);
        scale(_h_Minv_4jet[i], factor);
        scale(_h_NjetIncl[i], factor);
        scale(_h_SecondJetPt_2jet[i], factor);
        scale(_h_SecondJetPt_3jet[i], factor);
        scale(_h_SecondJetPt_4jet[i], factor);
        scale(_h_SumYElecJet[i], factor);
        scale(_h_ThirdJetPt_3jet[i], factor);
        scale(_h_ThirdJetPt_4jet[i], factor);
      }
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_h_DeltaPhi_2jet[2];
    AIDA::IHistogram1D *_h_DeltaR_2jet[2];
    AIDA::IHistogram1D *_h_DeltaY_2jet[2];
    AIDA::IHistogram1D *_h_DeltaYElecJet[2];
    AIDA::IHistogram1D *_h_FirstJetPt_1jet[2];
    AIDA::IHistogram1D *_h_FirstJetPt_2jet[2];
    AIDA::IHistogram1D *_h_FirstJetPt_3jet[2];
    AIDA::IHistogram1D *_h_FirstJetPt_4jet[2];
    AIDA::IHistogram1D *_h_FourthJetPt_4jet[2];
    AIDA::IHistogram1D *_h_Ht_1jet[2];
    AIDA::IHistogram1D *_h_Ht_2jet[2];
    AIDA::IHistogram1D *_h_Ht_3jet[2];
    AIDA::IHistogram1D *_h_Ht_4jet[2];
    AIDA::IHistogram1D *_h_JetRapidity[2];
    AIDA::IHistogram1D *_h_Minv_2jet[2];
    AIDA::IHistogram1D *_h_Minv_3jet[2];
    AIDA::IHistogram1D *_h_Minv_4jet[2];
    AIDA::IHistogram1D *_h_NjetIncl[2];
    AIDA::IDataPointSet *_h_RatioNjetIncl[2];
    AIDA::IHistogram1D *_h_SecondJetPt_2jet[2];
    AIDA::IHistogram1D *_h_SecondJetPt_3jet[2];
    AIDA::IHistogram1D *_h_SecondJetPt_4jet[2];
    AIDA::IHistogram1D *_h_SumYElecJet[2];
    AIDA::IHistogram1D *_h_ThirdJetPt_3jet[2];
    AIDA::IHistogram1D *_h_ThirdJetPt_4jet[2];
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WJETS_SYST_ANALYSIS);

}
