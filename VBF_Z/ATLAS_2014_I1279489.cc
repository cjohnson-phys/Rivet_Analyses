// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Particle.hh"

#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {

  struct Plots {
    std::string label;

    Histo1DPtr h_dy;
    Histo1DPtr h_mjj;
    Histo1DPtr h_njets;
    Histo1DPtr h_dphijj;
    Histo1DPtr h_ptbal;

    Histo1DPtr h_jetveto_mjj_veto;
    Histo1DPtr h_jetveto_mjj_inc;
    Histo1DPtr h_jetveto_dy_veto;
    Histo1DPtr h_jetveto_dy_inc;

    Histo1DPtr h_ptbaleff_mjj_veto;
    Histo1DPtr h_ptbaleff_mjj_inc;
    Histo1DPtr h_ptbaleff_dy_veto;
    Histo1DPtr h_ptbaleff_dy_inc;

    Profile1DPtr p_avgnjets_dy;
    Profile1DPtr p_avgnjets_mjj;
  };

  class Variables {
    public:
      double m_jet1pt;
      double m_jet2pt;
      double m_zpt;

      double m_deltay;
      double m_mjj;
      double m_deltaphijj;
      double m_ptbalance2;
      double m_ptbalance3;
      int m_ngapjets;

      double m_dilepton_dr;

      bool m_pass_jetveto;
      bool m_pass_ptbaleff;

      bool IsBetween(const Jet* probe, const Jet* boundary1, const Jet* boundary2)
      {
        double y_p = probe->momentum().rapidity();
        double y_b1 = boundary1->momentum().rapidity();
        double y_b2 = boundary2->momentum().rapidity();

        double y_min = std::min(y_b1, y_b2);
        double y_max = std::max(y_b1, y_b2);

        if(y_p > y_min && y_p < y_max){ return true; }
        else return false;
      }

      int GetGapJets(vector<const Jet*>& jets, FourMomentum& thirdJet)
      {
        if(jets.size() < 2) return 0;
        // The vector of jets is already sorted by pT. So the boundary jets will be the first two.
        const Jet* bj1 = jets.at(0);
        const Jet* bj2 = jets.at(1);

        int n_between = 0;
        // Start loop at the 3rd hardest pT jet
        for(unsigned int i=2; i<jets.size(); ++i)
        {
          const Jet* j = jets.at(i);

          // If this jet is between the boundary jets and is hard enough, increment counter
          if(IsBetween(j, bj1, bj2)) {
            if (n_between == 0) thirdJet = j->momentum();
            ++n_between;
          }
        }
        return n_between;
      }

      Variables(vector<const Jet*> jets, FourMomentum* lep1, FourMomentum* lep2) {
        FourMomentum j1 = jets.at(0)->momentum();
        FourMomentum j2 = jets.at(1)->momentum();
        m_jet1pt = j1.pT();
        m_jet2pt = j2.pT();
        assert(m_jet1pt > m_jet2pt);

        m_zpt = (*lep1 + *lep2).pT();

        m_deltay = fabs(j1.rapidity() - j2.rapidity());
        m_mjj = (j1 + j2).mass();
        m_deltaphijj = deltaPhi(j1, j2) / PI;

        FourMomentum gapjet(0., 0., 0., 0.);
        m_ngapjets = GetGapJets(jets, gapjet);

        double ptbal_vec = (j1 + j2 + *lep1 + *lep2).pT();
        double ptbal_sc = j1.pT() + j2.pT() + lep1->pT() + lep2->pT();
        m_ptbalance2 = ptbal_vec / ptbal_sc;

        double ptbal3_vec = (j1 + j2 + gapjet + *lep1 + *lep2).pT();
        double ptbal3_sc = j1.pT() + j2.pT() + gapjet.pT() + lep1->pT() + lep2->pT();
        m_ptbalance3 = ptbal3_vec / ptbal3_sc;

        m_pass_jetveto = gapjet.pT() < 25.0*GeV;
        m_pass_ptbaleff = m_ptbalance2 < 0.15;
      }
  };


  class ATLAS_2014_I1279489 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1279489()
      : Analysis("ATLAS_2014_I1279489")
    {    }


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(-5.0, 5.0);
      addProjection(fs, "FS");

      IdentifiedFinalState photon_fs(fs);
      photon_fs.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState electron_fs(fs);
      electron_fs.acceptIdPair(PID::ELECTRON);

      IdentifiedFinalState muon_fs(fs);
      muon_fs.acceptIdPair(PID::MUON);

      std::vector<std::pair<double, double> > etaRanges;
      etaRanges.push_back(make_pair(-2.47, 2.47));

      DressedLeptons dressed_electrons(photon_fs, electron_fs, 0.1, true, etaRanges, 25.0*GeV, false);
      addProjection(dressed_electrons, "DRESSED_ELECTRONS");

      DressedLeptons dressed_muons(photon_fs, muon_fs, 0.1, true, etaRanges, 25.0*GeV, false);
      addProjection(dressed_muons, "DRESSED_MUONS");

      FastJets jets(fs, FastJets::ANTIKT, 0.4);
      addProjection(jets, "JETS");

      InitialisePlots(baseline_plots, "baseline");
      InitialisePlots(highpt_plots, "highpt");
      InitialisePlots(search_plots, "search");
      InitialisePlots(control_plots, "control");
      InitialisePlots(highmass_plots, "highmass");

    }

    void InitialisePlots(Plots& plots, std::string phase_space){
      /****************************************
       * Plot labeling:                       *
       * format = d0_-x0_-y0_                 *
       * d01 = baseline fiducial region       *
       * d02 = high-pt fiducial region        *
       * d03 = search fiducial region         *
       * d04 = control fiducial region        *
       * d05 = high-mass fiducial region      *
       *                                      *
       * x01 = mjj on x-axis                  *
       * x02 = delta-y on x-axis              *
       * x03 = njets on x-axis                *
       * x04 = dphijj on x-axis               *
       * x05 = ptbalance on x-axis            *
       *                                      *
       * y01 = differential cross-section     *
       * y02 = jet veto efficiency            *
       * y03 = ptbalance efficiency           *
       * y04 = average njets                  *
       ****************************************/
      plots.label = phase_space;

      if(phase_space=="baseline") {
        plots.h_mjj = bookHisto1D(1, 1, 1);
        plots.h_dy = bookHisto1D(1, 2, 1);
        
        plots.h_jetveto_mjj_veto = bookHisto1D("jetveto_mjj_baseline_veto", refData(1,1,2));
        plots.h_jetveto_mjj_inc = bookHisto1D("jetveto_mjj_baseline_inc", refData(1,1,2));
        plots.h_jetveto_dy_veto = bookHisto1D("jetveto_dy_baseline_veto", refData(1,2,2));
        plots.h_jetveto_dy_inc = bookHisto1D("jetveto_dy_baseline_inc", refData(1,2,2));

        plots.h_ptbaleff_mjj_veto = bookHisto1D("ptbaleff_mjj_baseline_veto", refData(1,1,3));
        plots.h_ptbaleff_mjj_inc = bookHisto1D("ptbaleff_mjj_baseline_inc", refData(1,1,3));
        plots.h_ptbaleff_dy_veto = bookHisto1D("ptbaleff_dy_baseline_veto", refData(1,2,3));
        plots.h_ptbaleff_dy_inc = bookHisto1D("ptbaleff_dy_baseline_inc", refData(1,2,3));

        plots.p_avgnjets_mjj = bookProfile1D(1,1,4);
        plots.p_avgnjets_dy = bookProfile1D(1,2,4);
      }

      if(phase_space=="highpt") {
        plots.h_mjj = bookHisto1D(2, 1, 1);
        plots.h_dy = bookHisto1D(2, 2, 1);
        
        plots.h_jetveto_mjj_veto = bookHisto1D("jetveto_mjj_highpt_veto", refData(2,1,2));
        plots.h_jetveto_mjj_inc = bookHisto1D("jetveto_mjj_highpt_inc", refData(2,1,2));
        plots.h_jetveto_dy_veto = bookHisto1D("jetveto_dy_highpt_veto", refData(2,2,2));
        plots.h_jetveto_dy_inc = bookHisto1D("jetveto_dy_highpt_inc", refData(2,2,2));

        plots.h_ptbaleff_mjj_veto = bookHisto1D("ptbaleff_mjj_highpt_veto", refData(2,1,3));
        plots.h_ptbaleff_mjj_inc = bookHisto1D("ptbaleff_mjj_highpt_inc", refData(2,1,3));
        plots.h_ptbaleff_dy_veto = bookHisto1D("ptbaleff_dy_highpt_veto", refData(2,2,3));
        plots.h_ptbaleff_dy_inc = bookHisto1D("ptbaleff_dy_highpt_inc", refData(2,2,3));

        plots.p_avgnjets_mjj = bookProfile1D(2,1,4);
        plots.p_avgnjets_dy = bookProfile1D(2,2,4);
      }

      if(phase_space=="search") {
        plots.h_mjj = bookHisto1D(3,1,1);
        plots.h_dy = bookHisto1D(3,2,1);
      }

      if(phase_space=="control") {
        plots.h_mjj = bookHisto1D(4,1,1);
        plots.h_dy = bookHisto1D(4,2,1);
      }

      if(phase_space=="highmass") {
        plots.h_njets = bookHisto1D(5, 3, 1);
        plots.h_dphijj = bookHisto1D(5, 4, 1);
        plots.h_ptbal = bookHisto1D(5, 5, 1);
      }
         
      return;
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      vector<ClusteredLepton> muons = applyProjection<DressedLeptons>(event, "DRESSED_MUONS").clusteredLeptons();
      vector<ClusteredLepton> electrons = applyProjection<DressedLeptons>(event, "DRESSED_ELECTRONS").clusteredLeptons();

      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt(25.0*GeV, MAXDOUBLE, -4.4, 4.4, RAPIDITY);

      // Make sure that we have a Z-candidate:
      FourMomentum* lep1 = 0;
      FourMomentum* lep2 = 0;
      if(muons.size() == 2) {
        FourMomentum dimuon = muons.at(0).momentum() + muons.at(1).momentum();
        if( dimuon.mass() > 81.0 &&
            dimuon.mass() < 101.0 &&
            PID::charge(muons.at(0)) != PID::charge(muons.at(1)) ) {
          lep1 = new FourMomentum(muons.at(0).momentum());
          lep2 = new FourMomentum(muons.at(1).momentum());
        }
      }
      if(electrons.size() == 2) {
        FourMomentum dielectron = electrons.at(0).momentum() + electrons.at(1).momentum();
        if( dielectron.mass() > 81.0 &&
            dielectron.mass() < 101.0 &&
            PID::charge(electrons.at(0)) != PID::charge(electrons.at(1)) ) {
          if(!lep1 && !lep2){
            lep1 = new FourMomentum(electrons.at(0).momentum());
            lep2 = new FourMomentum(electrons.at(1).momentum());
          } else { MSG_INFO("Found Z-candidate using both electrons and muons! Continuing with the muon-channel candidate"); }
        }
      }

      // If there's no Z-candidate, we wont use this event:
      if(!lep1 || !lep2) vetoEvent;

      // Do lepton-jet overlap removal:
      vector<const Jet*> good_jets;
      foreach(const Jet& j, jets) {
        bool nearby_lepton = false;
        foreach(ClusteredLepton& m, muons) {
          double dR = deltaR(j.momentum(), m.momentum());
          if(dR < 0.3) nearby_lepton=true;
        }
        foreach(ClusteredLepton& e, electrons) {
          double dR = deltaR(j.momentum(), e.momentum());
          if(dR < 0.3) nearby_lepton=true;
        }
        if(!nearby_lepton) good_jets.push_back(&j);
      }

      // If we don't have at least 2 good jets, we wont use this event.
      if(good_jets.size() < 2) vetoEvent;

      vars = new Variables(good_jets, lep1, lep2);
      bool pass_baseline = (vars->m_jet1pt > 55.0*GeV && vars->m_jet2pt > 45.0*GeV);
      bool pass_highpt = (vars->m_jet1pt > 85.0*GeV && vars->m_jet2pt > 75.0*GeV);
      bool pass_highmass = (pass_baseline && vars->m_mjj > 1000.0*GeV);
      bool pass_search = (pass_baseline && vars->m_zpt > 20.0*GeV && vars->m_ngapjets == 0 && vars->m_ptbalance2 < 0.15 && vars->m_mjj > 250.0*GeV);
      bool pass_control = (pass_baseline && vars->m_zpt > 20.0*GeV && vars->m_ngapjets > 0 && vars->m_ptbalance3 < 0.15 && vars->m_mjj > 250.0*GeV);

      if(pass_baseline) FillPlots(vars, baseline_plots, "baseline", weight);
      if(pass_highpt) FillPlots(vars, highpt_plots, "highpt", weight);
      if(pass_highmass) FillPlots(vars, highmass_plots, "highmass", weight);
      if(pass_search) FillPlots(vars, search_plots, "search", weight);
      if(pass_control) FillPlots(vars, control_plots, "control", weight);

      delete vars;
      if(lep1) delete lep1;
      if(lep2) delete lep2;

    }

    void FillPlots(Variables* vars, Plots& plots, std::string phase_space, double weight) {
      if(phase_space == "baseline" || phase_space == "highpt" || phase_space == "search" || phase_space == "control") {
        plots.h_dy->fill(vars->m_deltay, weight);
        plots.h_mjj->fill(vars->m_mjj, weight);
      }

      if(phase_space == "baseline" || phase_space == "highpt") {
        if(vars->m_pass_jetveto) {
          plots.h_jetveto_dy_veto->fill(vars->m_deltay, weight);
          plots.h_jetveto_mjj_veto->fill(vars->m_mjj, weight);
        }
        plots.h_jetveto_dy_inc->fill(vars->m_deltay, weight);
        plots.h_jetveto_mjj_inc->fill(vars->m_mjj, weight);

        if(vars->m_pass_ptbaleff) {
          plots.h_ptbaleff_mjj_veto->fill(vars->m_mjj, weight);
          plots.h_ptbaleff_dy_veto->fill(vars->m_deltay, weight);
        }
        plots.h_ptbaleff_mjj_inc->fill(vars->m_mjj, weight);
        plots.h_ptbaleff_dy_inc->fill(vars->m_deltay, weight);

        plots.p_avgnjets_dy->fill(vars->m_deltay, vars->m_ngapjets, weight);
        plots.p_avgnjets_mjj->fill(vars->m_mjj, vars->m_ngapjets, weight);
      }

      if(phase_space == "highmass") {
        plots.h_njets->fill(vars->m_ngapjets, weight);
        plots.h_dphijj->fill(vars->m_deltaphijj, weight);
        plots.h_ptbal->fill(vars->m_ptbalance2, weight);
      }

      return;
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      FinalizePlots(baseline_plots);
      FinalizePlots(highpt_plots);
      FinalizePlots(search_plots);
      FinalizePlots(control_plots);
      FinalizePlots(highmass_plots);
      FinalizeEfficiencies(baseline_plots);
      FinalizeEfficiencies(highpt_plots);
    }

    void FinalizePlots(Plots& plots) {
      if(plots.h_dy) normalize(plots.h_dy);
      if(plots.h_mjj) normalize(plots.h_mjj);
      if(plots.h_dphijj) normalize(plots.h_dphijj);
      if(plots.h_njets) normalize(plots.h_njets);
      if(plots.h_ptbal) normalize(plots.h_ptbal);
    }

    void FinalizeEfficiencies(Plots& plots) {
      int region_index = 0;
      if(plots.label=="baseline") region_index = 1;
      else if(plots.label=="highpt") region_index = 2;
      else return;

      if(plots.h_jetveto_mjj_veto && plots.h_jetveto_mjj_inc) divide(plots.h_jetveto_mjj_veto, plots.h_jetveto_mjj_inc, bookScatter2D(region_index, 1, 2));
      getScatter2D(region_index, 1, 2)->addAnnotation("InclusiveSumWeights", plots.h_jetveto_mjj_inc->integral());
      removeAnalysisObject(plots.h_jetveto_mjj_veto); removeAnalysisObject(plots.h_jetveto_mjj_inc);

      if(plots.h_jetveto_dy_veto && plots.h_jetveto_dy_inc) divide(plots.h_jetveto_dy_veto, plots.h_jetveto_dy_inc, bookScatter2D(region_index, 2, 2));
      getScatter2D(region_index, 2, 2)->addAnnotation("InclusiveSumWeights", plots.h_jetveto_dy_inc->integral());
      removeAnalysisObject(plots.h_jetveto_dy_veto); removeAnalysisObject(plots.h_jetveto_dy_inc);

      if(plots.h_ptbaleff_mjj_veto && plots.h_ptbaleff_mjj_inc) divide(plots.h_ptbaleff_mjj_veto, plots.h_ptbaleff_mjj_inc, bookScatter2D(region_index, 1, 3));
      getScatter2D(region_index, 1, 3)->addAnnotation("InclusiveSumWeights", plots.h_ptbaleff_mjj_inc->integral());
      removeAnalysisObject(plots.h_ptbaleff_mjj_veto); removeAnalysisObject(plots.h_ptbaleff_mjj_inc);

      if(plots.h_ptbaleff_dy_veto && plots.h_ptbaleff_dy_inc) divide(plots.h_ptbaleff_dy_veto, plots.h_ptbaleff_dy_inc, bookScatter2D(region_index, 2, 3));
      getScatter2D(region_index, 2, 3)->addAnnotation("InclusiveSumWeights", plots.h_ptbaleff_dy_inc->integral());
      removeAnalysisObject(plots.h_ptbaleff_dy_veto); removeAnalysisObject(plots.h_ptbaleff_dy_inc);
    }

    //@}


  private:

    Variables* vars;

    Plots baseline_plots;
    Plots highpt_plots;
    Plots search_plots;
    Plots control_plots;
    Plots highmass_plots;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1279489);

}
