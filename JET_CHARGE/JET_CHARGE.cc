// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"

size_t event_num = 0;
namespace Rivet {


  class JET_CHARGE : public Analysis {
  public:

    /// Constructor
    JET_CHARGE()
      : Analysis("JET_CHARGE")
    {    }


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

      /// Book histograms
      _h_pdgid = bookHisto1D("PID", 501, -0.5, 500.5);
      _h_partP = bookHisto1D("partP", 100, 0, 100);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      event_num+=1;
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets_w");

      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(30.0*GeV)) {
        if (fabs(jet.momentum().rapidity()) > 4.4) continue;
        jets.push_back(jet.momentum());
      }

      if (jets.size() < 2) {
        MSG_INFO("Not a dijet event!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        vetoEvent;
      }
      
      size_t jet_num = 0;
      MSG_INFO("Event Number: " << event_num);
      foreach (const Jet& jet, jetpro.jetsByPt(30.0*GeV)) {
        jet_num+=1;
        const double jet_pt = jet.pT();
        MSG_INFO("Jet Number: " << jet_num << " Jet pT: " << jet_pt);
        GenParticle* common_ancestor = NULL;
        size_t qg_index = 0;
        size_t part_num = 0;
        double part_pt_sum = 0;
        foreach (const Particle& p, jet.particles()) {
          part_num += 1;
          const double part_pT = p.momentum().pT();
          part_pt_sum += part_pT;
          int curr_pid = p.pdgId();
          MSG_INFO("Particle Number: " << part_num << " Particle pT: " << part_pT << " PID: " << curr_pid);
          std::vector<GenParticle*> ancestors = particles_in(p);
          int anc_pid_below = 0;
          int anc_pid_above = 0;
          double max_below_jet_pt = 0.0;
          double min_above_jet_pt = 10000000000000;
          for (std::vector<int>::size_type i = ancestors.size() - 1; i != (std::vector<int>::size_type) -1; i--){
              double ancestor_pt = ancestors[i]->momentum().perp();
              if (ancestor_pt > max_below_jet_pt && ancestor_pt < jet_pt) {
                  max_below_jet_pt = ancestor_pt;
                  anc_pid_below = ancestors[i]->pdg_id();
              }
              if (ancestor_pt < min_above_jet_pt && ancestor_pt > jet_pt) {
                  min_above_jet_pt = ancestor_pt;
                  anc_pid_above = ancestors[i]->pdg_id();
              }
              // if (ancestors[i]->momentum().perp() == jet_pt){
              //     size_t ancestor_pid = ancestors[i]->pdg_id();
              //     MSG_INFO("Parent Particle: " << ancestor_pid);
              // }
              // if (fabs(ancestors[i]->pdg_id()) < 7 || ancestors[i]->pdg_id() == 21) {
              //     size_t ancestor_pid = ancestors[i]->pdg_id();
              //     MSG_INFO("Quark/Gluon: " << ancestor_pid);
              //     qg_index = i;
              //     break;
              // } 
          }
          MSG_INFO("PID(below): " << anc_pid_below << " " << max_below_jet_pt << " PID(above): " << anc_pid_above << " " << min_above_jet_pt);
        //   if (qg_index == 0) MSG_INFO("Ancestor: " << ancestors.back());
        //   else MSG_INFO("Ancestor: " << ancestors[qg_index]);
        //   size_t ancestor_pid = ancestors.back()->pdg_id();
        //   MSG_INFO("Ancestor PID: " << ancestor_pid);
        //   GenParticle* latest_ancestor = ancestors.back();
        }
        MSG_INFO("Sum of particle pT: " << part_pt_sum);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    Profile1DPtr _h_XXXX;
    Histo1DPtr _h_pdgid;
    Histo1DPtr _h_partP;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JET_CHARGE);

}
