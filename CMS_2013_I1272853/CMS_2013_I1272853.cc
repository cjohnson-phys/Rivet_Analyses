// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {


  class CMS_2013_I1272853 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1272853() : Analysis("CMS_2013_I1272853") {}


    /// Book histograms and initialise projections before the run
    void init() {
 
       const FinalState fs(-MAXRAPIDITY,MAXRAPIDITY);
       addProjection(fs, "FS");
    
       vector<pair<PdgId,PdgId> > vidsW;
       vidsW.push_back(make_pair(PID::MUON, PID::NU_MUBAR));
       vidsW.push_back(make_pair(PID::ANTIMUON, PID::NU_MU));

       FinalState fsW(-MAXRAPIDITY,MAXRAPIDITY);
       InvMassFinalState invfsW(fsW, vidsW, 20*GeV, 99999*GeV);
       addProjection(invfsW, "INVFSW");

       VetoedFinalState vfs(fs);
       vfs.addVetoOnThisFinalState(invfsW);
       addProjection(vfs, "VFS");
       addProjection(FastJets(vfs, FastJets::ANTIKT, 0.5), "Jets");


       _h_deltaS_eq2jet_Norm = bookHisto1D(1,1,1);
       _h_rel_deltaPt_eq2jet_Norm = bookHisto1D(2,1,1);

    } 


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      const InvMassFinalState& invMassFinalStateW = applyProjection<InvMassFinalState>(event, "INVFSW");
      bool isW(false);
      isW = (!(invMassFinalStateW.empty()));

      const ParticleVector& WDecayProducts = invMassFinalStateW.particles();
      if (WDecayProducts.size() < 2) vetoEvent;
      
      double pt1=-9999.,  pt2=-9999.;
      double phi1=-9999., phi2=-9999.;
      double eta1=-9999.;

      double mt = -9999; 

      int iNU_MU=-9999, iAN_MU=-9999;

      if (isW) {
        iNU_MU = (fabs(WDecayProducts[1].pdgId()) == PID::NU_MU) ? 1 : 0;
        iAN_MU = 1 - iNU_MU;
        pt1  = WDecayProducts[iAN_MU].momentum().pT();
        pt2  = WDecayProducts[iNU_MU].momentum().Et();
        eta1 = WDecayProducts[iAN_MU].momentum().eta();
        phi1 = WDecayProducts[iAN_MU].momentum().phi();
        phi2 = WDecayProducts[iNU_MU].momentum().phi();
        mt   = sqrt(2.0*pt1*pt2*(1.0-cos(phi1-phi2)));
      }
      
      if (!isW || mt < 50. || pt1 < 35. || fabs(eta1) > 2.1 || pt2 < 30.) vetoEvent; 


      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(20)) {
        if (fabs(jet.momentum().rapidity()) < 2.0) {
          jets.push_back(jet.momentum());
        }
      }
      
      if (jets.size() != 2) vetoEvent;

      double mupx    = pt1*cos(phi1);
      double mupy    = pt1*sin(phi1);
      double met_x   = pt2*cos(phi2);
      double met_y   = pt2*sin(phi2);

      double dpt = ((jets[0].px() + jets[1].px())*(jets[0].px() + jets[1].px()) + \
                    (jets[0].py() + jets[1].py())*(jets[0].py() + jets[1].py())); 
      double rel_dpt = sqrt(dpt)/ (jets[0].pT() + jets[1].pT());
         
      double pT2 = (mupx + met_x) * (mupx + met_x) + \
                   (mupy + met_y) * (mupy + met_y); 
      double Px       = (mupx + met_x)*(jets[0].px() + jets[1].px());
      double Py       = (mupy + met_y)*(jets[0].py() + jets[1].py());
      double p1p2_mag = sqrt(dpt) * sqrt(pT2);
      double dS       = acos((Px+Py) / p1p2_mag);

      _h_rel_deltaPt_eq2jet_Norm->fill(rel_dpt,weight);
      _h_deltaS_eq2jet_Norm->fill(dS,weight); 

    } 




    /// Normalise histograms etc., after the run
    void finalize() {

      double rel_dpt_bw = 1.0002 / 30.0;
      double dphi_bw = 3.14160 / 30.0;

      normalize(_h_rel_deltaPt_eq2jet_Norm, 1.*rel_dpt_bw);
      normalize(_h_deltaS_eq2jet_Norm, 1.*dphi_bw);

    }




  private:
    
    Histo1DPtr _h_rel_deltaPt_eq2jet_Norm;
    Histo1DPtr _h_deltaS_eq2jet_Norm;

  };




  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1272853);

}
