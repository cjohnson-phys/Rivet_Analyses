// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Spherocity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/FParameter.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {

  using namespace Cuts;
  
  class EVENT_SHAPE_VAR : public Analysis {
  public:

    /// Constructor
    EVENT_SHAPE_VAR()
      : Analysis("EVENT_SHAPE_VAR") {    }

    /// Book histograms and initialise projections before the run
    void init() {
		//Should I project out the things that should be vetoed from the fastjet calc?
		//Probably not since I am not using an Identified Final State
	    // Projections
        FinalState fs;
        WFinder wfinder(fs, etaIn(-2.5, 2.5) & (pT >= 25.0*GeV), PID::MUON, 0.0*GeV, 100000000.0*GeV, 0.0*GeV, 0.2);
        addProjection(wfinder, "WFinder");
		
	    const FastJets jets( wfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
		//jets.useInvisibles();
	    addProjection(jets, "Jets");

	    // Book histograms
	    _hist_Thrust_180_350  = bookHisto1D("Thrust_180_350",20,0.0,1.0);
	    _hist_ThrustM_180_350  = bookHisto1D("ThrustMinor_180_350",20,0.0,1.0);
		_hist_Spherocity_180_350 = bookHisto1D("Spherocity_180_350",20,0.0,1.0);
		_hist_Sphericity_180_350 = bookHisto1D("Sphericity_180_350",20,0.0,1.0);
		_hist_F_180_350 = bookHisto1D("F_180_350",20,0.0,1.0);
		_hist_Sum_Mass_180_350 = bookHisto1D("Sum_Mass_180_350",20,0.0,1.0);
		_hist_Heavy_Mass_180_350 = bookHisto1D("Heavy_Mass_180_350",20,0.0,1.0);
		_hist_Tot_Broadening_180_350 = bookHisto1D("Tot_Broadening_180_350",20,0.0,1.0);
		_hist_Wide_Broadening_180_350 = bookHisto1D("Wide_Broadening_180_350",20,0.0,1.0);

	    _hist_Thrust_350_750 = bookHisto1D("Thrust_350_750",20,0.0,1.0);
	    _hist_ThrustM_350_750 = bookHisto1D("ThrustMinor_350_750",20,0.0,1.0);
		_hist_Spherocity_350_750 = bookHisto1D("Spherocity_350_750",20,0.0,1.0);
		_hist_Sphericity_350_750 = bookHisto1D("Sphericity_350_750",20,0.0,1.0);
		_hist_F_350_750 = bookHisto1D("F_350_750",20,0.0,1.0);
		_hist_Sum_Mass_350_750 = bookHisto1D("Sum_Mass_350_750",20,0.0,1.0);
		_hist_Heavy_Mass_350_750 = bookHisto1D("Heavy_Mass_350_750",20,0.0,1.0);
		_hist_Tot_Broadening_350_750 = bookHisto1D("Tot_Broadening_350_750",20,0.0,1.0);
		_hist_Wide_Broadening_350_750 = bookHisto1D("Wide_Broadening_350_750",20,0.0,1.0);

	    _hist_Thrust_750_up = bookHisto1D("Thrust_750_up",20,0.0,1.0);
	    _hist_ThrustM_750_up = bookHisto1D("ThrustMinor_750_up",20,0.0,1.0);
		_hist_Spherocity_750_up = bookHisto1D("Spherocity_750_up",20,0.0,1.0);
		_hist_Sphericity_750_up = bookHisto1D("Sphericity_750_up",20,0.0,1.0);
		_hist_F_750_up = bookHisto1D("F_750_up",20,0.0,1.0);
		_hist_Sum_Mass_750_up = bookHisto1D("Sum_Mass_750_up",20,0.0,1.0);
		_hist_Heavy_Mass_750_up = bookHisto1D("Heavy_Mass_750_up",20,0.0,1.0);
		_hist_Tot_Broadening_750_up = bookHisto1D("Tot_Broadening_750_up",20,0.0,1.0);
		_hist_Wide_Broadening_750_up = bookHisto1D("Wide_Broadening_750_up",20,0.0,1.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    	const double weight = event.weight();
		const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV);	//p_T threshold of 30.0*GeV

		const double YStarMax = 0.5;
		const double YBoostMax = 0.5;
		const double DeltaY = 1.0;
		double YStar = fabs(jets[0].momentum().rapidity()-jets[1].momentum().rapidity())/2;
		double YBoost = (jets[0].momentum().rapidity()+jets[1].momentum().rapidity())/2;

    	if ( (jets.size() < 2) ||
			 (jets[0].momentum().pT() < 60.0*GeV) ||
			 (jets[1].momentum().pT() < 30.0*GeV) ||
			 (YStar > YStarMax) ||
			 (YBoost > YBoostMax) ) { vetoEvent; }


		std::vector<Vector3> momenta;
		vector<FourMomentum> Jets;
		double Ht = 0.0;
		foreach (const Jet& j, jets) {
			if (fabs(j.momentum().eta() - YBoost) <= (YStarMax + DeltaY)) {
				Jets.push_back(j.momentum());
				Ht += j.momentum().pT();
				Vector3 mom = j.momentum().vector3();
				mom.setZ(0.0);
				momenta.push_back(mom); // Vector with just Px and Py -> to calculate transverse ESVs.
			}
		}

		// Leading jet should account for more than one third of the Ht in the central region.
		if (Jets[0].pT() < Ht/3.0) {
			cout << "pT_1 was less than one third of H_T. H_T: " << Ht <<  ", pT_1: " << Jets[0].pT() << endl;
			vetoEvent;
		}

		if (momenta.size() == 2) {
			// We need to use a ghost so that Thrust.calc() doesn't return 1.
			momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
		}
		
		Thrust thrust;
		thrust.calc(momenta);						// Calculates thrust and thrust axis
		//const double T = 1.0 - thrust.thrust();
		//const double TM = thrust.thrustMajor();
		
		Sphericity sphericity;
		sphericity.calc(momenta);
		//const double Spheri = sphericity.transSphericity();
		
		Spherocity spherocity;
		spherocity.calc(momenta);					// Calculates spherocity and spherocity axis
		//const double Sphero = spherocity.transSpheroctity();
		
		Vector3 thrust_axis = thrust.axis1();		// Retrieves thrust axis
		
		// Define Sphericity, F-Parameter, Masses, Broadenings by hand.
		
		// Split the event in upper and lower hemispheres and prepare weights for jet broadenings.
		double qup = 0.0, qdown = 0.0;
		double yup = 0.0, ydown = 0.0;
		double phiup = 0.0, phidown = 0.0;
		double bup = 0.0, bdown = 0.0;
		for (size_t index=0; index<Jets.size(); ++index) {
			if (thrust_axis.dot(Jets[index].vector3()) > 0.0) {
				qup += Jets[index].pT();
				yup += Jets[index].pT()*Jets[index].rapidity();
				phiup += Jets[index].pT()*Jets[index].phi();
			}
			else {
				qdown += Jets[index].pT();
				ydown += Jets[index].pT()*Jets[index].rapidity();
				phidown += Jets[index].pT()*Jets[index].phi();
			}
		}
		
		// Initialize M, M_lin matrices and Invariant Mass Vectors to zero
		FourVector Mass_Inv_Up, Mass_Inv_Down;
		double M[2][2] = {{0}}, M_lin[2][2] = {{0}};
		for (size_t index=0; index<Jets.size(); ++index) {
			// Fill hists for Px, Py, Pz, E, pT, Y, Phi (to be defined)
			double E = Jets[index].E();
			double Px = Jets[index].px();
			double Py = Jets[index].py();
			double Pz = Jets[index].pz();
			double Pt = Jets[index].pT();
			double Y = Jets[index].rapidity(); //maybe use eta here...probably not a big deal
			double Phi = Jets[index].phi();
			
			// Fill M, M_lin and Inv_Mass vectors
			if (thrust_axis.dot(Jets[index].vector3()) > 0.0) {
				Mass_Inv_Up[0] += E;
				Mass_Inv_Up[1] += Px;
				Mass_Inv_Up[2] += Py;
				Mass_Inv_Up[3] += Pz;
				bup += (1.0/(2.0*Ht))*Pt*sqrt((Y - yup/qup)*(Y - yup/qup) + (Phi - phiup/qup)*(Phi - phiup/qup));
			}
			else {
				Mass_Inv_Down[0] += E;
				Mass_Inv_Down[1] += Px;
				Mass_Inv_Down[2] += Py;
				Mass_Inv_Down[3] += Pz;
				bdown += (1.0/(2.0*Ht))*Pt*sqrt((Y - ydown/qdown)*(Y - ydown/qdown) + (Phi - phidown/qdown)*(Phi - phidown/qdown));
			}
			
			// M_lin Calculation
			M_lin[0][0] += Px*Px/Pt;
			M_lin[0][1] += Px*Py/Pt;
			M_lin[1][0] += Px*Py/Pt;
			M_lin[1][1] += Py*Py/Pt;
			// M Calculation
			M[0][0] += Px*Px;
			M[0][1] += Px*Py;
			M[1][0] += Px*Py;
			M[1][1] += Py*Py;
			
		}

		// The lowest bin also includes the underflow:
		double shapes[9];
		// Retrieves thrust value
		const double T = 1-thrust.thrust();
		// Retrieves thrust minor value
		const double T_minor = thrust.thrustMinor();
		// Retrieves spherocity value
		const double Sphero = spherocity.spherocity();
		// Calculates Rho (jet masses)
		const double rhoup = Mass_Inv_Up.invariant()/(Ht*Ht);
		const double rhodown = Mass_Inv_Down.invariant()/(Ht*Ht);
		const double RHOsc = rhoup + rhodown;
		const double RHOhc = max(rhoup,rhodown);
		// Calculates Jet Broadenings
		const double BTC = bup + bdown;
		const double BWC = max(bup,bdown);
		// Calculates F-parameter
		const double TrF = M_lin[0][0] + M_lin[1][1];
		const double DetF = M_lin[0][0]*M_lin[1][1] - M_lin[0][1]*M_lin[1][0];
		const double Lambda_1_F = TrF/2.0 + sqrt((TrF*TrF)/4.0 - DetF);
		const double Lambda_2_F = TrF/2.0 - sqrt((TrF*TrF)/4.0 - DetF);
		const double F = min(Lambda_1_F,Lambda_2_F)/max(Lambda_1_F,Lambda_2_F);
		// Calculate Sphericity
		const double TrS = M[0][0] + M[1][1];
		const double DetS = M[0][0]*M[1][1] - M[0][1]*M[1][0];
		const double Lambda_1_S = TrS/2.0 + sqrt((TrS*TrS)/4.0 - DetS);
		const double Lambda_2_S = TrS/2.0 - sqrt((TrS*TrS)/4.0 - DetS);
		const double Spheri = 2.0*min(Lambda_1_S,Lambda_2_S)/(Lambda_1_S+Lambda_2_S);
		shapes[0]=T,shapes[1]=T_minor,shapes[2]=Sphero,shapes[3]=RHOsc,shapes[4]=RHOhc;
		shapes[5]=BTC,shapes[6]=BWC,shapes[7]=F,shapes[8]=Spheri;
		cout << shapes << endl;
		
		
		if (180.0*GeV <= Ht && Ht < 350.0*GeV) {
			_hist_Thrust_180_350->fill(T, weight);
			_hist_ThrustM_180_350->fill(T_minor, weight);
			_hist_Spherocity_180_350->fill(Sphero, weight);
			_hist_Sphericity_180_350->fill(Spheri, weight);
			_hist_F_180_350->fill(F,weight);
			_hist_Sum_Mass_180_350->fill(RHOsc,weight);
			_hist_Heavy_Mass_180_350->fill(RHOhc,weight);
			_hist_Tot_Broadening_180_350->fill(BTC,weight);
			_hist_Wide_Broadening_180_350->fill(BWC,weight);
		} else if (Ht < 750.0*GeV) {
			_hist_Thrust_350_750->fill(T, weight);
			_hist_ThrustM_350_750->fill(T_minor, weight);
			_hist_Spherocity_350_750->fill(Sphero, weight);
			_hist_Sphericity_350_750->fill(Spheri, weight);
			_hist_F_350_750->fill(F,weight);
			_hist_Sum_Mass_350_750->fill(RHOsc,weight);
			_hist_Heavy_Mass_350_750->fill(RHOhc,weight);
			_hist_Tot_Broadening_350_750->fill(BTC,weight);
			_hist_Wide_Broadening_350_750->fill(BWC,weight);
		} else if (750.0*GeV <= Ht) {
			_hist_Thrust_750_up->fill(T, weight);
			_hist_ThrustM_750_up->fill(T_minor, weight);
			_hist_Spherocity_750_up->fill(Sphero, weight);
			_hist_Sphericity_750_up->fill(Spheri, weight);
			_hist_F_750_up->fill(F,weight);
			_hist_Sum_Mass_750_up->fill(RHOsc,weight);
			_hist_Heavy_Mass_750_up->fill(RHOhc,weight);
			_hist_Tot_Broadening_750_up->fill(BTC,weight);
			_hist_Wide_Broadening_750_up->fill(BWC,weight);
		}
	}

    void finalize() {
		double factor = crossSection()/sumOfWeights();
		
		scale(_hist_Thrust_180_350, factor);
		scale(_hist_ThrustM_180_350, factor);
		scale(_hist_Spherocity_180_350, factor);
		scale(_hist_Sphericity_180_350, factor);
		scale(_hist_F_180_350, factor);
		scale(_hist_Sum_Mass_180_350, factor);
		scale(_hist_Heavy_Mass_180_350, factor);
		scale(_hist_Tot_Broadening_180_350, factor);
		scale(_hist_Wide_Broadening_180_350, factor);

		scale(_hist_Thrust_350_750, factor);
		scale(_hist_ThrustM_350_750, factor);
		scale(_hist_Spherocity_350_750, factor);
		scale(_hist_Sphericity_350_750, factor);
		scale(_hist_F_350_750, factor);
		scale(_hist_Sum_Mass_350_750, factor);
		scale(_hist_Heavy_Mass_350_750, factor);
		scale(_hist_Tot_Broadening_350_750, factor);
		scale(_hist_Wide_Broadening_350_750, factor);

		scale(_hist_Thrust_750_up, factor);
		scale(_hist_ThrustM_750_up, factor);
		scale(_hist_Spherocity_750_up, factor);
		scale(_hist_Sphericity_750_up, factor);
		scale(_hist_F_750_up, factor);
		scale(_hist_Sum_Mass_750_up, factor);
		scale(_hist_Heavy_Mass_750_up, factor);
		scale(_hist_Tot_Broadening_750_up, factor);
		scale(_hist_Wide_Broadening_750_up, factor);
	}

  private:
	
	Histo1DPtr _hist_Thrust_180_350;
    Histo1DPtr _hist_ThrustM_180_350;
	Histo1DPtr _hist_Spherocity_180_350;
    Histo1DPtr _hist_Sphericity_180_350;
	Histo1DPtr _hist_F_180_350;
	Histo1DPtr _hist_Sum_Mass_180_350;
	Histo1DPtr _hist_Heavy_Mass_180_350;
	Histo1DPtr _hist_Tot_Broadening_180_350;
	Histo1DPtr _hist_Wide_Broadening_180_350;

    Histo1DPtr _hist_Thrust_350_750;
    Histo1DPtr _hist_ThrustM_350_750;
	Histo1DPtr _hist_Spherocity_350_750;
    Histo1DPtr _hist_Sphericity_350_750;
	Histo1DPtr _hist_F_350_750;
	Histo1DPtr _hist_Sum_Mass_350_750;
	Histo1DPtr _hist_Heavy_Mass_350_750;
	Histo1DPtr _hist_Tot_Broadening_350_750;
	Histo1DPtr _hist_Wide_Broadening_350_750;

    Histo1DPtr _hist_Thrust_750_up;
    Histo1DPtr _hist_ThrustM_750_up;
	Histo1DPtr _hist_Spherocity_750_up;
    Histo1DPtr _hist_Sphericity_750_up;
	Histo1DPtr _hist_F_750_up;
	Histo1DPtr _hist_Sum_Mass_750_up;
	Histo1DPtr _hist_Heavy_Mass_750_up;
	Histo1DPtr _hist_Tot_Broadening_750_up;
	Histo1DPtr _hist_Wide_Broadening_750_up;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EVENT_SHAPE_VAR);

}
