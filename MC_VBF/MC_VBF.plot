#### Lepton Histos ------------------------------------------------------------------------

# BEGIN PLOT /MC_VBF/Lepton_phi
FullRange=1
Title=$\phi$ of highest $p_{\perp}$ $e^{\pm}$
XLabel=$\phi_{lepton}$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Lepton_eta
FullRange=1
Title=$\eta$ of highest $p_{\perp}$ $e^{\pm}$ 
XLabel=$\eta_{lepton}$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Lepton_pT
FullRange=1
Title=$p_{\perp}$ of highest $p_{\perp}$ $e^{\pm}$ 
XLabel=$lepton_{p_{\perp}}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}lepton_{p_{\perp}}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/SecLep_phi
FullRange=1
Title=$\phi$ of second highest $p_{\perp}$ $e^{\pm}$
XLabel=$\phi_{lepton}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{Second lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/SecLep_eta
FullRange=1
Title=$\eta$ of second highest $p_{\perp}$ $e^{\pm}$
XLabel=$\eta_{lepton}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{Second lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/SecLep_pT
FullRange=1
Title=$p_{\perp}$ of second highest $p_{\perp}$ $e^{\pm}$
XLabel=$lepton_{p_{\perp}}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}Second lepton_{p_{\perp}}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Neutrino_phi
FullRange=1
Title=$\phi$ of 2nd highest $p_{\perp}$ $e^{\pm}$ or ${\nu}_{e}$ 
XLabel=$\phi_{lepton}$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Neutrino_eta
FullRange=1
Title=$\eta$ of 2nd highest $p_{\perp}$ $e^{\pm}$ or ${\nu}_{e}$ 
XLabel=$\eta_{lepton}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{lepton}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Neutrino_pT
FullRange=1
Title=$p_{\perp}$ of 2nd highest $p_{\perp}$ $e^{\pm}$ or ${\nu}_{e}$ 
XLabel=$lepton_{p_{\perp}}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}lepton_{p_{\perp}}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Missing_eT
FullRange=1
Title=Missing Transverse Energy
XLabel=$E_{T}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}met$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Missing_eT_phi
FullRange=1
Title=$\phi$ of Missing Transverse Energy
XLabel=$\phi$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}met_{\phi}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dilepton_mT
FullRange=1
Title=Transverse Mass of di-lepton system
XLabel=$m^T_{ll}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m^T_{ll}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dilepton_mass
FullRange=1
Title=Di-lepton mass
XLabel=$m_{ll}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{ll}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/boson_pT
FullRange=1
Title=Boson $p_{\perp}$ 
XLabel=$p^W_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p^W_{\perp}$ [pb/GeV]
LogX=0
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/RATIOboson_pT
MainPlot=0
RatioPlotYMin=0
RatioPlotYMax=12
Title=RATIO of Boson $p_{\perp}$ 
XLabel=$p^W_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p^W_{\perp}$ [pb/GeV]
LogX=0
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/boson_pT_50
FullRange=1
Title=Boson $p_{\perp}$ alternate binning 
XLabel=$p^W_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p^W_{\perp}$ [pb/GeV]
LogX=0
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/njet_inc
FullRange=1
Title=Inclusive jet multiplicity
XLabel=$N_{\text{jet}}$
YLabel=$\sigma$ [pb]
XMajorTickMarks=10
XMinorTickMarks=0
RatioPlotYMin=0.45
RatioPlotYMax=1.57
##LogY=0
# END PLOT

#Jet Histos ------------------------------------------------------------------------ 

# BEGIN PLOT /MC_VBF/Jet1_phi
FullRange=1
Title=Jet 1 $\phi$
XLabel=$\phi$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{J1}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet1_eta
FullRange=1
Title=Jet 1 $\eta$
XLabel=$\eta$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{J1}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet1_pT
FullRange=1
Title=Jet 1 $p_{\perp}$
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}/J1 p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet2_phi
FullRange=1
Title=Jet 2 $\phi$
XLabel=$\phi$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{J2}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet2_eta
FullRange=1
Title=Jet 2 $\eta$
XLabel=$\eta$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{J2}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet2_pT
FullRange=1
Title=Jet 2 $p_{\perp}$
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}/J2 p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet3_phi
FullRange=1
Title=Jet 3 $\phi$
XLabel=$\phi$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\phi_{J3}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet3_eta
FullRange=1
Title=Jet 3 $\eta$
XLabel=$\eta$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_{J3}$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet3_pT
FullRange=1
Title=Jet 3 $p_{\perp}$
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}/J3 p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_Mass
FullRange=1
Title=Di-jet mass
XLabel=$m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{jj}$ [pb/GeV]
#LogY=0
GofType=chi2
GofLegend=1
Rebin=5
LegendXPos=0.3
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_Mass_50
FullRange=1
Title=Di-jet mass alternate binning
XLabel=$m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{jj}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_Mass_zoomed
FullRange=1  
Title=Di-jet mass
XLabel=$m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{jj}$ [pb/GeV]
#LogY=0
# END PLOT 

# BEGIN PLOT /MC_VBF/RATIODijet_Mass
MainPlot=0
RatioPlotYMin=0
RatioPlotYMax=2
Title=Di-jet mass ratio
XLabel=$m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{jj}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Scalar_sum_pT
FullRange=1
Title=Scalar $\Sigma$ $p_{\perp}$
XLabel=$\Sigma p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Vector_sum_pT
FullRange=1
Title=Vector $\Sigma$ $p_{\perp}$
XLabel=$\Sigma p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet_delta_rapidity
FullRange=1
Title=$\Delta$ Rapidity Jets
XLabel=$\Delta$ Rapidity
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta$ Rapidity [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet_delta_Phi
FullRange=1
Title=$\Delta \phi$ Jets
XLabel=$\Delta \phi$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta \phi$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Jet_delta_R
FullRange=1
Title=$\Delta R$ Jets
XLabel=$\Delta R$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta R$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Eta_product
FullRange=1
Title=Eta Product Jets 1,2
XLabel=$\eta_1 * \eta_2$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta_1 * \eta_2$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_pT
FullRange=1
Title=Dijet $p_{\perp}$
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_pT_50
FullRange=1
Title=Dijet $p_{\perp}$ alternate binning
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/RATIODijet_pT
MainPlot=0
RatioPlotYMin=0
RatioPlotYMax=7
Title=Dijet $p_{\perp}$
XLabel=$p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/GapRef_PtAvg
FullRange=1
Title=GapRef (avg jet $p_{\perp}$)
XLabel=$AVG p_{\perp}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}$ AVG $p_{\perp}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Opposite_Hemisphere
FullRange=1
Title=Jets in opposite hemispheres? (2 = yes)
XLabel=0-same hemi, 2-opposite
YLabel=$\sigma$ [pb]
LogY=0
# END PLOT

#Other Histos  ------------------------------------------------------------------------

# BEGIN PLOT /MC_VBF/Total_mT
FullRange=1
Title=Transverse Mass of entire final state
XLabel=$m^T$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m^T$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/_h_mllmjj
Title=Di-lepton mass vrs Di-jet mass
XLabel=$m_{jj} [GeV]$
YLabel=$m_{ll} [GeV]$
#LogY=0
LogX=0
RatioPlot=0
LineStyle=none
LineWidth=0pt
LogZ=1
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Delta_Phi_dijet_dilep
FullRange=1
Title=$\phi$ between di-jet and di-lepton
XLabel=$\Delta \phi$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta \phi$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Total_visible_Mass
FullRange=1
Title=Total jet mass
XLabel=$m_{\Sigma j_i}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{\Sigma j_i}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Total_VBF_Mass
FullRange=1
Title=Total VBF mass
XLabel=$m_{VBF}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}m_{VBF}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Delta_boson_Jet_Mass
FullRange=1
Title=Difference between boson mass and closest di-jet mass
XLabel=$\Delta m_{jj} m_{V}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta m_{jj} m_{V}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Closest_to_boson_Mass
FullRange=1
Title=Mass of di-jet closest to boson Mass
XLabel=$\Delta m_{jj}$ [GeV]
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta m_{jj}$ [pb/GeV]
#LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/Dijet_Eta*
FullRange=1
Title=Dijet $\eta$*
XLabel=$\eta$ 
YLabel=$\mathrm{d}\sigma/\mathrm{d}\eta$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/W_or_Z
FullRange=1
Title=Event considered as $W^-=0, Z=2, W^+=4$
XLabel=W$\pm$ or Z 
YLabel=$\sigma$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/dRLepJet
FullRange=1
Title=Delta R between leading Lepton and leading Jet
XLabel=$\Delta R$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Delta R$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/MomFacTot
FullRange=1
Title=MomFacTot (X*Y)
XLabel=Momentum Variable
YLabel=$\mathrm{d}\sigma/\mathrm{d}MomVar$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/MomFacTot350
FullRange=1
Title=MomFacTot (X*Y) dijet_m $\gt$ 350 GeV
XLabel=Momentum Variable
YLabel=$\mathrm{d}\sigma/\mathrm{d}MomVar$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/MomFacTot1000
FullRange=1
Title=MomFacTot (X*Y) dijet_m $\gt$ 1 TeV
XLabel=Momentum Variable
YLabel=$\mathrm{d}\sigma/\mathrm{d}MomVar$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/BosonDijetPtSum
FullRange=1
Title=Scalar sum of dijet and boson pT
XLabel=$\Sigma p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/BosonDijetPtSum350
FullRange=1
Title=Scalar sum of dijet and boson pT, dijet_m $\gt$ 350GeV
XLabel=$\Sigma p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/BosonDijetPtSum1000
FullRange=1
Title=Scalar sum of dijet and boson pT, dijet_m $\gt$ 1 TeV
XLabel=$\Sigma p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/BosonDijetPtSum2Jets
FullRange=1
Title=Scalar sum of dijet and boson pT, 2 JETS, dijet_m $\gt$ 350 GeV
XLabel=$\Sigma p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d}\Sigma p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/Zeppenfeld
FullRange=1
Title=Zeppenfeld Variable
XLabel=$\eta \ast_3$
YLabel=$\mathrm{d}\sigma/\mathrm{d}Zeppenfeld$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/DPtInW
FullRange=1
Title=Dijet p_T in direction of Boson
XLabel=$p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d} p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/DPtInW350
FullRange=1
Title=Dijet p_T	in direction, dijet_m $\gt$ 350GeV
XLabel=$p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d} p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/DPtInW1000
FullRange=1
Title=Dijet p_T	in direction, dijet_m $\gt$ 1 TeV
XLabel=$p_{perp}$
YLabel=$\mathrm{d}\sigma/\mathrm{d} p_{perp}$ [pb]
# END PLOT

# BEGIN PLOT /MC_VBF/VBFAngle_cuts
FullRange=1
Title=Lepton Gap Centrality (VBF Cuts)
XLabel=VBF Angle
YLabel=$\mathrm{d}\sigma/\mathrm{d} VBF Angle$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/VBFAngle
FullRange=1
Title=Lepton Gap Centrality
XLabel=VBF Angle
YLabel=$\mathrm{d}\sigma/\mathrm{d} VBF Angle$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/VBFAngle_cuts_B
FullRange=1
Title=Boson Gap Centrality (VBF Cuts)
XLabel=VBF Angle
YLabel=$\mathrm{d}\sigma/\mathrm{d} VBF Angle_B$ [pb]
LogY=0
# END PLOT

# BEGIN PLOT /MC_VBF/VBFAngle_B
FullRange=1
Title=Boson Gap Centrality
XLabel=VBF Angle
YLabel=$\mathrm{d}\sigma/\mathrm{d} VBF Angle_B$ [pb]
LogY=0
# END PLOT
