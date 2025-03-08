# LC-pericyte-project
# Below are notes of MATLAB codes for cerebral blood flow modeling;

# Funtions include 1) bfmodel.m that calculates cerebral blood flow based on the diameter dynamics of vessels;
#                  2) dfvsfo.m that calculates changes in blood flow;
#                  3) fwhm.m that calculates the full width at half maximum of blood flow changes;
#                  4) shadedErrorBar.m (from https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar) to plot the SEM of population data;

# Main files:
# 1) CBFmodeling_LPCcomponent.m for dissecting the contribution of direct and indirect LPC components to NVC regulation; load 2024_workspace_CBFmodeling.mat that contains the diameter dynamics of NVC and non-NVC vessels under LC ablation or RG ablation conditions;
# 2) CBFmodeling_LPCintensity.m for studying how LPC of different intensities modulate NVC flow; load 2024_workspace_LPCintensity.mat that contains the diameter dynamics of NVC and non-NVC vessels under mild, medium, and strong LPC;
# 3) CBFmodeling_LPCpattern.m for studying how LPC of different patterns (varied latency between capillary and arteriole-end vessel constrictions); load 2024_workspace_LPCpattern.mat that contains the diameter dynamics of NVC and non-NVC vessels under LPC of different patterns.
