# scs_frequency_collaterals_and_dcn

Code to generate models of dorsal column nucleus (DCN) neurons and primary afferent terminations within the dorsal horn, 
as described in "Biophysics of Frequency-Dependent Variation in Paresthesia and Pain Relief during Spinal Cord Stimulation" by
Rogers, Capogrosso, and Lempka (JNeurosci, 2024) (https://www.jneurosci.org/content/44/26/e2199232024).

To create a DCN cell (in python):
from DCN import * 
dcn_Cell = DCN()

##### To create a collateral: 
xyz positions of 25 collaterals are given in .xlsx format (here, zipped). Note, some collaterals were divided into 2 or 3 separate folders, 
///// and should be recombined before use. The accompanying python script can be used to build collaterals from these xyz data.
