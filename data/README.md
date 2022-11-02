# Data for the article "Latent functional diversity may accelerate microbial community responses to temperature fluctuations"

## Dataset Attribution and Usage

Dataset to reproduce the results in our eLife article "Latent functional diversity may accelerate microbial community responses to temperature fluctuations". Here we provide the raw trait data estimates from our experiment, alongside thermal performance curve parameter estimates.

* Persistent Identifier: https://doi.org/10.5061/dryad.f1vhhmh0g

* Dataset Contributors: Thomas P. Smith, Shorok Mombrikotb, Emma Ransome, Dimitrios-Georgios Kontopoulos, Samraat Pawar, Thomas Bell.

 * Methodological Information: Bacterial isolates were collected from a species sorting experiment. Isolates were grown at different temperatures, their growth, respiration and ATP content measured in order to produce thermal performance curves for these traits. See manuscript for full details.

-----------------

## Description of the Data and file structure

### Table of Contents:

* ATP-trait-data.csv
* growth_summary.csv
* summary_OD_R.csv
* summary_resp_R_biomass.csv

### Details for: ATP-trait-data.csv

* Description: raw trait measurements from all strains at all temperatures. These trait measurements were used for the Schoolfield-Sharpe model fitting.

* Format: .csv

* Size: 2.0 MB

* Variables (note: measurements were taken at two time-points, the start and end of the experment, these are denoted _start and _end respectively):
  * Well: well of 96-well plate used.
  * Plate: plate number, for where multiple plates were required for a given temperature condition.
  * Strain: unique identifier for bacterial strain (see manuscript for genbank accession numbers).
  * Temperature: experimental temperature (celsius).
  * Hours: duration of the experiment in hours.
  * cells_per_ml: bacterial cell concentration estimated by flow cytometry.
  * FSCH: forward scatter by flow cytometry.
  * Microresp: absorbance values for microresp plate.
  * ATP: lumiesence measurements from bactiter-glo assay.
  * OD: Optical density at 600nm.
  * cell_diameter: cell diameter (µm), estimated from FSCH measurements.
  * cell_volume: cell volume (µm^3), estimated from cell diameter measurements.
  * mean_diameter: mean diameter across replicates.
  * mean_volume: mean volume across replicates.
  * biomass: estimated biomass (µg carbon) based on cell number and mean size.
  * biomass_increase: change in biomass across duration of experiment.
  * CO2_respired: quantity of CO2 released (mg), per microresp measurements.
  * C_respired: quanity of carbon respired (µg).
  * ATP_per_ml: concentration of ATP in solution, based on luminesence measurements.
  * ATP_per_biomass: concentration of ATP per unit of biomass.
  * growth_rate; growth_rate_biomass: growth rate based on cell numbers, or biomass estimates.
  * Resp_rate; Resp_rate_per_cell; Resp_rate_per_biomass: respiration rate (carbon respired over time), expressed per cell and per unit biomass.
  * Phylum; Class; Order; Family; Genus; Species: taxonomic information.


### Details for: growth_summary.csv

* Description: Schoolfield-Sharpe model fits to growth-rate thermal performance curves, derived from flow cytometry measurements (see manuscript methods for full details).

* Format: .csv

* Size: 18.5 kB

* Variables:
  * strain: unique identifier for strain of bacteria (see manuscript for genbank accession numbers).
  * volume: estimated cell volume (µm^3).
  * E_sch; B0_sch; E_D_sch; T_pk_est_sch: E, B0, E_D and T_pk schoolfield model fitted parameter values - see manuscript for details. 
  * T_pk_sch: Temperature (Kelvin) with highest recorded growth, used as starting value for Tpk in model fitting.
  *	P_pk_sch: growth rate at fitted Tpk.
  * r_sq_sch: model R^2.
  * n_temps; temps_before_peak; temps_after_peak: number of temperatures with datapoints that the model was fitted to, number of temperatures before and after the fitted curve peak.
  * Phylum; Class; Order; Family; Genus; Species: taxonomic information.
  * iso_temp; incu_temp: isolation and incubation temperature (celsius), see manuscript for methodological details. RT is room temperature, 22C.
  * niche_width: width of Schoolfield curve (degress celsius) where growth is 50% of maximum.

### Detials for: summary_OD_R.csv

* Description: Schoolfield-Sharope model fits to growth-rate thermal performance curves derived from optical density measurements (see manuscript for details).

* Format: .csv

* Size: 944 bytes

* Variables:  
  * Variables as described for growth_summary.csv


### Details for: summary_resp_R_biomass.csv

* Description: Schoolfield-Sharpe model fits to mass-specific respiration rate thermal performance curves, derived from MicroResp measurements (see manuscript for full methodological details).

* Format: .csv

* Size: 14.5 kB

* Variables:  
  *  Variables as described for growth_summary.csv. Where NAs are present, this is because mass-specific respiration rates were unable to be measured for these strains (see manuscript for details).



## Sharing/access Information

Data and all code to perform our analyses are available in a github repository: https://github.com/smithtp/latent-diversity