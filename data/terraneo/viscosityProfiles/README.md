# Viscosity Profiles 

Within the directory `data/terraneo/viscosityProfiles` two viscosity profiles can be found in order to read in background viscosity profiles for simulation runs with a convection simulation application. This will allow the user to parse a background viscosity profile and estimate temperature dependent viscosity changes at depth. 
The viscosity profiles files are in `.json` format and contain two columns, the first one being `Radius [m]` and the second one being `Viscosity [Pa s]`.

For temperature (and depth) depdent viscosity estimations using the TerraNeo application a user can choose between 4 definitions within the corresponding parameters file: Frank-Kamenetskii type 1, Frank-Kamenetskii type 2, Arrhenius type and with respect to the mean temperature. 
The four temperature dependent laws are defined as follows:

## Frank-Kamenetskii type 1:
$\eta = \eta_0  e^{-E_A \cdot T + df \cdot (\frac{r_{Max} - r }{r_{Max} - r_{Min}})}$

## Frank-Kamenetskii type 2:
$\eta = \eta_0  e^{E_A \cdot (0.5 - T) + df \cdot (\frac{r_{Max} - r }{r_{Max} - r_{Min}})}$

## Arrhenius:
$\eta = \eta_0  e^{E_A \cdot (\frac{1}{T + 0.5} - 1) + df \cdot (\frac{r_{Max} - r }{r_{Max} - r_{Min}})}$

## With respect to the mean temperature $\overline{T}$:
$\eta = \eta_0  e^{-E_A \cdot ( T - \overline{T} )}$


With $\eta_0$ as the reference (background) viscosity, $E_A$ being the activation energy that determines the sensitivity of the viscosity to the temperature T (non-dimensionalised), `df` being a depth dependence factor accounting for a pressure-induced decrease in viscosity with depth, similar to the activation volume $V^*$ as mentioned in Lin et al. 2022. $r_{Max}$ and $r_{Min}$ describe the maximum and minimum radius of the domain, respectively. In case of a 3D Earth mantle convection simulation: $r_{Max} = r_{Surface}$ and $r_{Min} = r_{CoreMantleBoundary}$. 

## Viscosity Profile after Lin et al. 2022

The first viscosity profile `ViscosityProfile_Lin_et_al_2022.json` is similar to the one used by Lin et al. 2022 in combination with a Frank-Kamenetskii type 1 temperature dependent viscosity law with an activation energy of $E_A = 4.610$ and a depth dependence factor $df = 3$. 
This viscosity profile can be used to determine background viscosity values at depth in order to estimate the temperature dependent viscosity following the specified law.
A detailed description can be found in Lin et al. 2022: `DOI: 10.1029/2022GC010514`!

![Viscosity Profile Lin et al. 2022](doc/images/ViscosityProfile_Lin_et_al_2022.png "Viscosity profile Lin et al. 2022"){width=50%}


## Viscosity Profile after Stotz et al. 2017

The second viscosity profile `ViscosityProfile_Stotz_et_al_2017.json` corresponds to the viscosity profile utilized by Stotz et al. 2017 together with a Frank-Kamenentskii type 1 temperature dependence for $\eta$. The activation energy was defined with $E_A = 4.605$ and a depth dependence factor of $df = 3.976$ was chosen.
This viscosity profile accounts for a low viscous Asthenosphere of $\eta = 5 \cdot 10^{-19}$ Pa s with a channel thickness of approximately $150$ km. 
An elaboration of the utilized viscosity profile can be found in Stotz et al. 2017: `DOI: 10.1002/2017GL075697`

![Viscosity Profile Stotz et al. 2017 ](doc/images/ViscosityProfile_Stotz_et_al_2017.png "Viscosity profile Stotz et al. 2017"){width=50%}

