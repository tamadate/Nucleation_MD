# Nucleation_MD
This is a molecular dynamics (MD) simulation code for the modeling of the heterogeneous nucleation of vapor molecules with a low concentration molecule or ion existing in gas phase.  One of the feature of this code is the feasibility of the calculation at extremely low vapor pressure (we tested ~10 Pa of vapor pressure with MeOH and water).
## Overview
Figure 1 is a schematic diagram of this MD simulation code.  A molecular specie is arranged at the center of cubic calculation domain as a seed of the heterogeneous vapor nucleation, while the gas molecules and the vapor moecules are arranged entire of the cubic domain. These moeluces are treated as a all-atoms model in a spherical "effective domain" which is centerd on a seed molecule (which diameter should be large enough and 10 nm is currently employed).  The effective domain moves with centered seed molecule during whole simulation time, hence the seed molecule is always calclated as all-atom model. On the other hand, the gas molecule and the vapor molecule are hadlied as a point of mass at the outside of the effective domain.  In addition, every interactions at the outside of effective domain are ignored.  
## Usage

## Input script commands
### Ion input file
>**Syntax**:  `Input  fileName`  
**Example**:  `Input  angiotensinII+1.atom`<br><br>
>Ion input file name `fileName` is mentioned via the `Input` command, where the file include following components.
> - What atoms compose the ion
> - The atoms charges *i.e., particle charges*
> - The atoms position *i.e., ion initial structure*
> - Types of the atoms (potential parametes)
> - Bond types (potential parametes)
> - Angle types (potential parametes)
> - Dihedral types (potential parametes)
> - Bond list
> - Angle list
> - Dihedral list  
>  
> This section explain more detail about the ion input file and you can see some expamples from here.

### Vapor input file
>**Syntax**:  `Vapor	fileName	N`  
**Example**:  `Vapor	MeOH	100`<br><br>
>Vapor input file name (`fileName`) and number of vapors (`N`) are mentioned via `Vapor` command.  The input file structure is same as ion input.  Example arange 100 MeOH molecules in calculation doamin.
### Gas setting
>**Syntax**:  `Gas	type	N`  
**Example**:  `Gas	N2	100000`<br><br>
>Gas type (`type`) and number of its gas molecules (`N`) are mentioned via `Gas` command.  Currently, following four types of gases are available.  
>`He`: helium  
>`Ar`: Argon  
>`N2`:Nitrogen (diatomic)  
>`N2mono`: Nitrogen (monoatomic)  
Example arange 100,000 diatomic nitrogen molecules in calculation domain.
### Temperature
>**Syntax**:  `Temperature	T`  
**Example**:  `Temperature  300`<br><br>
>`Temperature` command set the calculation temperature (`T`) in the unit of Kelvin.  In this simulator, only gas & vapor temperatures are controled by velocity resampling method at the boundary.  The ion temperature is also supposed to be near this temperature by heat transfer with gas & vapor molecules.  Please refer here about the detail of this method.
### Pressure setting
>**Syntax**:  `Pressure	p`  
**Example**:  `Pressure 100000`<br><br>
>`Pressure` command set the gas pressure however this calculation does not control the pressure.  Also, virial is not considered in this pressure setting since this simulator is designed for diluted system, i.e., gas phase calculation.  If the molecule concentration is sufficiently low, the molecular interactions and molecular volumes are negligible, which means the perfect gas theory $PV=Nk_bT$ is avairable.  When the pressure, temperature, and number of gas molecule are given, the simulation volume is decided.  The simulation performed with that calculated volume (with constant T & constant N) is supposed to show the setting pressure.
### Time step
>**Syntax**:  `dt value`  
**Example**:  `dt 1`<br><br>
>`dt` command set the time step `value` in unit of fs.
### Total time steps
>**Syntax**:  `TotalSteps N`  
**Example**:  `TotalSteps 1000000000`<br><br>
>`Totalsteps` command set the total number of calculation steps (iteration), `N`.
### Relaxation time steps
>**Syntax**:  `RelaxSteps N`  
**Example**:  `RelaxSteps 1000000`<br><br>
>`RelaxSteps` command set the total number of relaxation steps which performed in prior of the main calculation to make ion structure thermally relax.
### Thermal bath settings (for relaxation)
>**Syntax**:  `NVTion	value`  
**Example**:  `NVTion	value`<br><br>
>`NVTion` command set the ion temperature, `T` at the thermal relaxation step.  If the relaxation process is performed without thermal bath, `OFF` is sbstituted in `value` in stead of temperature value.
### Output file setting
>**Syntax**:  `Output fileName  N`  
**Example**:  `output NaCl  100000`<br><br>
>`Output` command set a output file name, `fileName` and output interval, `N`.  Output file is generated with readable format via visualization tool [OVITO](https://www.ovito.org/) and the atomic positions/velocities of atoms in the effective domain are recorded.  The extension `.dump` is automatically given to `fileName`. In this example, the simulation generate `NaCl.dump` and write the atomic informaiton in each `N` steps.
### Atomic interactions
>**Syntax**:
>```
>Interactions
>   ion	potential
>   gg	potential
>   gi	potential
>   gv	potential
>   vi	potential
>   vv	potential
>```
>**Example**:
>```
>Interactions
>   ion	Born-Mayer-Huggins-NaCl
>   gg	OFF
>   gi	LJ
>   gv	LJ
>   vi	LJcoul
>   vv	LJcoul
>```
>`Interactions` command set the atomic interactions.


## Author
* Dr. Tomoya Tamadate
* [LinkedIn](https://www.linkedin.com/in/tomoya-tamadate-953673142/)/[ResearchGate](https://www.researchgate.net/profile/Tomoya-Tamadate)/[Google Scholar](https://scholar.google.com/citations?user=XXSOgXwAAAAJ&hl=ja)
* University of Minnesota
* tamalab0109[at]gmail.com
