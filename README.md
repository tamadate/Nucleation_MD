# Nucleation_MD
This is a molecular dynamics (MD) simulation code for the modeling of the heterogeneous nucleation of vapor molecules with a low concentration molecule or ion existing in gas phase.  One of the feature of this code is the feasibility of the calculation at extremely low vapor pressure (we tested ~10 Pa of vapor pressure with MeOH and water).
## Overview
Figure 1 is a schematic diagram of this MD simulation code.  A molecular specie is arranged at the center of cubic calculation domain as a seed of the heterogeneous vapor nucleation, while the gas molecules and the vapor moecules are arranged entire of the cubic domain. These moeluces are treated as a all-atoms model in a spherical "effective domain" which is centerd on a seed molecule (which diameter should be large enough and 10 nm is currently employed).  The effective domain moves with centered seed molecule during whole simulation time, hence the seed molecule is always calclated as all-atom model. On the other hand, the gas molecule and the vapor molecule are hadlied as a point of mass at the outside of the effective domain.  In addition, every interactions at the outside of effective domain are ignored.  
## Usage

## Input script commands
### Ion input file
>**Syntax**:  `Input  fileName`  
**Example**:  `Input  angiotensinII+1.atom`<br><br>
>Ion input file is mentioned via the `Input` command, where the file include following components.
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
>
### Gas setting
>```
>Input  fileName
>```
### Temperature setting
>```
>Temperature T
>```
### Pressure setting
>```
>Pressure	p
>```
* Time step
>```
>dt  value
>```
* Total time steps
>```
>TotalSteps	N
>```
* Relaxation time steps
>```
>RelaxSteps	N
>```
* Thermal bath settings (for relaxation)
>```
>NVTion	Tion or OFF
>```
* Output interval
>```
>Output	fileName	N
>```
* Atomic interactions
>```
>Interactions
>	ion	potential
>	gg	potential
>	gi	potential
>	gv	potential
>	vi	potential
>	vv	potential
>```
Born-Mayer-Huggins-NaCl

## Author
* Dr. Tomoya Tamadate
* [LinkedIn](https://www.linkedin.com/in/tomoya-tamadate-953673142/)/[ResearchGate](https://www.researchgate.net/profile/Tomoya-Tamadate)/[Google Scholar](https://scholar.google.com/citations?user=XXSOgXwAAAAJ&hl=ja)
* University of Minnesota
* tamalab0109[at]gmail.com
