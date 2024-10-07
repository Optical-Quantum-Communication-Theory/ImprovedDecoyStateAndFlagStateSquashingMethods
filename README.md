# Paper release Template (Replace with the title of the paper.)

This is a public version of the code used in Improved Decoy-state and Flag-state Squashing Methods \[[Arxiv](https://arxiv.org/abs/2405.05069)\]. This was built for a modified version of openQKDsecurity v1.0 and all changes have been included in this repository.


## BB84: 
`Main_BB84`: main file for any version of the BB84 protocol

### no decoy:
- `pmBB84WCP_nodecoy`: preset for decoy analysis without decoy intensities
- `pmBB84WCPChannel_Choi_nodecoy`: channel for no decoy intensities using improved decoy methods
- `pmBB84WCPChannel_constr_nodecoy`: channel for no decoy intensities using results from Rusca et. al., i.e. e_0 = 1/2

### decoy:
- `pmBB84WCP_decoy`: preset for BB84 protocol with decoy states
- `pmBB84WCPChannel_Choi`: channel with decoy intensities using improved decoy methods 

## Six-State:

### Qubit squasher from Gittsovich et. al.:
`Main_SixState`: main file for the six-state protocol with the squashing map from Gittsovich et. al.

- `pmSixStateWCP_decoy`: preset for six-state protocol with decoy state
- `pmSixStateLossyDescription`: Description file for six-state protocol with decoy state
- `pmSixStateWCPChannel`: Channel model for six-state protocol implementing decoy analysis from Wang et.al.
- `pmSixStateWCPChannel_improved_decoy`: Channel model for six-state protocol implementing improved decoy methods

### Flag-state squasher:
`Main_SixState_Flag`: main file for the six-state protocol with the flag-state squasher from Theorem 1 together with Eq. (78)

#### equal intensities:
- `pmSixStateWCP_decoy_Flag`: preset for six-state protocol with decoy state with equal intensities (clean and correct intervals)
- `pmSixStateLossyDescription_Flag_Reduced`: description file for six-state protocol with decoy states and the flag-state squasher from Theorem 1 together with Eq. (78)
- `pmSixStateWCPChannel_Flag_Reduced`: Channel model for six-state protocol implementing decoy analysis from Wang et. al. with flag-state squasher
- `pmSixStateWCPChannel_Flag_improved`: Channel model for six-state protocol using improved decoy methods with flag-state squasher


#### different intensities:
- `pmSixStateWCP_decoy_Flag_dif_int`: preset for six-state protocol with decoy state with different intensities (used with and without bit bias)
- `pmSixStateLossyDescription_Flag_dif_int`: description file for six-state protocol with decoy states and different intensities
- `pmSixStateLossyDescription_Flag_bit_bias`: description file for six-state protocol with decoy states and different intensities applying a bit bias
- `pmSixStateWCPChannel_Flag_dif_int`: Channel model for six-state protocol using improved decoy methods or decoy analysis from Wang et. al. for different intensities
- `pmSixStateWCPChannel_Flag_bit_bias`: Channel model for six-state protocol using improved decoy methods or decoy analysis from Wang et. al. for different intensities and applies a bit bias


## Install instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

1. Clone or download the repository.
2. Follow all further [install instructions](/openQKDsecurityV1/README.md).
3. Also follow the additional Mosek install instructions if you want an exact match.
3. \<Install directions for this repository. For example, add this folder to the Matlab path and save. Run this test function. Etc.\>
