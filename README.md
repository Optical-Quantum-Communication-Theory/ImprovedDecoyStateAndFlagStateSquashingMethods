# Improved Decoy-state and Flag-state Squashing Methods

This is a public version of the code used in Improved Decoy-state and Flag-state Squashing Methods \[[Arxiv](https://arxiv.org/abs/2405.05069)\]. This was built for a modified version of openQKDsecurity v1.0 and all changes have been included in this repository.

## Installation Instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

1. Clone or download the repository.
2. Follow all further [install instructions](/openQKDsecurityV1/README.md).
3. Also follow the additional Mosek install instructions if you want an exact match.



## Figure Instructions

The functions below are required to generate the plots in the paper. For each plot we describe below which settings need to be chosen. For all plots the change from two to one decoy intensity is done as follows
1. Comment out the parameters with two decoy intensities in the preset file,
2. Comment out the variable names with two decoy intesities in the channel function,
3. Comment out the variable 'decoys' containing two decoy intensities.

Here are the specific minor modificaitons required to reproduce the figures in the paper.

### Figure 1
Requires functions `Main_BB84` and `pmBB84WCP_nodecoy` for the correct parameters. For comparing different decoy methods, in the `pmBB84WCP_nodecoy` preset change the channel model to `pmBB84WCPChannel_constr_nodecoy` for the yellow triangles, and to `pmBB84WCPChannel_Choi_nodecoy` for the red circles. 

### Figure 2
Requires functions `Main_BB84` and `pmBB84WCP_decoy` for the correct parameters and `pmBB84WCPChannel_Choi` (yellow stars), `pmBB84WCPChannel` (green diamonds and red circles) for different decoy-state methods.

### Figure 3
This uses the same data as in figure 2.

### Figure 7
Requires functions `Main_SixState` and `pmSixStateWCP_decoy` for the correct parameters of the simple squasher [2, 3] and `pmSixStateWCPChannel` or `pmSixStateWCPChannel_improved_decoy` for the decoy methods. Both will create the same (up to small differences) key rate curve (green crosses). `Main_SixState_Flag` and `pmSixStateLossyDescription_Flag_Reduced` with $p_z^A = \dots = p_y^A = p_z^B = \dots = p_y^B = \frac{1}{3}$ and the decoy state methods from `pmSixStateWCPChannel_Flag_Reduced` or `pmSixStateWCPChannel_Flag_improved` create the curve corresponding to the flag-state squasher (red diamonds).

### Figure 8
`Main_SixState_Flag` and `pmSixStateLossyDescription_Flag_Reduced` with $p_z^A= p_z^B = 0.8, p_y^A = p_y^B = p_x^A = p_x^B = 0.1$ and the decoy state methods from `pmSixStateWCPChannel_Flag_improved` create the curve corresponding to the flag-state squasher (red stars). The file `pmSixStateWCPChannel_Flag_Reduced` is used for the standard decoy methods, creating the curve with green circles.

### Figure 9
`Main_SixState_Flag` and `pmSixStateLossyDescription_Flag_Reduced` with $p_z^A= p_z^B = 0.8, p_y^A = p_y^B = p_x^A = p_x^B = 0.1$, commenting out all unncessary decoy intensities. The files `pmSixStateWCPChannel_Flag_improved` and  `pmSixStateWCPChannel_Flag_nodecoy` create the curves using the improved decoy methods from this work and applying $e_0 = 1/2$ from Rusca et. al. [1], respectively.

### Figure 10
`Main_SixState_Flag`, `pmSixStateWCP_decoy_Flag_dif_int` with $p_z^A = \dots = p_y^A = p_z^B = \dots = p_y^B = \frac{1}{3}$. The description file `pmSixStateLossyDescription_Flag_dif_int` and the channel `pmSixStateWCPChannel_Flag_dif_int` are required for the curves **without** a bit bias. The description file `pmSixStateLossyDescription_Flag_bit_bias` and the channel `pmSixStateWCPChannel_Flag_bit_bias` are required for the curves **with** a bit bias.

### Figure 11
Same as figure 10, but each key rate curve is divided by the sigle photon probability.

## Functions
### BB84: 
`Main_BB84`: main file for any version of the BB84 protocol

#### no decoy:
- `pmBB84WCP_nodecoy`: preset for decoy analysis without decoy intensities
- `pmBB84WCPChannel_Choi_nodecoy`: channel for no decoy intensities using improved decoy methods
- `pmBB84WCPChannel_constr_nodecoy`: channel for no decoy intensities using results from Rusca et. al. [1], i.e. e_0 = 1/2

#### decoy:
- `pmBB84WCP_decoy`: preset for BB84 protocol with decoy states
- `pmBB84WCPChannel_Choi`: channel with decoy intensities using improved decoy methods 

### Six-State:

#### Simple squasher from Refs [2, 3]:
`Main_SixState`: main file for the six-state protocol with the squashing map from Refs [2, 3]

- `pmSixStateWCP_decoy`: preset for six-state protocol with decoy state
- `pmSixStateLossyDescription`: Description file for six-state protocol with decoy state
- `pmSixStateWCPChannel`: Channel model for six-state protocol implementing decoy analysis from Wang et.al. [4]
- `pmSixStateWCPChannel_improved_decoy`: Channel model for six-state protocol implementing improved decoy methods

#### Flag-state squasher:
`Main_SixState_Flag`: main file for the six-state protocol with the flag-state squasher from Theorem 1 together with Eq. (78)

##### equal intensities:
- `pmSixStateWCP_decoy_Flag`: preset for six-state protocol with decoy state with equal intensities (clean and correct intervals)
- `pmSixStateLossyDescription_Flag_Reduced`: description file for six-state protocol with decoy states and the flag-state squasher from Theorem 1 together with Eq. (78)
- `pmSixStateWCPChannel_Flag_Reduced`: Channel model for six-state protocol implementing decoy analysis from Wang et. al. [4] with flag-state squasher
- `pmSixStateWCPChannel_Flag_improved`: Channel model for six-state protocol using improved decoy methods with flag-state squasher
- `pmSixStateWCPChannel_Flag_nodecoy`: Channel model for six-state protocol with no decoy intensities using results from Rusca et. al. [1], i.e. e_0 = 1/2


##### different intensities:
- `pmSixStateWCP_decoy_Flag_dif_int`: preset for six-state protocol with decoy state with different intensities (used with and without bit bias)
- `pmSixStateLossyDescription_Flag_dif_int`: description file for six-state protocol with decoy states and different intensities
- `pmSixStateLossyDescription_Flag_bit_bias`: description file for six-state protocol with decoy states and different intensities applying a bit bias
- `pmSixStateWCPChannel_Flag_dif_int`: Channel model for six-state protocol using improved decoy methods or decoy analysis from Wang et. al. [4] for different intensities
- `pmSixStateWCPChannel_Flag_bit_bias`: Channel model for six-state protocol using improved decoy methods or decoy analysis from Wang et. al. [4] for different intensities and applies a bit bias


## References
[1] D. Rusca, A. Boaron, F. Grünenfelder, A. Martin, and H. Zbinden, Finite-key analysis for the 1-decoy state QKD protocol, Applied Physics Letters 112, 171104 (2018).

[2] O. Gittsovich, N. J. Beaudry, V. Narasimhachar, R. R. Alvarez, T. Moroder, and N. Lütkenhaus, Squashing model for detectors and applications to quantum-key-distribution protocols, Physical Review A 89, 012325 (2014).

[3] N. J. Beaudry, T. Moroder, and N. Lütkenhaus, Squashing Models for Optical Measurements in Quantum Communication, Physical Review Letters 101, 093601 (2008).

[4] W. Wang and N. Lütkenhaus, Numerical security proof for the decoy-state BB84 protocol and measurement-device-independent quantum key distribution resistant against large basis misalignment, Phys. Rev. Res. 4, 043097 (2022).

