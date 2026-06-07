# Direction of Arrival (DOA), Large Aperture Array (LAA) and Multiple-Input Multiple-Output (MIMO) Radar MATLAB Simulations
## Project Overview
### With the growing popularity of array-based localisation, different localisation frameworks use different signal models, array configurations and signal processing approaches, resulting in differences in localisation performance in different environments. Each technique may be better suited to particular operating environments and applications. Thus, this project focuses on conducting a comparative analysis of three array-based localisation techniques: Direction of Arrival (DOA) localisation, Large Aperture Array (LAA) localisation, and MIMO Radar. The simulations on this Github repository can be easily run locally. 

## Usage instructions
### DOA Localisation Simulation
`doa_single.m`: Performs single-target DOA Localisation using MUSIC-based DOA estimation
| Input Parameter | Description |
|-----------|-------------|
| `SNR_dB` | Signal-to-noise ratio in dB |
| `target` | 2x1 vector of true target position given by `[x;y]` |
| `n_elements` | Number of array elements in the Uniform Circular Array (UCA) at each receiver |
| `N_snapshots` | Number of snapshots |
| `r` | Receiver coordinates, where each column is the position of one receiver |

### LAA Localisation Simulation
`laa_single.m`: Performs single-target LAA Localisation using spherical wave model
| Parameter | Description |
|-----------|-------------|
| `SNR_dB` | Signal-to-noise ratio in dB |
| `target` | 2x1 vector of true target position given by `[x;y]` |
| `r` | Receiver coordinates, where each column is the position of one receiver |
| `N` | Number of snapshots |
| `alpha` | Path-loss exponent (In free space, `alpha` = 2) |

### MIMO Radar Localisation Simulation
`main_MIMO.m`: Performs MIMO Radar localisation of one or more targets 
| Parameter | Description |
|-----------|-------------|
| `num_PRIs` | Number of Pulse Repetition Intervals |
| `SNR_dB` | Signal-to-noise ratio in dB |
| `target_range` | True target range(s), specified as a vector |
| `target_angles_deg` | True target azimuth angle(s), must contain same number of elements as `target_range` |
| `target_velocity` | True target velocity(ies) |
| `virtual` | `true` uses Virtual MIMO mode, `false` uses Classical MIMO mode |
| `N_tx` | Number of transmit antennas |
| `N_rx` | Number of receive antennas |
