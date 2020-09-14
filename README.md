# SAMM
Spectral Analysis Modal Methods (SAMMs) using Non-Time-Resolved PIV

Similar to SPOD (Towne et al. 2018), SAMMs can obtain the dynamically important coherent structures from time-resolved pressure measurements and non-time-resolved PIV measurements. This is based on a linear multi-input multi-output model, as an alteration of the spectral linear stochastic estimation by Tinney et al. (2006).

SAMM includes two variants, first one is 'SAMM-SPOD' and the other one is 'SAMM-RR'. 'SAMM-SPOD' combines spectral-LSE and SPOD, while 'SAMM-RR' is based on the conditioned analysis that provides the rank-1 approximation of the system input. Both of the methods give the same dominant modes associated with a particular frequency.

The sample data given here is downsampled from the original data due to the file size limitaion. Only a portion of the PIV snapshots and two out of ten Kulite sensors are given. The test code gives the mode shapes of the v-component of a cavity flow at Mach 0.6 for the first four Rossiter mode frequencies. 
# Files
| File        | Description     | 
| ------------- |:-------------:| 
| vel_sample.mat    | velocity sample data | 
| kulite1_sample.mat   | kulite 1 sample data   |  
| kulite2_sample.mat   | kulite 2 sample data   |  
| SAMM_sample.m  | script for the sample data   |  
| SAMM_methods.m   | SAMM function   |  

# Usage
Examples:

[freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs) 
returns the POD modes for all frequencies and all snapshots using SAMM-RR method.<br/>
Y: velocity matrix, which dimension is grids x snapshots;<br/>X: pressure structure;<br/>X(i).pp is the ith sensor matrix, which dimension is snapshots times time delays;<br/>nfft: number of samples for fft calculation;<br/>fs: sampling rate

[freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs,f_index) <br/>
f_index: specified frequnecy indices (calculate all frequencies if not specified)

[freq,eigenvalues,modes] = SAMM_methods(Y,X,nfft,fs,f_index,'SPOD',blocks) <br/>
blocks: specified number of blocks (use all blocks for SPOD if not specified);<br/>
specify method 'RR' (default) or 'SPOD'

[freq,RR_eigenvalues,RR_modes,SPOD_eigenvalues,SPOD_modes] = SAMM_methods(Y,X,nfft,fs,f_index,'Both',blocks) <br/>
This retunrs POD modes for both SAMM-RR and SAMM-SPOD methods

# Reference
Yang Zhang, Louis N. Cattafesta III, Lawrence Ukeiley, *Spectral Analysis Modal Methods (SAMMs) using Non-Time-Resolved PIV*, J. of Experiments in Fluids, 2020. (accepted)


