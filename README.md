# Dynamics Fluidity Toolbox

This toolbox is made to compute the instantaneous dynamic fluidity of multivariate time series in MATLAB.
This measure tracks temporal clustering of system’s configurations in dynamical phase space without prior assumptions about the system’s structure.
Originally developed for atmospheric sciences (Faranda et al., 2017), it was recently applied to multichannel EEG recordings to study early alterations in Alzheimer's disease (Aguilera et al., 2025).

## Usage example: EEG pipeline

This example shows how to apply the toolbox to EEG recordings, following the procedure from Aguilera et al., 2025.

### Data
Here we provide an example of 100 seconds of pre-processed 30 channels EEG recording with a sampling frequency of 1000Hz.

```
% Loading the EEG data
load('EEG_Data.mat'); % Variable should be named "Data", with format Channels x Time

% Visualizing the EEG data
figure
for ch = 1:30
plot(Data(ch,:)-ch);
hold on;
end
```

### Coarsegraining the EEG
This step can be optionnal, however, fluidity computation time can rise exponentially with the number of point. In the Aguilera et al., 2025 study, the EEG signal was coarsegrained with a 40ms time window (i.e. the signal was averaged over 40ms time windows) to reduce the size of the signal and fit to the 25 fps frame rate of the video recording.
The coaresegraining can be performed with the *coarsegrain* function as follow

```
% Coarsegraining the data with a 40ms non-overlapping window

Coarse_EEG = coarsegrain(Data, 40, 'unit', 'ms');
```

### Computing Fluidity
The *dynamics fluidity* is then computed on the coarsegrained EEG time series with the *dynamics_fluidity* function to obtain one single vector of fluidity.

```
% Computing fluidity of the Coarsegrained EEG

Fluidity = dynamics_fluidity(Coarse_EEG);
```

### Extracting Results
From the fluidity time series, global measures can be extracted as the mean or the median fluidity.
```
% Mean and Median fluidity

Mean_Fluidity = mean(Fluidity);
Median_Fluidity = median(Fluidity);
```
However, for a better sensitivity, one could observe the distributions of fludity values using *histogram* MATLAB function and defined bins.
Here we will use default standard bin vector but should be adapted depending on the values, however, for group comparison, bins should be the same for everyone.

```
% Set the bins for distribution

% Create a vector of 500 bins between 0.5 and 1 as for this data set the minimal fluidity is around 0.8
Bins = linspace(0.5, 1, 500); 

% Compute distribution
h = histogram(Fluidity, Bins);
Distribution = h.Values

% Plot Distribution
figure;
plot(Bins(2:end), Distribution);
xlabel('Fluidity');
ylabel('Counts');
```
Then statistical tools can be used for comparing distributions

## References
Faranda, D., Messori, G., & Yiou, P. (2017). Dynamical proxies of North Atlantic predictability and extremes. Scientific Reports, 7, 41278. 
https://doi.org/10.1038/srep41278

Aguilera, M., Mathis, C., Herbeaux, K., Isik, A., Faranda, D., Battaglia, D., Goutagny, R. (2025). 40 Hz light stimulation restores early brain dynamics alterations and associative memory in Alzheimer’s disease model mice. 
Imaging Neuroscience, 3, IMAG.a.70. https://doi.org/10.1162/IMAG.a.70
