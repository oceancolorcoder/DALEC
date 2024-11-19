# DALEC
## dalec_L1A.m
 Read raw dalec file into arrays

 * Screen for NaNs in GPS Lat/Lon
 * Screen for QFlag "This flag for your current DALEC firmware version (v3.5-96-g130bc2f) relates
   to the GearPos movement. 0 = stationary and 1=moving during the spectrometer integration. This
   column will contain extra bits of info with the next firmware update and will let you know of the details when the time comes."
 * Screen for relAz
 * Screen for SZA
 * Screen for tilt
 * Screen for spectral outliers and NaNs

 Input: Raw text files from DALEC acquisition software (duration unclear)
 Output: L1A matlab structures
        L1A Plots

 D. Aurin NASA/GSFC November 2024
 ## dalec_L1B.m
 Process L1A DALEC to L1B
 
 * Apply calibrations and dark offsets
 Calibration file usage:
```
    K1=d0*(V-DC)+d1
    K2=e0*(V-DC)+e1
    K3=f0*(V-DC)+f1
    Ed=a0*((V-DC)/(Inttime+DeltaT_Ed)/K1)/(Tempco_Ed*(Temp-Tref)+1)
    Lu=b0*((V-DC)/(Inttime+DeltaT_Lu)/K2)/(Tempco_Lu*(Temp-Tref)+1)
    Lsky=c0*((V-DC)/(Inttime+DeltaT_Lsky)/K3)/(Tempco_Lsky*(Temp-Tref)+1)
```
 * Interpolate to common timestamps and wavebands

 Inputs: L1A files from dalec_L1A.m, calibration file from IMO
 Output: L1B files matlab structures
        L1B Plots
D. Aurin NASA/GSFC November 2024

## dalec_L2.m
 Process L1B DALEC to L2

* Add ancillary data from field log file
* Break into hourly file groups for output
* Run the spectral filter (again, now that it's one-hour files)
* Extract 300s ensembles
* Drop brightest 90 of Lt(780)
* Take the slice mean of the ensemble for Lt,Li,Es
* Calculate the Zhang et al. 2017 glint correction
* Calculate Lw and Rrs for the ensemble with uncertainty
* Calculate the NIR residual correction
   Only implemented for flat offset
   TO DO: Improve by switching between flat offset and SimSpec based on AVW
* Export L2 ensembles as hourly files (mat) and SeaBASS files

 Inputs: L1B files from dalec_L1A.m, calibration file from IMO
 Output: L12 files matlab structures
        L2 SeaBASS files
        L2 Plots of Es, Lt, Lw, Rrs with uncertainty