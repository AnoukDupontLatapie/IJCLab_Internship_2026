# IJCLab_Internship_2026
CMB spectral distortions internship at IJCLab for 7 weeks.

The internship was divided into two parts : 

- First 4 weeks : theoretical approach (what characterize a spectral distortion, when and why do they happen, how can we define them mathematicaly...)
- Second 4 weeks : foreground modelling. The goal was to quantify the small fluctuations between model and theory so no mixing could occur between real spectral distortions and fluctuations of the model. The models used were solely dust and synchrotron. 


Supervisors : Thibaut LOUIS and Léo VACHER

Informations for each file :

plot_SD.py : - plot of $\mu$ and $y$ distortion spectrum along with a blackbody and a shifted blackbody.
             - animated plot of the distortion spectrum from $2\ 10^6$ to $10^2$ redshift.

func_theo.py : analytical expressions of blackbody, modified blackbody (Pysm-defined), temperature shift, $\mu$ and $y$ distortion

bruit_fossil.py : frequential binning method for FOSSIL sensibility.

conversion.py : conversions necessary to switch beetween units

mask_file_download.py : downloading of the mask file

dust_SED.py : Masked intensity of dust generated with PySM (execution : ~10 min)

fit_dust.py : Modified blackbody fit + residues

sync_SED.py : Masked intensity of synchrotron emission generated with PySM (execution : ~10 min)

fit_sync.py : Power law fit + residues