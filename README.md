# Conformal Models Repository 🌌

This repository contains all the files related with the project **Exploring Redshift-Scale Factor Relations in Cosmological Models with DESI**. Check these two papers for more information:
-  [Wojtak, R. & Prada, F. 2016, MNRAS, 458, 3331](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3331W/abstract)
- [Wojtak, R. & Prada, F. 2017, MNRAS, 470, 4493](https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4493W/abstract)


---

## Structure of the respository 📁

### **1. conformal_model_equations/** 📚
Here you can find the files with the equations of the distances conformally transformed for three different cases:

- `redshift_remapping.py`: file that contains all the conformally transformed distances. This can be used when we are able to constrain rdrag.
- `redshift_remapping_pantheon.py`: file that contains all the conformally transformed distances. This can be used for SNIa, where rdrag and Ha0 cannot be constrained (we need Cepehids or CMB as callibrators)
- `redshift_remapping_DESI.py`: file that contains all the conformally transformed distances. This can be used for BAO distance,s where we can constrain rdrag*Ha0, but not rdrag or Ha0 independently.


### **2. data/** 📊
Here you can find all the data used in this proyect:

- `DES-SN5YR_HD.csv/`: Hubble diagram from DES-SN-Y5 CITA (can also be downloaded [here](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT)) 
-  `STAT+SYS.txt.gz`: covariance matrix DES (download it [here](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT)) 
- `pantheon_sn.txt/`: Hubble diagram from Pantheon+ (can also be downloaded [here](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)).
- `Pantheon+SH0ES_STAT+SYS.cov/`: covariance matrix Pantheon+ (can also be downloaded [here](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)).

### **3. likelihoods/** 📔
Here you can find the likelihoods defined for each of the data mentioned above:

- `conformal_models/DESI.py`: file that calculates the likelihood from DESI with conformal transformation.
- `conformal_models/DES.py`: file that calculates the likelihood from DES with conformal transformation.
- `conformal_models/pantheon.py`: file that calculates the likelihood from Pantheon+ with conformal transformation.
- `standard_models/DESI.py`: file that calculates the likelihood from DESI using the standard FLRW metric. 
- `standard_models/DES.py`: file that calculates the likelihood from DES using the standard FLRW metric. 
- `standard_models/pantheon.py`: file that calculates the likelihood from Pantheon+ using the standard FLRW metric. 

### **4. yaml_files/** 📈
Here you can find the configuration files (.yaml) necessary in order to run the MCMC chains using Cobaya for different combinations of the datasets:

- `conformal_models/DESI_Planck_Pantheon.yaml`: configuration file needed to compute the MCMC chain using DESI, Planck 2018 and Pantheon+ with conformal transformation in the LCDM model.
- `conformal_models/DESI_Planck_Pantheon_CDM.yaml`: configuration file needed to compute the MCMC chain using DESI, Planck 2018 and Pantheon+ with conformal transformation in the CDM model.
- `conformal_models/DESI_planck_DES.yaml`: configuration file needed to compute the MCMC chain using DESI, Planck 2018 and DES with conformal transformation in the LCDM model.
- `conformal_models/full_w.yaml`: configuration file needed to compute the MCMC chain using DESI, Planck 2018 and Pantheon+ with conformal transformation in the wCDM model.
- `conformal_models/full_w0wa.yaml`: configuration file needed to compute the MCMC chain using DESI, Planck 2018 and Pantheon+ with conformal transformation in the w0waCDM model.
- `conformal_models/pantheon.yaml`: configuration file needed to compute the MCMC chain using Pantheon+ with conformal transformation in the LCDM model.
- `conformal_models/planck.yaml`: configuration file needed to compute the MCMC chain Planck 2018 with conformal transformation in the LCDM model.
- `conformal_models/planck_CDM.yaml`: configuration file needed to compute the MCMC chain Planck 2018 with conformal transformation in the CDM model.
- `standard_models/DESI_Planck_pantheon.yaml`: configuration file needed to compute the MCMC chain DESI, Planck and Pantheon+  with conformal transformation in the LCDM model.
- `standard_models/DESI_pantheon.yaml`: configuration file needed to compute the MCMC chain DESI and Pantheon+  with conformal transformation in the LCDM model.
- `standard_models/pantheon.yaml`: configuration file needed to compute the MCMC chain Pantheon+ with conformal transformation in the LCDM model.

---

## Installation 🧑‍🔬

1. Clone the repository:
   ```bash
   git clone https://github.com/elenaferga/Conformal-Models.git

2. Installing the dependencies:
   ```bash
   pip install -r requirements.txt


# Contributions 🤝

If you would like to contrubute, please open an **issue** or send a **pull request**. All the contributions are welcome. You can also send me an e-mail (see below).

# Contact ✉️

For any consult, please contact me at **efdez@iaa.es**
