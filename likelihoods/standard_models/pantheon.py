import redshift_remappingLZ as redshift_remapping
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import pandas as pd

def your_likelihood_function(omegam, w, wa, Ha0, rdragT:
    omegal = 1-omegam
    a=0
    b=0
    c=0
    d=0
    d1=0
    d2=0
    if not (a <= b <= c <= d):
        return -np.inf
    else:
        model = redshift_remapping.Redshift_Remapping(omegam, omegal, w, wa, Ha0, rdragT, a, b, c, d, d1, d2)

        # Load data
        file_path = 'pantheon_sn.txt'
        df = pd.read_csv(file_path, sep=' ')
        df = df[df['zHD'] >= 0.01]
        mag = df['MU_SH0ES'].values
        zobsmag = df['zHD'].values
        indices_to_keep = df.index

        # Load and process covariance matrix
        covariance_matrix = np.loadtxt('Pantheon+SH0ES_STAT+SYS.cov')
        matrix_dimension = 1701
        CM = covariance_matrix.reshape((matrix_dimension, matrix_dimension))
        CM = CM[np.ix_(indices_to_keep, indices_to_keep)]
        CM_inv = np.linalg.inv(CM)

        # Predicted magnitudes from the model
        mag_pred = model.magnitude(zobsmag)

        # Marginalize over Msn
        deriv = np.ones_like(mag)[:, None]  # Derivative with respect to Msn is just ones
        deriv_prime = CM_inv @ deriv  # Project into the inverse covariance space
        fisher = deriv.T @ deriv_prime  # Compute the Fisher matrix

        # Update inverse covariance matrix to marginalize over Msn
        CM_inv_marginalized = CM_inv - deriv_prime @ np.linalg.inv(fisher) @ deriv_prime.T

        # Compute the chi-squared with marginalized covariance matrix
        residuals = mag - mag_pred
        error_mag = residuals.T @ CM_inv_marginalized @ residuals

        etotal = error_mag
        return -etotal / 2  # Return the log-likelihood
def get_H0rdrag(H0, rdrag):
    return H0*rdrag
    
