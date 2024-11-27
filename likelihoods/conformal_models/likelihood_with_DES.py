import redshift_remappingLZ as redshift_remapping
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import pandas as pd

def your_likelihood_function(omegam, w0, wa, Ha0, rdragT, a, b, c, d, d1, d2):
    omegal = 1-omegam
    if not (a <= b <= c <= d):
        return -np.inf
    else:
        model = redshift_remapping.Redshift_Remapping(omegam, omegal, w0, wa, Ha0, rdragT, a, b, c, d, d1, d2)

        file_path = 'DES-SN5YR_HD.csv'
        df = pd.read_csv(file_path, sep=',')
        df = df[df['zHD']>=0.01]
        mag =  df['MU'].values
        zobsmag =  df['zHD'].values
        indices_to_keep = df.index


        covariance_matrix = np.loadtxt('covsys_000.txt')
    
    
        matrix_dimension = 1829
        CM = covariance_matrix.reshape((matrix_dimension, matrix_dimension))
        CM = CM[np.ix_(indices_to_keep, indices_to_keep)]
        CM = np.linalg.inv(CM)

        mag_pred = model.magnitude(zobsmag)

        error_mag = 0
        rows, cols = np.shape(CM)
        for i in range(rows):
            for j in range(cols):
                aux = (mag[i]-mag_pred[i])*(mag[j]-mag_pred[j])*CM[i][j]
                error_mag = error_mag + aux


        etotal=0
        etotal = etotal + error_mag 
        return -etotal/2
    
def get_H0rdrag(H0, rdrag):
    return H0*rdrag
    