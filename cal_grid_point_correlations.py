import numpy as np
from scipy.stats import linregress

def cal_grid_point_correlations(X,Y,return_pvals=False):
    """ Calculate correlations between two fields at each grid-point. Assumes
    fields have time as first dimension then two spatial dimensions."""

    # check that time series and field have the same time dimension
    for i in np.arange(X.ndim):
        if X.shape[i] != Y.shape[i]:
            print("Fields 1 and 2 do not have the same dimensions")
            print("Field 1 size:",X.shape,"Field 2 size:",Y.shape)
            raise ValueError

    corrs_map = np.zeros_like(X[0,:,:])
    pvals_map = np.zeros_like(X[0,:,:])

    # calculate anomalies
    X_anom = X - np.mean(X,axis=0)
    Y_anom = Y - np.mean(Y,axis=0)

    print("Calculating grid-point correlations...")
    milestone = 20
    for i in np.arange(X.shape[1]):
        if 100*float(i+1)/float(X.shape[1]) >= milestone:
            print("%d %% complete" % (milestone))
            milestone = milestone + 20

        for j in np.arange(X.shape[2]):
            try:
                slope, intercept, r_value, p_value, std_err = linregress(X_anom[:,i,j],Y_anom[:,i,j])
            except ValueError:
                #print(X_anom[:,i,j])
                slope, intercept, r_value, p_value, std_err = 0, 0, 0, 0, 0
            corrs_map[i,j] = r_value
            pvals_map[i,j] = p_value

    if return_pvals == True: return corrs_map, pvals_map
    else: return corrs_map
