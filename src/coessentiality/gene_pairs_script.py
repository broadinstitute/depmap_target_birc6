#!/usr/bin/env python
# coding: utf-8


import numpy as np, statsmodels.api as sm
from load_screens import load_screens
from tqdm import tqdm
import sys

# Load batch-corrected screens

screens = load_screens()
genes = screens.index

# Compute p-values and effect directions with GLS



def calculate_pairs(start_index, end_index, chuck_num, screens = screens):

    end_index = min(end_index, len(screens.index))
    subset_screens = screens.iloc[start_index:end_index,:]

    cholsigmainv = np.linalg.cholesky(np.linalg.pinv(screens.cov())).T
    GLS_beta = np.empty((len(subset_screens), len(screens)))
    GLS_p = np.empty((len(subset_screens), len(screens)))
    for A_index, A_row in tqdm(enumerate(subset_screens.values)):
        input = cholsigmainv.dot(sm.add_constant(A_row))
        for B_index, B_row in enumerate(screens.values):
            output = cholsigmainv.dot(B_row)
            model = sm.OLS(output, input).fit()
            GLS_beta[A_index, B_index] = model.params[1]
            GLS_p[A_index, B_index] = model.pvalues[1]
    #np.fill_diagonal(GLS_p, 1)
    GLS_sign = np.sign(GLS_beta)

    # Save everything

    np.save('GLS_p_{}.npy'.format(chuck_num), GLS_p)
    np.save('GLS_sign_{}.npy'.format(chuck_num), GLS_sign)
    subset_screens.index.to_series().to_csv('genes_{}.txt'.format(chuck_num), index=False, header=False)

print(sys.argv)
calculate_pairs(start_index = int(sys.argv[1]), end_index= int(sys.argv[2]), chuck_num = int(sys.argv[3]))
