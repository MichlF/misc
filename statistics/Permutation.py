import cv2
import numpy as np
from scipy.stats import ttest_rel, ttest_ind


def clusterbased_permutation(X1, X2, p_val=0.05, cl_p_val=0.05, paired=True, tail='both', nr_perm=1000, mask=None, conn=None):
    '''
    Implements Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG- data. 
    Journal of Neurosience Methods, 164(1), 177?190. http://doi.org/10.1016/J.Jneumeth.2007.03.024

    Arguments
    - - - - - 

    X1 (array): subject X dim1 X dim2 (optional), where dim1 and dim2 can be any type of dimension 
                            (time, frequency, electrode, etc). Values in array represent some dependent
                            measure (e.g classification accuracy or power)
    X2 (array | float): either a datamatrix with same dimensions as X1, or a single value 
                            against which X1 will be tested
    p_val (float): p_value used for inclusion into the cluster
    cl_p_val (float): p_value for evaluation overall cluster significance
    paired (bool): paired t testing (True) or independent t testing (False)
    tail (str): apply one- or two- tailed t testing
    nr_perm (int): number of permutations
    mask (array): dim1 X dim2 array. Can be used to restrict cluster based test to a specific region. 
    conn (array): outlines which dim1 points are connected to other dim1 points. Usefull
                              when doing a cluster based permutation test across electrodes 

    Returns
    - - - -

    cl_p_vals (array): dim1 X dim2 with p-values < cl_p_val for significant clusters and 1's for all other clusters

    '''

    # if no mask is provided include all datapoints in analysis
    if mask == None:
        mask = np.array(np.ones(X1.shape[1:]), dtype=bool)
        print(
            f'\nUsing all {mask.size} datapoints in cluster based permutation')
    elif mask.shape != X1[0].shape:
        print('\nMask does not have the same shape as X1. Adjust mask!')
    else:
        print(
            f'\nThere are {int(mask.sum())} out of {mask.size} datapoints in your mask during cluster based permutation')

    # check whether X2 is a chance variable or a data array
    if isinstance(X2, (float, int)):
        X2 = np.tile(X2, X1.shape)

    # compute observed cluster statistics
    pos_sizes, neg_sizes, pos_labels, neg_labels, sig_cl = compute_clustersizes(
        X1, X2, p_val, paired, tail, mask, conn)
    cl_p_vals = np.ones(sig_cl.shape)

    # iterate to determine how often permuted clusters exceed the observed cluster threshold
    c_pos_cl = np.zeros(np.max(np.unique(pos_labels)))
    c_neg_cl = np.zeros(np.max(np.unique(neg_labels)))

    # initiate random arrays
    X1_rand = np.zeros(X1.shape)
    X2_rand = np.zeros(X1.shape)

    for p in range(nr_perm):

        print(f"{(float(p)/nr_perm)*100}% of permutations\r")

        # create random partitions
        if paired:  # keep observations paired under permutation
            rand_idx = np.random.rand(X1.shape[0]) < 0.5
            X1_rand[rand_idx, :] = X1[rand_idx, :]
            X1_rand[~rand_idx, :] = X2[~rand_idx, :]
            X2_rand[rand_idx, :] = X2[rand_idx, :]
            X2_rand[~rand_idx, :] = X1[~rand_idx, :]
        else:  # fully randomize observations under permutation
            all_X = np.vstack((X1, X2))
            all_X = all_X[np.random.permutation(all_X.shape[0]), :]
            X1_rand = all_X[:X1.shape[0], :]
            X2_rand = all_X[X1.shape[0]:, :]

        # compute cluster statistics under random permutation
        rand_pos_sizes, rand_neg_sizes, _, _, _ = compute_clustersizes(
            X1_rand, X2_rand, p_val, paired, tail, mask, conn)
        max_rand = np.max(np.hstack((rand_pos_sizes, rand_neg_sizes)))

        # count cluster p values
        c_pos_cl += max_rand > pos_sizes
        c_neg_cl += max_rand > neg_sizes

    # compute cluster p values
    p_pos = c_pos_cl / nr_perm
    p_neg = c_neg_cl / nr_perm

    # remove clusters that do not pass threshold
    if tail == 'both':
        # 0 is not a cluster
        for i, cl in enumerate(np.unique(pos_labels)[1:]):
            if p_pos[i] < cl_p_val/2:
                cl_p_vals[pos_labels == cl] = p_pos[i]
            else:
                pos_labels[pos_labels == cl] = 0

        # 0 is not a cluster
        for i, cl in enumerate(np.unique(neg_labels)[1:]):
            if p_neg[i] < cl_p_val/2:
                cl_p_vals[neg_labels == cl] = p_neg[i]
            else:
                neg_labels[neg_labels == cl] = 0

    elif tail == 'right':
        # 0 is not a cluster
        for i, cl in enumerate(np.unique(pos_labels)[1:]):
            if p_pos[i] < cl_p_val:
                cl_p_vals[pos_labels == cl] = p_pos[i]
            else:
                pos_labels[pos_labels == cl] = 0

    elif tail == 'left':
        # 0 is not a cluster
        for i, cl in enumerate(np.unique(neg_labels)[1:]):
            if p_neg[i] < cl_p_val:
                cl_p_vals[neg_labels == cl] = p_neg[i]
            else:
                neg_labels[neg_labels == cl] = 0

    # ADD FUNCTION TO GET

    return cl_p_vals


def compute_clustersizes(X1, X2, p_val, paired, tail, mask, conn):
    '''

    Helper function for clusterbased_permutation (see documentation)

    # ! NOTE:
    At the moment only supports two tailed tests
    At the moment does not support connectivity
    '''

    # STEP 1: determine 'actual' p value
    # apply the mask to restrict the data
    X1_mask = X1[:, mask]
    X2_mask = X2[:, mask]

    p_vals = np.ones(mask.shape)
    t_vals = np.zeros(mask.shape)

    if paired:
        t_vals[mask], p_vals[mask] = ttest_rel(X1_mask, X2_mask)
    else:
        t_vals[mask], p_vals[mask] = ttest_ind(X1_mask, X2_mask)

    # initialize clusters and use mask to restrict relevant info
    sign_cl = np.mean(X1, 0) - np.mean(X2, 0)
    sign_cl[~mask] = 0
    p_vals[~mask] = 1

    # STEP 2: apply threshold and determine positive and negative clusters
    cl_mask = p_vals < p_val
    pos_cl = np.zeros(cl_mask.shape)
    neg_cl = np.zeros(cl_mask.shape)
    pos_cl[sign_cl > 0] = cl_mask[sign_cl > 0]
    neg_cl[sign_cl < 0] = cl_mask[sign_cl < 0]

    # STEP 3: label clusters
    if conn == None:
        nr_p, pos_labels = cv2.connectedComponents(np.uint8(pos_cl))
        nr_n, neg_labels = cv2.connectedComponents(np.uint8(neg_cl))
        # hack to control for onedimensional data (CHECK whether correct)
        pos_labels = np.squeeze(pos_labels)
        neg_labels = np.squeeze(neg_labels)
    else:
        print('Function does not yet support connectivity')

    # STEP 4: compute the sum of t stats in each cluster (pos and neg)
    pos_sizes, neg_sizes = np.zeros(nr_p - 1), np.zeros(nr_n - 1)
    for i, label in enumerate(np.unique(pos_labels)[1:]):
        pos_sizes[i] = np.sum(t_vals[pos_labels == label])

    for i, label in enumerate(np.unique(neg_labels)[1:]):
        neg_sizes[i] = abs(np.sum(t_vals[neg_labels == label]))

    if sum(pos_sizes) == 0:
        pos_sizes = 0

    if sum(neg_sizes) == 0:
        neg_sizes = 0

    return pos_sizes, neg_sizes, pos_labels, neg_labels, p_vals


def cluster_plot(X1, X2, times, y, p_val=.05, cl_p_val=.05, color='black', ls='-', linewidth=2):

    sig_cl = clusterbased_permutation(X1, X2, p_val=p_val, cl_p_val=cl_p_val)
    mask = np.where(sig_cl < 1)[0]
    sig_cl = np.split(mask, np.where(np.diff(mask) != 1)[0]+1)
    for cl in sig_cl:
        plt.plot(times[cl], np.ones(cl.size) * y,
                 color=color, ls=ls, linewidth=linewidth)
