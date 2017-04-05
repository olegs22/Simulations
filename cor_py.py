import numpy as np


def point2c(data, b_size, bins):
    n_data = len(data[:,0])

    r_x = np.random.rand(n_data)*b_size
    r_y = np.random.rand(n_data)*b_size
    r_z = np.random.rand(n_data)*b_size
    dat_r = np.array([r_x,r_y,r_z]).T

    max_d = np.sqrt(3.)*b_size
    bin_width = max_d / bins
    spacing = np.linspace(0,max_d, bins)
    hist0 = np.zeros(bins)
    hist1 = np.zeros(bins)
    hist2 = np.zeros(bins)
    corr = np.zeros(bins)

    for i in range(n_data):
        for j in range(i,n_data):
            DD_s = pow(data[j,:]-data[i,:], 2)
            RR_s = pow(dat_r[j,:]-dat_r[i,:], 2)

            DD = np.sqrt(np.sum(DD_s, axis=0))
            RR = np.sqrt(np.sum(RR_s, axis=0))

            nbin_d = int(round(DD / bin_width))
            nbin_r = int(round(RR / bin_width))

            hist0[nbin_d] += 1.0
            hist1[nbin_r] += 1.0

        for k in range(n_data):
            DR_s = pow(data[i,:]-dat_r[k,:], 2)

            DR = np.sqrt(np.sum(DR_s, axis=0))

            n_bin_rd = int(round(DR / bin_width))

            hist2[n_bin_rd] += 1.0

    for i in range(bins):
        if hist1[i] == 0.0:
            corr[i] = 0.0
    corr = (hist0 - hist2 + hist1) / hist1

    return corr
