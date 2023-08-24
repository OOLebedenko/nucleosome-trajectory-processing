import numpy as np


def calc_radial_autocorr_fft(vectors, limit=-1):

    def calc_autocorr_fft(a):
        a = np.concatenate((a,np.zeros(len(a)-1))) # added zeros to your signal
        A = np.fft.fft(a)
        S = np.conj(A)*A
        c_fourier = np.fft.ifft(S)
        c_fourier = c_fourier[:(c_fourier.size//2)+1]

        return np.real(c_fourier)

    vectors = vectors.to_numpy()
    vectors_norm = np.linalg.norm(vectors, axis=1)
    n_lag_points = np.arange(vectors_norm.size, 0, -1)

    mean_acorr = calc_autocorr_fft(vectors_norm ** (-3)) / n_lag_points

    return mean_acorr[:limit]
