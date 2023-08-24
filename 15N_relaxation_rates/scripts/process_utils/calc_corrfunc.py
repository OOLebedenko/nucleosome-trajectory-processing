import numpy as np


def Y2m2(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] - 1j * vectors[:, 1]) ** 2


def Y2m1(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 2.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] - 1j * vectors[:, 1]) * vectors[:, 2]


def Y20(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(5.0 / np.pi)

    return coef * (2.0 * vectors[:, 2] ** 2 - vectors[:, 0] ** 2 - vectors[:, 1] ** 2)


def Y2p1(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = -1.0 / 2.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] + 1j * vectors[:, 1]) * vectors[:, 2]


def Y2p2(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    coef = 1.0 / 4.0 * np.sqrt(15.0 / 2.0 / np.pi)

    return coef * (vectors[:, 0] + 1j * vectors[:, 1]) ** 2


def autocorr(x):
    # return real part of autocorrelation function
    f = np.fft.fft(np.pad(x, len(x), mode='constant'))
    result = np.fft.ifft(f * np.conj(f))
    result = result[:len(x)]
    result /= np.linspace(len(x), 1, len(x))
    return np.real(result)


def autocorr_all_harmonics(vectors):
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
    res = [autocorr(vectors) for vectors in [Y2m2(vectors), Y2m1(vectors), Y20(vectors)]]
    res[0] = 2.0 * res[0]
    res[1] = 2.0 * res[1]
    return 4.0 * np.pi / 5.0 * np.sum(res, axis=0)


def autocorr_different_length_arrays(arrays):
    def autocorr_all_harmonics_doesnt_averaged(vect):

        def autocorr_doesnt_averaged(x):
            # return real part of autocorrelation function
            f = np.fft.fft(np.pad(x, len(x), mode='constant'))
            result = np.fft.ifft(f * np.conj(f))
            result = result[:len(x)]
            # result /= np.linspace(len(x), 1, len(x))
            return np.real(result)

        Y_func_v1 = np.array([Y2m2(vect), Y2m1(vect), Y20(vect)])

        res = [autocorr_doesnt_averaged(f1) for f1 in Y_func_v1]
        res[0] = 2.0 * res[0]
        res[1] = 2.0 * res[1]
        return 4.0 * np.pi / 5.0 * np.sum(res, axis=0)

    # number_of_points = np.zeros(largest_size)
    delta_arr = []
    for arr in arrays:
        delta_arr.append(autocorr_all_harmonics_doesnt_averaged(arr))

    delta_arr = np.array(delta_arr, dtype=object)
    result = []
    number_of_points = []
    largest_size = np.array([arr.size for arr in delta_arr]).max()
    for delta in delta_arr:
        result.append(np.concatenate((delta, np.zeros(largest_size - delta.size))))
        number_of_point_delta = np.linspace(delta.size, 1, delta.size)
        number_of_points.append(
            np.concatenate((number_of_point_delta, np.zeros(largest_size - number_of_point_delta.size))))

    result = np.array(result, dtype=object)
    number_of_points = np.array(number_of_points, dtype=object)

    return result.sum(axis=0) / number_of_points.sum(axis=0)

