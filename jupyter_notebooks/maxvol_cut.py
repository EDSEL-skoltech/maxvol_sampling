import numpy as np
from scipy.linalg import solve_triangular, get_lapack_funcs, get_blas_funcs
from maxvolpy.maxvol import maxvol


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def f_no_cut(idx, i, copy=False):
    if copy:
        idx = np.copy(idx)
    idx[i] = 0
    return idx

def f_cut_eps(idx, i, X, eps=0.1, copy=False):
    if copy:
        idx = np.copy(idx)

    #print(np.abs(X - X[i]) < eps)
    # print(idx.shape, X.shape)
    idx[np.abs(X - X[i]) < eps] = 0
    return idx

def rect_maxvol_cut(A, tol = 1., maxK = None, min_add_K = None, minK = None, start_maxvol_iters = 10, identity_submatrix = True, top_k_index = -1, cut_fun=None, penalty=None):
    """Python implementation of rectangular 2-volume maximization. For information see :py:func:`rect_maxvol` function"""
    # tol2 - square of parameter tol
    tol2 = tol**2
    # N - number of rows, r - number of columns of matrix A
    N, r = A.shape
    if N <= r:
        return np.arange(N, dtype = int), np.eye(N, dtype = A.dtype)
    if maxK is None or maxK > N:
        maxK = N
    if maxK < r:
        maxK = r
    if minK is None or minK < r:
        minK = r
    if minK > N:
        minK = N
    if min_add_K is not None:
        minK = max(minK, r + min_add_K) 
    if minK > maxK:
        minK = maxK
    if top_k_index == -1 or top_k_index > N:
        top_k_index = N
    if top_k_index < r:
        top_k_index = r

    if cut_fun is None:
        cut_fun = f_no_cut

    if penalty is None:
        #penalty_fun = np.ones(top_k_index, dtype=int)
        chosen = np.ones(top_k_index, dtype=int)
    else:
        chosen = np.copy(penalty)


    index = np.zeros(N, dtype = int)
    tmp_index, C = maxvol(A, tol = 1, max_iters = start_maxvol_iters, top_k_index = top_k_index)
    # -- 
    index[:r] = tmp_index
    #chosen[tmp_index] = 0 -- replaced
    for ti in tmp_index:
        cut_fun(chosen, ti)
    C = np.asfortranarray(C)
    # compute square 2-norms of each row in matrix C
    row_norm_sqr = np.array([chosen[i]*np.linalg.norm(C[i], 2)**2 for i in range(top_k_index)])
    # find maximum value in row_norm_sqr
    i = np.argmax(row_norm_sqr)
    K = r
    # set cgeru or zgeru for complex numbers and dger or sger for float numbers
    try:
        ger = get_blas_funcs('geru', [C])
    except:
        ger = get_blas_funcs('ger', [C])
    while (row_norm_sqr[i] > tol2 and K < maxK) or K < minK:
        # add i to index and recompute C and square norms of each row by SVM-formula
        index[K] = i
        #chosen[i] = 0 -- replaced by the next line
        #print(chosen)
        cut_fun(chosen, i)
        if (chosen == 0).all():
            print('Failed')
        c = C[i].copy()
        v = C.dot(c.conj())
        l = 1.0/(1+v[i])
        ger(-l,v,c,a=C,overwrite_a=1)
        C = np.hstack([C, l*v.reshape(-1,1)])
        row_norm_sqr -= (l*v[:top_k_index]*v[:top_k_index].conj()).real
        row_norm_sqr *= chosen
        # find maximum value in row_norm_sqr
        i = row_norm_sqr.argmax()
        K += 1
    if identity_submatrix:
        C[index[:K]] = np.eye(K, dtype = C.dtype)
    return index[:K].copy(), C


def make_dist(X):
    n = len(X)
    A = np.empty((n, n), dtype=X.dtype)
    for ix, x in enumerate(X):
        for iy, y in enumerate(X):
            A[ix, iy] = np.abs(x - y)

    return A

def f_penal(X, bnd, level=0.0):
    Xmin = np.min(X)
    Xmax = np.max(X)
    bnd_abs = (Xmax - Xmin)*bnd
    dist = np.minimum(np.abs(X - Xmin), np.abs(Xmax - X))
    def lin_func(x):
        if bnd == 0:
            return x*0.0 + 1.0  # crookedly, but it works. Ann, never do like this!
        else:
            return (1.0 - level)*np.minimum(x, bnd_abs)/bnd_abs + level

    return lin_func(dist)


def f_penal_2D(X, Y, bnd, level=0.0):
    return f_penal(X, bnd=bnd, level=level)*f_penal(Y, bnd=bnd, level=level)


if __name__ == '__main__':
    import matplotlib.pyplot as plt


    def get_x_mat(n, m):
        A = np.empty((n, m), dtype=float)
        x = np.linspace(-1, 1, n, endpoint=True)
        A[:, 0] = 1.0
        for i in range(1, m):
            A[:, i] = A[:, i-1]*x
        
        return A, x


    A, x = get_x_mat(10000, 11)

    expr_type = 2

    if expr_type == 1:
        f_cut = lambda i, idx : f_cut_eps(i, idx, X=x, eps=0.05)
        num_row = 20
        num_col = A.shape[1]
        idx, _ = rect_maxvol_cut(A, minK=num_row, maxK=num_row, cut_fun=f_cut)
        plt.plot(x[idx[:num_col]], np.ones(len(idx[:num_col])), 'b.');
        plt.plot(x[idx[num_col:]], np.ones(len(idx[num_col:])), 'r.');
        plt.show()

    if expr_type == 2:

        num_row = 20
        num_col = A.shape[1]
        bnd = 0.3
        idx, _ = rect_maxvol_cut(A, minK=num_row, maxK=num_row, penalty=f_penal(x, bnd=bnd))
        plt.plot(x[idx[:num_col]], np.ones(len(idx[:num_col])), 'b.');
        plt.plot(x[idx[num_col:]], np.ones(len(idx[num_col:])), 'r.');
        plt.show()

