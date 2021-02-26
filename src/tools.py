import numpy as np
# from maxvol_cut import rect_maxvol_cut, f_no_cut, f_penal_2D
# from sklearn.naive_bayes import GaussianNB

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
    print(idx.shape, X.shape)
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




def norm_data(X, bounds=(-1.0, 1.0), copy=True):
    X = np.array(X, copy=copy).T
    for i, x in enumerate(X):
        # print(len(x))
        min_v, max_v = np.min(x), np.max(x)
        b = (bounds[0]*max_v - bounds[1]*min_v)/(max_v-min_v)
        k = float(bounds[1] - bounds[0])/(max_v-min_v)
        X[i] *= k
        X[i] += b
        

    return X.T

def points_selection(X, max_n_pnts, min_n_pnts, cut_fun=None, penalty = None):
    
    """Function for selecting optimal parameters for dimentionality reduction method and for clustering.
    
    Parameters 
    ----------------
    X: array with shape (number_of_pixels*number_of_features)
            Initial data
           
    """
    #MaxVol
    
    res = rect_maxvol_cut(X, maxK=max_n_pnts, minK=min_n_pnts, cut_fun=cut_fun, penalty=penalty)[0]

    return res


def add_coords(X=None, size=(285, 217), order='C', idx_good_mask=None):
    """
    order can by 'C' or 'F'
    """
    w, h = size
    x_coord, y_coord = np.meshgrid(np.arange(h), np.arange(w))
    
    
    if idx_good_mask is None:
        idx_good_mask = np.arange(x_coord.size)
    
    if X is None:
        return np.hstack((
            x_coord.flatten(order=order)[idx_good_mask, None],
            y_coord.flatten(order=order)[idx_good_mask, None]))
    else:
        return np.hstack((np.array(X, copy=False),
                          x_coord.flatten(order=order)[idx_good_mask, None],
                          y_coord.flatten(order=order)[idx_good_mask, None]))
    
def gen_input(mode, data, shapes,mask):
    modes = ['usual', 'normed',
         'XY', 'XY_normed']
    fn_X_embedded = modes[mode]
    return [
        lambda x: np.array(x),
        lambda x: norm_data(x),
        lambda x: add_coords(
            x, size=shapes[0], idx_good_mask=mask),
        lambda x: norm_data(gen_input(2, x, shapes, mask)[0], copy=False),
    ][mode](data), fn_X_embedded


def my_score(a, b):
    a = np.array(a, copy=False)
    b = np.array(b, copy=False)
    n = len(a)
    assert len(b) == n, 'Arrays of differnet shapes :((('
    m = len(a[a==b])
    return float(m)/float(n)


def extend_score(S_tru, S_obt):
    vals, cnts = np.unique(S_tru, return_counts=True)

    part_scores = np.array([np.count_nonzero(S_obt[S_tru==v]==v) for v in vals], dtype=float)
    part_scores /= cnts
    wts = np.ones_like(vals)/len(vals)
    return np.average(part_scores, weights=wts)
    

def points_selection_accuracy(X,y, n_pnts, cut_fun=None, penalty = None, to_ret_pred=False):
    
    """Function for selecting optimal parameters for dimentionality reduction method and for clustering.
    
    Parameters 
    ----------------
    X: array with shape (number_of_pixels*number_of_features)
            Initial data
    
    data_mode: 0,1,2,3
            Mode for preparing data.
            0 - raw data
            1 - normed data
            2 - unnormed data with coordinates as features
            3 - normed data with coordinates as feature           
 
    """
   
    #MaxVol
    num_pnts = n_pnts
    res = rect_maxvol_cut(X, maxK=num_pnts, minK=num_pnts, cut_fun=cut_fun, penalty=penalty)[0]

    #Naive Bayes
    gnb = GaussianNB()
    gnb_model = gnb.fit(X[res], y[res])
    scores = gnb_model.score(X, y)
    if to_ret_pred:
       
        my_sc = extend_score(y, gnb_model.predict(X))
        print("Accuracy: %0.2f" % (my_sc))
        return my_sc, res
    
    else:
        print("Accuracy: %0.2f" % (scores.mean()))
        return scores, res

# enforce df_drop as X
def points_selection_accuracy_tail(data_mode, cut_fun=None, penalty = None):
    return points_selection_accuracy(df_5features, data_mode, cut_fun=cut_fun, penalty = f_penal)

def f_no_cut(idx, i, copy=False):
    if copy:
        idx = np.copy(idx)
    idx[i] = 0
    return idx

def f_cut_eps(idx, i, X, eps=0.1, copy=False):
    if copy:
        idx = np.copy(idx)
    xx = X[:, -2] 
    yy = X[:, -1]   
    #idx[i] = 0
    idx[(xx - xx[i])**2 + (yy-yy[i])**2 <= eps**2] = 0
    return idx


def calc_score(idx, X, y, to_ret_pred=True):
    gnb = GaussianNB()
    gnb_model = gnb.fit(X[idx], y[idx])
    
    if to_ret_pred:
        scores = extend_score(y, gnb_model.predict(X))
    else:
        scores = gnb_model.score(X, y)

    return scores


def good_points_brute_force(idx, num, X, y):
    sc = -1
    cmb_good = None
    for comb in combinations(idx, num):
        comb = np.array(comb)
        #print(comb)
        sc_curr = calc_score(comb, X=X, y=y, to_ret_pred=True)
        if sc_curr > sc:
            sc = sc_curr
            cmb_good = comb
            
    return cmb_good, sc

def idx_to_idx(idx_big, idx):
    hass = dict()
    for i, elem in enumerate(idx_big):
        hass[elem] = i
        
    return np.array([hass[i] for i in idx])
