import numpy as np
from maxvol_cut import rect_maxvol_cut, f_no_cut, f_penal_2D
from sklearn.naive_bayes import GaussianNB

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