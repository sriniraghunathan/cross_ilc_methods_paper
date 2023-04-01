#For fisher forecasting
import numpy as np, sys, os, glob
from scipy.io import readsav
from scipy import interpolate as intrp
import scipy as sc
from scipy import integrate, stats

def get_cl_ksz(els, Aksz = 1., dl_ksz_amp_total = 3.):
    dl_fac = els * (els+1)/2/np.pi
    dl_ksz = np.tile(dl_ksz_amp_total, len(els))
    cl_ksz = dl_ksz / dl_fac
    cl_ksz[np.isnan(cl_ksz)] = 0.
    cl_ksz[np.isinf(cl_ksz)] = 0.
    return Aksz * cl_ksz

def get_homo_ksz(els, Aksz_h = 1., alphaksz_h = 0., el_norm = 3000):
    dl_fac = els * (els+1)/2/np.pi
    dl_ksz_homo = Aksz_h * (els/el_norm)**alphaksz_h
    cl_ksz_homo = dl_ksz_homo/dl_fac

    return cl_ksz_homo

def get_knox_errors_parent(els, cl_dic, nl11_dic, fsky, nl22_dic = None, nl12_dic = None):

    delta_cl_dic = {}
    for XX in cl_dic:
        if XX == 'TT':
            nl11 = nl11_dic['TT']
            if nl22_dic is not None:
                nl22 = nl22_dic['TT']
                nl12 = nl12_dic['TT']
        elif XX == 'EE' or XX == 'BB':
            nl11 = nl11_dic['EE']
            if nl22_dic is not None:
                nl22 = nl22_dic['EE']
                nl12 = nl12_dic['EE']
        elif XX == 'TE':
            nl11 = np.copy(nl11_dic['TT']) * 0.
            if nl22_dic is not None:
                nl22 = np.copy(nl11) * 0.
                nl12 = np.copy(nl11) * 0.

        cl11 = cl_dic[XX] + nl11
        if nl22_dic is not None:
            cl22 = cl_dic[XX] + nl22
            cl12 = cl_dic[XX] + nl12
        else:
            cl22, cl12 = None, None

        delta_cl_dic[XX] = get_knox_errors(els, cl11, fsky, cl22 = cl22, cl12 = cl12)

    return delta_cl_dic

def get_knox_errors(els, cl11, fsky, cl22 = None, cl12 = None):

    """
    get Knox bandpower errors.

    els: multipoles
    cl11: signal + noise in the first map.
    cl22: signal + noise in the second map.
    cl12: signal + noise cross spectrum of the two maps.
    """

    delta_el = np.diff(els)[0]

    if cl22 is not None and cl12 is not None:
        cl_total = np.sqrt( (cl12**2. + cl11 * cl22)/2. )
    else:
        cl_total = cl11

    cl_knox_err = np.sqrt(2./ (2.*els + 1.) / fsky / delta_el ) * (cl_total)
    cl_knox_err[np.isnan(cl_knox_err)] = 0.
    cl_knox_err[np.isinf(cl_knox_err)] = 0.

    return cl_knox_err


def fix_params(F_mat, param_names, fix_params):

    #remove parameters that must be fixed    
    F_mat_refined = []
    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 in fix_params or p2 in fix_params: continue
            F_mat_refined.append( (F_mat[pcntr2, pcntr1]) )

    totparamsafterfixing = int( np.sqrt( len(F_mat_refined) ) )
    F_mat_refined = np.asarray( F_mat_refined ).reshape( (totparamsafterfixing, totparamsafterfixing) )

    param_names_refined = []
    for p in param_names:
        if p in fix_params: continue
        param_names_refined.append(p)


    return F_mat_refined, param_names_refined

def add_priors(F_mat, param_names, prior_dic):

    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 == p2 and p1 in prior_dic:
                prior_val = prior_dic[p1]
                F_mat[pcntr2, pcntr1] += 1./prior_val**2.

    return F_mat

def get_fisher_matrix(els, cl_deriv_dic, delta_cl_dic, params, pspectra_to_use, min_l_temp = None, max_l_temp = None, min_l_pol = None, max_l_pol = None):

    if min_l_temp is None: min_l_temp = 0
    if max_l_temp is None: max_l_temp = 10000

    if min_l_pol is None: min_l_pol = 0
    if max_l_pol is None: max_l_pol = 10000

    npar = len(params)
    F = np.zeros([npar,npar])
    #els = np.arange( len( delta_cl_dic.values()[0] ) )

    with_lensing = 0
    if 'PP' in pspectra_to_use:
        with_lensing = 1

    all_pspectra_to_use = []
    for tmp in pspectra_to_use:
        if isinstance(tmp, list):      
            all_pspectra_to_use.extend(tmp)
        else:
            all_pspectra_to_use.append(tmp)

    for lcntr, l in enumerate( els ):

        TT, EE, TE = 0., 0., 0.
        Tphi = Ephi = PP = 0.
        if 'TT' in delta_cl_dic:
            TT = delta_cl_dic['TT'][lcntr]
        if 'EE' in delta_cl_dic:
            EE = delta_cl_dic['EE'][lcntr]
        if 'TE' in delta_cl_dic:
            TE = delta_cl_dic['TE'][lcntr]
        if with_lensing:
            Tphi, Ephi, PP = delta_cl_dic['Tphi'][lcntr], delta_cl_dic['Ephi'][lcntr], delta_cl_dic['PP'][lcntr]

        null_TT, null_EE, null_TE = 0, 0, 0
        if l<min_l_temp or l>max_l_temp:
            null_TT = 1
        if l<min_l_pol or l>max_l_pol: 
            null_EE = 1
            null_TE = 1
        null_PP = 0 #Lensing noise curves already have pretty large noise outside desired L range
        #if l<min_l_TE or l>max_l_TE:  
        #    null_TE = 1

        #20200611
        if 'TT' not in all_pspectra_to_use:
            null_TT = 1
        if 'EE' not in all_pspectra_to_use:
            null_EE = 1
        if 'TE' not in all_pspectra_to_use:
            #if 'TT' not in pspectra_to_use and 'EE' not in pspectra_to_use:
            #    null_TE = 1
            if 'TT' in pspectra_to_use and 'EE' in pspectra_to_use:
                null_TE = 0
            else:
                null_TE = 1
        if ['TT', 'EE', 'TE'] in pspectra_to_use:
            null_TT = 0
            null_EE = 0
            null_TE = 0
        #20200611

        ##if (null_TT and null_EE and null_TE): continue# and null_PP): continue
        if (null_TT and null_EE and null_TE and null_PP): continue

        param_combinations = []
        for pcnt,p in enumerate(params):
            for pcnt2,p2 in enumerate(params):
                ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
                param_combinations.append([p,p2, pcnt, pcnt2])

        def get_cov(TT, EE, TE, PP, TP, EP):

            C = np.zeros( (3,3) ) #TT, EE, PP
            C[0,0] = TT
            C[1,1] = EE
            C[0,1] = C[1,0] = TE

            C[2,2] = PP
            C[0,2] = C[2,0] = TP
            C[1,2] = C[2,1] = EP ##0. ##EP

            return np.mat( C )

        #nulling unwanted fields
        if null_TT and null_TE: TT = 0
        if null_EE and null_TE: EE = 0
        #if null_TE and (null_TT and null_EE): TE = 0
        if null_TE: 
            if not null_TT and not null_EE:
                pass
            else:
                TE = 0
        if null_PP: PP = Tphi = EPhi = 0
        if null_TT: Tphi = 0
        if null_EE: Ephi = 0
        #nulling unwanted fields

        COV_mat_l = get_cov(TT, EE, TE, PP, Tphi, Ephi)
        inv_COV_mat_l = sc.linalg.pinv2(COV_mat_l)
        #inv_COV_mat_l = np.linalg.inv(COV_mat_l)

        if (0):##l%500 == 0: 
            from IPython import embed; embed()
            print(l, null_TT, null_EE, null_TE, null_PP)
            print(COV_mat_l)

        for (p,p2, pcnt, pcnt2) in param_combinations:


            TT_der1, EE_der1, TE_der1 = 0., 0., 0.
            TT_der2, EE_der2, TE_der2 = 0., 0., 0.

            if 'TT' in cl_deriv_dic[p]:
                TT_der1 = cl_deriv_dic[p]['TT'][lcntr]
                TT_der2 = cl_deriv_dic[p2]['TT'][lcntr]
            if 'EE' in cl_deriv_dic[p]:
                EE_der1 = cl_deriv_dic[p]['EE'][lcntr]
                EE_der2 = cl_deriv_dic[p2]['EE'][lcntr]
            if 'TE' in cl_deriv_dic[p]:
                TE_der1 = cl_deriv_dic[p]['TE'][lcntr]
                TE_der2 = cl_deriv_dic[p2]['TE'][lcntr]


            if with_lensing:
                PP_der1, TPhi_der1, EPhi_der1 = cl_deriv_dic[p]['PP'][lcntr], cl_deriv_dic[p]['Tphi'][lcntr], cl_deriv_dic[p]['Ephi'][lcntr]
                PP_der2, TPhi_der2, EPhi_der2 = cl_deriv_dic[p2]['PP'][lcntr], cl_deriv_dic[p2]['Tphi'][lcntr], cl_deriv_dic[p2]['Ephi'][lcntr]
            else:
                PP_der1 = PP_der2 = 0.
                TPhi_der1 = TPhi_der2 = 0. 
                EPhi_der1 = EPhi_der2 = 0.


            if null_TT: TT_der1 = TT_der2 = TPhi_der1 = TPhi_der2 = 0
            if null_EE: EE_der1 = EE_der2 = EPhi_der1 = EPhi_der2 = 0
            if null_TE: TE_der1 = TE_der2 = 0
            if null_PP: PP_der1 = PP_der2 = 0

            fprime1_l_vec = get_cov(TT_der1, EE_der1, TE_der1, PP_der1, TPhi_der1, EPhi_der1)
            fprime2_l_vec = get_cov(TT_der2, EE_der2, TE_der2, PP_der2, TPhi_der2, EPhi_der2)

            curr_val = np.trace( np.dot( np.dot(inv_COV_mat_l, fprime1_l_vec), np.dot(inv_COV_mat_l, fprime2_l_vec) ) )

            F[pcnt2,pcnt] += curr_val

    return F   

def get_fisher_inv(F_mat):

    F_mat = np.asarray(F_mat)

    Flen = len(F_mat)
    all_inds = np.arange(Flen)

    F_mat_diag = np.diag(F_mat)
    good_inds = np.where(F_mat_diag > 0)[0]
    Flen_refined = len(good_inds)

    #from IPython import embed; embed(); sys.exit()

    F_mat_refined = []
    used_i, used_j = [], []
    for i in all_inds:
        for j in all_inds:
            if i in good_inds and j in good_inds: 
                F_mat_refined.append( (F_mat[j, i]) )
                used_i.append(i)
                used_j.append(j)
    used_i = np.asarray(used_i)
    used_j = np.asarray(used_j)
    F_mat_refined = np.asarray( F_mat_refined ).reshape( (Flen_refined, Flen_refined) )
    C_mat_refined = np.linalg.inv(F_mat_refined)

    C_mat = np.zeros(F_mat.shape)
    C_mat[used_j, used_i] = C_mat_refined.flatten()

    return C_mat
