###########################################################################
###########################################################################
###########################################################################
import os, copy, pickle
from pylab import *
if (1):
    from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
    rcParams['font.family'] = 'serif'
    rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')
    #from matplotlib import transforms

import tools
import scipy as sc, warnings, argparse
np.seterr(divide='ignore', invalid='ignore')
import warnings, matplotlib.cbook
warnings.filterwarnings('ignore', category=UserWarning)

from matplotlib.patches import Ellipse

#simple line fitting case

########################################################
########################################################
#parameters
afid, bfid = 1., 0.
astep, bstep = 0.01, 0.01
params = ['a', 'b']
fix_params = []
use_parameter_shifts_for_sys_errors = 0
if use_parameter_shifts_for_sys_errors:
    asys = -0.05 #-0.05 #-0.05 #-0.05
    bsys = 0. #-.2
else:
    std_error_frac_for_sys = 1.
    asys = bsys = None

which_param_to_constrain = params
for p in params:
    if p in which_param_to_constrain: continue
    fix_params.append(p)
param_dict = {}
param_dict['a'] = [afid, 0.3] #fid, delta for axis limits
param_dict['b'] = [bfid, 0.3] 
########################################################
########################################################
reclen = 5#2 #00 #100
x = np.arange(reclen)#-reclen/2, reclen/2)
if (1):
    x = np.asarray( [-1, 1, 0] ) #np.arange(reclen)
    x = np.asarray( [-1, 1, 1] ) #np.arange(reclen)
    x = np.asarray( [-2, -1, 0, 1, 2, 3] ) #np.arange(reclen)
    reclen = len(x)
std_error = 0.1 #sqrt(variance) in each observed poin
errors = np.zeros(reclen) + std_error
########################################################
########################################################
#get derivatives now
der_dic = {}
der_dic['a'] = np.zeros(reclen) + x
der_dic['b'] = np.zeros(reclen) + 1.
'''
for p in params:
    if p == 'a':
        delta = astep
        a1, a2 = afid - astep, afid + astep
        b1, b2 = bfid, bfid
    elif p == 'b':
        delta = bstep
        a1, a2 = afid, afid
        b1, b2 = bfid - bstep, bfid + bstep
    der_dic[p] = ( get_model(x, a2, b2) - get_model(x, a1, b1) ) / (2 * delta)
'''
if (0):
    for p in params:
        plot(der_dic[p], label = r'%s' %(p))
    legend(loc = 1)
    show()
########################################################
########################################################
#Fisher approach
param_combinations = []
for pcnt,p in enumerate(params):
    for pcnt2,p2 in enumerate(params):
        ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
        param_combinations.append([p,p2, pcnt, pcnt2])

nparams = len(params)
F_mat = np.zeros((nparams, nparams))
for i in range(reclen):
    curr_sigma = errors[i]
    for (p1,p2, pcnt1, pcnt2) in param_combinations:
        der1_vec = np.asarray( [der_dic[p1][i]] )
        der2_vec = np.asarray( [der_dic[p2][i]] )
        inv_cov_mat = curr_sigma**-2.
        curr_val = np.dot(der1_vec, np.dot( inv_cov_mat, der2_vec ))
        F_mat[pcnt2,pcnt1] += curr_val

F_mat, params = tools.fix_params(F_mat, params, fix_params)
fisher_C_mat = np.linalg.inv(F_mat)
print(fisher_C_mat)
#print(np.diag(fisher_C_mat)**0.5); #sys.exit()

#plot Fisher constraints
tr = tc = nparams
F_dic = {}
F_dic['Simple'] = F_mat
#make_triangle_plot(F_dic, param_dict, tr, tc, params)

########################################################
########################################################
#calculate bias vector
#Eq.(8) of https://arxiv.org/pdf/0710.5171.pdf
if use_parameter_shifts_for_sys_errors:
    ysys = tools.get_data(x, asys, bsys)
else:
    ysys = np.ones(reclen) * std_error * std_error_frac_for_sys
bias_vector = []
for pcntr, p in enumerate( params ):
    curr_bias_vector = []
    for i in range(reclen):
        der1_vec = np.asarray( [ysys[i]] )
        der2_vec = np.asarray( [der_dic[p][i]] )
        inv_cov_mat = curr_sigma**-2.
        curr_bias_val = np.dot(der1_vec, np.dot( inv_cov_mat, der2_vec ))
        curr_bias_vector.append(curr_bias_val)
    bias_vector.append( np.sum(curr_bias_vector) )
bias_vector = np.asarray(bias_vector)
final_bias_vector = np.asarray( np.dot( np.mat(bias_vector), fisher_C_mat ) )[0]
final_bias_dic = {}
for p, b in zip(params, final_bias_vector):
    final_bias_dic[p] = b
print(final_bias_vector)
########################################################
########################################################
#MLE approach
debug = 0 ##1
totsims = 10 #10 #100 #25
noofsims_for_cov = int(reclen**2 * 100)
anal_cov = np.eye(reclen) * std_error**2.

if (1):
    #generate covariance
    sim_arr = []
    for n in range(noofsims_for_cov):
        errors = np.random.standard_normal(reclen) * std_error
        ysim = tools.get_data(x, afid, bfid, err = errors ) #generate sim
        sim_arr.append( ysim )
    sim_mean = np.mean(sim_arr, axis = 0)
    sim_arr = np.asarray(sim_arr) - sim_mean
    sim_arr = np.mat(sim_arr).T    
    cov = (sim_arr * sim_arr.T) / (noofsims_for_cov - 1) #compute pixel-pixel covarince between different sims: cov should have dimensions (reclen x reclen)
    if debug:
        subplot(121); imshow(anal_cov); colorbar(); title(r'Analytical')
        subplot(122); imshow(cov); colorbar(); title(r'Using sims')
        show(); sys.exit()
else: #analytic covariance
    cov = anal_cov

if len(which_param_to_constrain)>1:
    mina, maxa, dela = afid-.2, afid+.2, 0.005
    minb, maxb, delb = bfid-.3, bfid+.3, 0.005
    a_arr = np.arange(mina, maxa, dela)
    b_arr = np.arange(minb, maxb, delb)
    alen, blen = len(a_arr), len(b_arr)
else:
    if which_param_to_constrain[0] == 'a':
        minfit, maxfit, delfit = afid-.2, afid+.2, 0.001
        bval = bfid
    elif which_param_to_constrain[0] == 'b':
        minfit, maxfit, delfit = bfid-.2, bfid+.2, 0.001
        aval = afid
    fit_arr = np.arange(minfit, maxfit, delfit)


master_L_arr, master_logL_arr = [], []
for t in range(totsims):
    print(t)
    ##generate mock data first
    ytrue = tools.get_data(x, afid, bfid)
    errors = np.random.standard_normal(reclen) * std_error
    yobs = tools.get_data(x, afid, bfid, err = errors )
    yobs += ysys
    if debug:
        plot(x, yobs, 'o', label = r'Data'); plot(x, ytrue, label = r'True');
        legend(loc = 1); show(); sys.exit()

    #perform likelihood fit for slope a
    if len(which_param_to_constrain) == 2:
        logL_arr = []
        for bval in b_arr:
            for aval in a_arr:
                ymodel = tools.get_data(x, aval, bval)
                logLval = tools.get_likelihood(yobs, ymodel, cov)
                logL_arr.append( logLval )
        logL_arr = np.asarray(logL_arr).reshape(blen, alen)
    else:
        logL_arr = []
        for fval in fit_arr:
            if which_param_to_constrain[0] == 'a':
                aval = fval
            elif which_param_to_constrain[0] == 'b':
                bval = fval
            ymodel = tools.get_data(x, aval, bval)
            logLval = tools.get_likelihood(yobs, ymodel, cov)
            logL_arr.append( logLval )

    master_logL_arr.append( np.copy(logL_arr) )

    #logL to likelihoods
    L_arr = np.exp( logL_arr - np.max(logL_arr) )
    L_arr /= np.max(L_arr)

    if (0):
        imshow(L_arr, extent = [mina, maxa, minb, maxb]); colorbar(); show(); sys.exit()

    master_L_arr.append(L_arr)

master_L_arr_copy = np.copy(master_L_arr)
master_logL_arr_copy = np.copy(master_logL_arr)

master_logL_arr = np.sum(master_logL_arr, axis = 0) #logL combined from all sims

master_logL_arr_copy = master_logL_arr_copy.tolist()
master_logL_arr_copy.append(master_logL_arr)
master_logL_arr_copy = np.asarray(master_logL_arr_copy)

master_L_arr = np.prod(master_L_arr, axis = 0)
master_L_arr = master_L_arr/np.max(master_L_arr) #Likelihoods combined from all sims
#master_recov_fit, master_recov_fit_low_err, master_recov_fit_high_err = np.asarray( tools.get_width_from_sampling(fit_arr, L_arr) )
#master_recov_fit_width = (master_recov_fit_high_err + master_recov_fit_low_err)/2.
if (0):
    if len(which_param_to_constrain)>1:
        X, Y = np.meshgrid(a_arr, b_arr)
        master_logL_arr -= np.max(master_logL_arr)
        imshow(master_logL_arr, extent = [mina, maxa, minb, maxb], origin = 'lower', aspect = 'auto'); colorbar(); 
        #pcolor(X, Y, master_logL_arr); colorbar(); 
        sigarr = [-1.]#, 5.]
        for lev in sigarr:
            contour(X, Y, master_logL_arr, levels = lev, colors = ['k'])#, extent = [xmin_A, xmax_A, xmin_alpha, xmax_alpha])#, alpha = 0.4
        xlabel(r'$a$'); ylabel(r'$b$'); show(); sys.exit()

########################################################
########################################################
#show plots
fsval = 12

#sys.exit()
if len(params) > 1:

    import copy

    #Fisher corner plot
    expname = 'Fisher with bias'
    tr = tc = nparams
    F_dic = {}
    F_dic[expname] = F_mat
    ls_dic = {}
    ls_dic[expname] = '-.'
    color_dic = {}
    color_dic[expname] = 'k'
    alpha_dic = {}
    alpha_dic[expname] = 0.2

    param_dict_copy = copy.deepcopy(param_dict)

    #show likelihoods now
    X, Y = np.meshgrid(a_arr, b_arr)
    for lcntr, logL_arr in enumerate( master_logL_arr_copy ):


        #logL for a and b
        logL_arr = logL_arr - np.max(logL_arr)        
        for pntr, p in enumerate( params ):

            if lcntr == 0 and pcntr == 0:
                labval = r'MLE (1 sim)'
                labval_fisher = r'Fisher'                
            else:
                labval = None
                labval_fisher = None

            if p == 'a':
                axval = 0
                fit_arr = a_arr
                pstr = r'$a$'
                sbpl = 1
                minx, maxx = mina, maxa
                fidval = afid
            elif p == 'b':
                axval = 1
                fit_arr = b_arr
                pstr = r'$b$'
                sbpl = 4
                minx, maxx = minb*2, maxb*2
                fidval = bfid

            #logL_arr = logL_arr - np.max(logL_arr)
            #logL_arr = logL_arr / np.mean(logL_arr)

            logL_arr_p = np.mean(np.copy(logL_arr), axis = axval)
            L_arr_p = np.exp(logL_arr_p - np.max(logL_arr_p)); 
            L_arr_p/=np.max(L_arr_p)

            recov_fit, recov_fit_low_err, recov_fit_high_err = np.asarray( tools.get_width_from_sampling(fit_arr, L_arr_p) )
            recov_fit_width = (recov_fit_high_err + recov_fit_low_err)/2.

            param_dict_copy[p][0] = recov_fit

            ax=subplot(tr, tc, sbpl)
            if lcntr == len(master_logL_arr_copy) - 1:
                colroval = 'darkred'
                lwval = 2.
                zorderval = 1000.
                alphaval = 1.
            else:
                colroval = 'goldenrod'
                lwval = 0.5
                zorderval = -1000.
                alphaval = 0.5
            plot(fit_arr, L_arr_p, color = colroval, lw = lwval, label = labval, zorder = zorderval, alpha = alphaval)#, lw = 0.25, alpha = 0.25); 
            if p == 'b':
                xlabel(pstr, fontsize = fsval)
            axvline(fidval, lw = 0.5)
            #xlim(minx, maxx)
            x, deltax = param_dict_copy[p]
            x1, x2 = x - deltax/2., x + deltax/2.
            print(x1, x2, p)
            xlim(x1, x2)
            setp(ax.get_yticklabels(), visible=False); tick_params(axis='y', which='both', left=False, right = False)

        #2D contour
        subplot(tr, tc, 3)
        sigarr = [-1.]#, 5.]
        for lev in sigarr:
            #imshow(logL_arr, extent = [mina, maxa, minb, maxb], origin = 'lower', aspect = 'auto'); colorbar(); 
            contour(X, Y, logL_arr, levels = [lev], colors = [colroval], lw = lwval, linestyle = 'solid', zorder = zorderval, alpha = alphaval)
        axvline(afid, lw = 0.5); axhline(bfid, lw = 0.5)
        xlabel(r'$a$', fontsize = fsval)
        ylabel(r'$b$', fontsize = fsval)
        #xlim(mina, maxa); ylim(minb, maxb)

        x, deltax = param_dict_copy['a']
        y, deltay = param_dict_copy['b']

        x1, x2 = x - deltax/2., x + deltax/2.
        y1, y2 = y - deltay/2., y + deltay/2.

        xlim(x1, x2); ylim(y1, y2)

        if (0):
            #param_dict_copy[p][0] = recov_fit        
            #recov_fit_inds = np.unravel_index(np.argmax(logL_arr), logL_arr.shape)
            #param_dict_copy['a'] = [a_arr[recov_fit_inds[1]], maxa-mina]
            #param_dict_copy['b'] = [b_arr[recov_fit_inds[0]], maxb-minb]
            param_dict_copy['a'][1] = (maxa-mina)/2.
            param_dict_copy['b'][1] = (maxb-minb)/2.


        if (0):#lcntr < len(master_logL_arr_copy) - 1:
            tools.make_triangle_plot(F_dic, param_dict_copy, tr, tc, params, color_dic = color_dic, ls_dic = ls_dic, alpha_dic = alpha_dic, lwval = 0.5, show_one_sigma_lab = 0)

        #xlabel(r'$a$'); ylabel(r'$b$'); 

    #show final bias lines
    if (asys ==0. and bsys == 0.) and (asys is None and bsys is None):
        titstr = r'$y=ax + b$ ($a = %g, b = %g, \sigma_{\rm stat} = %g, \sigma_{\rm sys} = {\rm None}$)' %(afid, bfid, std_error)
    else:
        if use_parameter_shifts_for_sys_errors:
            axvline(sysval, label = r'$%s_{\rm true} + %s_{\rm sys}$' %(pstr, pstr), color = 'darkgreen', ls = '-', lw = 0.5)
            titstr = r'$y=ax + b$ ($a = %g, b = %g, \sigma_{\rm stat} = %g, a_{\rm sys} = %g, b_{\rm sys} = %g$)' %(afid, bfid, std_error, asys, bsys)
        else:
            titstr = r'$y=ax + b$ ($a = %g, b = %g, \sigma_{\rm stat} = %g, \sigma_{\rm sys} = %g\sigma_{\rm stat}$)' %(afid, bfid, std_error, std_error_frac_for_sys)
        for p in final_bias_dic:
            if p == 'a':
                subplot(tr, tc, 1)
                axvline(afid+final_bias_dic[p], label = r'Fisher bias prediction', color = 'royalblue', ls = ':', lw = 2., zorder = -1000)
                subplot(tr, tc, 3)
                axvline(afid+final_bias_dic[p], label = r'Fisher bias prediction', color = 'royalblue', ls = ':', lw = 2., zorder = -1000)
            elif p == 'b':
                subplot(tr, tc, 4)
                axvline(bfid+final_bias_dic[p], label = r'Fisher bias prediction', color = 'royalblue', ls = ':', lw = 2., zorder = -1000)
                subplot(tr, tc, 3)
                axhline(bfid+final_bias_dic[p], label = r'Fisher bias prediction', color = 'royalblue', ls = ':', lw = 2., zorder = -1000)

    suptitle(titstr, fontsize = fsval + 2)
    if use_parameter_shifts_for_sys_errors:
        plname = 'fisher_2params_bias_calc_linefitting_%.2f%ssys.png' %(asys, which_param_to_constrain)
    else:
        plname = 'fisher_2params_bias_calc_linefitting_syserroris%.2fsigma.png' %(std_error_frac_for_sys)
    print(plname)
    savefig(plname, dpi = 200.)

    show(); sys.exit()

else:
    #Fisher
    fisher_widthval = fisher_C_mat[0, 0]**0.5
    x, deltax = param_dict['a']
    x1, x2 = x-deltax, x+deltax
    hor_fisher, ver_fisher = tools.get_Gaussian(x, fisher_widthval, x1, x2)

    titstr = r'$y=ax + b$ ($a = %g, b = %g, a_{\rm sys} = %g, b_{\rm sys} = %g$)' %(afid, bfid, asys, bsys)
    title(titstr, fontsize = fsval + 2)
    #MLE
    #plot(fit_arr, L_arr, label = r'MLE', color = 'orangered'); 
    for lcntr, L_arr in enumerate( master_L_arr_copy ):
        if lcntr == 0:
            labval = r'MLE (1 sim)'
            labval_fisher = r'Fisher'
        else:
            labval = None
            labval_fisher = None
        plot(fit_arr, L_arr, color = 'goldenrod', lw = 0.5, label = labval)#, lw = 0.25, alpha = 0.25); 

        #get best-fit value and width
        recov_fit, recov_fit_low_err, recov_fit_high_err = np.asarray( tools.get_width_from_sampling(fit_arr, L_arr) )
        recov_fit_width = (recov_fit_high_err + recov_fit_low_err)/2.

        xshift = recov_fit - x
        plot(hor_fisher + xshift, ver_fisher, label = labval_fisher, color = 'black', lw = 0.2, ls = '-.')
    #hor_mle, ver_mle = tools.get_Gaussian(x, recov_fit_width, x1, x2)
    #hor_mle+=recov_fit
    #plot(hor_mle, ver_mle, label = r'MLE', color = 'orangered'); 


    #combined likelihood from many sims
    plot(fit_arr, master_L_arr, label = r'MLE (%s sims: combined)' %(totsims), color = 'darkred', lw = 1.5); 

    if which_param_to_constrain[0] == 'a':
        fidval = afid
        sysval = asys
        pstr = 'a'
    elif which_param_to_constrain[0] == 'b':
        fidval = bfid
        sysval = bsys
        pstr = 'b'

    axvline(fidval, label = r'$%s_{\rm true}$' %(pstr), color = 'gray', ls = '-', alpha = 0.5); 
    #axvline(master_recov_fit, color = 'darkred'); 
    axvline(sysval, label = r'$%s_{\rm true} + %s_{\rm sys}$' %(pstr, pstr), color = 'darkgreen', ls = '-', lw = 0.5)
    #axvline(fidval+final_bias_val[0], label = r'Fisher bias prediction', color = 'royalblue', ls = ':', lw = 2.)

    xlabel(r'$%s$' %(pstr), fontsize = fsval)
    ylabel(r'Normalised $\mathcal{L}$', fontsize = fsval)
    legend(loc = 1, fancybox = True, fontsize = fsval-2)
    xlim( minfit, maxfit )

    show()
