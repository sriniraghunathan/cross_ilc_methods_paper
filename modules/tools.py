import numpy as np, sys, os, glob
from scipy.io import readsav
from scipy import interpolate as intrp
import scipy as sc
from scipy import integrate, stats
import ilc


h, k, c, temp=6.62607e-34, 1.38065e-23, 2.9979e8, 2.725

################################################################
################################################################
def Jy_Sr_K(band, Jy_K = True, K_Jy = False):
    """
    Conversion from Jy/Sr to Kelvin and vice versa.

    Parameters
    ----------
    band: float
        freqeuncy in GHz.
    Jy_K: bool
        Default is True. Convert from flux to temperature.
    k_Jy: bool
        Default is False. Convert from temperature to flux.

    Returns
    -------
    conv: float
        conversion factor from Jy/Sr or vice versa.
    """
    band *= 1e9

    x=h*band/(k*temp)
    
    dbnu = 2.*k*(band**2/c**2)*(x**2*np.exp(x))/(np.exp(x)-1.)**2.0
    conv=1e26*dbnu #Jy to K

    if Jy_K:
        return 1./conv
    else:
        return conv

def interpolate_dn_ds(s, dnds, increasing_spacing_by = 10):
    """
    interpolate dN/ds to increase the resolution.
    This is particularly important for high flux (S150>=2 mJy) sources.

    Parameters
    ----------
    s: array
        flux bins.
    dnds: array
        source count in each flux bin.
    increasing_spacing_by: int
        interploation factor. Default is x10.
    
    Returns
    -------
    s_ip: array
        interpolated flux bins.
    dnds_ip: array
        interpolated source counts in each flux bin.
    """

    lns = np.log(s)
    dlns =  (lns[1]-lns[0])

    #get a new closely separated range using interpolation
    dlns_ip = dlns / increasing_spacing_by
    #lns_ip = np.logspace(min(lns), max(lns), dlns_ip)
    s_ip = np.exp( np.arange(min(lns), max(lns), dlns_ip) )
    dnds_ip = np.interp(s_ip, s, dnds)
        
    return s_ip, dnds_ip

def get_poisson_source_counts(band, min_flux_mJy = -1., max_flux_mJy = 6.4e-3, band0 = 150., which_dnds = 'lagache', spec_index_radio = -0.76, dnds_ip_factor = 10):

    """
    get radio source population dN/ds for the desired band.

    Parameters
    ----------
    band: float
        freqeuncy in GHz.
    min_flux_mJy: float
        min. flux required in mJy.
        default is -1 --> no minimum threshold.
    max_flux_mJy: float
        min. flux required in mJy.
    band0: float
        default band in which the source count file is defined.
        default is 150 GHz. Does not work for others.
    which_dnds: str
        dn/ds model to be used.
        default in lagache.
        options are lagache and dezotti.
    spec_index_radio: float
        radio spectral index to be used to scale from band0 to other bands.
        default is -0.7 (R21 SPT result).
    dnds_ip_factor: int
        interploation factor for dN/ds. Default is x10.
    
    Returns
    -------
    s: array
        flux array.
    nsources: array
        total source in each flux bin.
    dndlns: array
        logrithminc source counts in each flux bin.
    """
    if which_dnds == 'dezotti':
        assert band0 == 150
        de_zotti_number_counts_file = 'publish/data/counts_150GHz_de_Zotti_radio.res'
        radio_log_s,radio_x,radio_y,radio_logs25n = np.loadtxt(de_zotti_number_counts_file, skiprows = 12, unpack = True)

        #first get the number of sources in each flux range
        s = 10**radio_log_s

        s150_flux = np.copy(s)
        ###print(s); sys.exit()
        #perform scaling
        s = (band/band0)**spec_index_radio * s
        s25n = 10**radio_logs25n
        dnds = s25n / s**2.5

        #20221101 - perform interpolation
        s_ip, dnds_ip = interpolate_dn_ds(s, dnds, increasing_spacing_by = dnds_ip_factor)
        s, dnds = s_ip, dnds_ip
        s150_flux = np.copy(s)
        s = (band/band0)**spec_index_radio * s

        ###print(s[-50:]); sys.exit()

        #get dn/ds and number of sources
        lns = np.log(s)
        dlns = lns[1] - lns[0]
        dndlns = dnds * s
        nsources = dndlns * dlns  # number of sources obtained

    elif which_dnds == 'lagache':
        assert band0 == 150
        lagache_number_counts_file = 'publish/data/lagache_2019_ns150_radio.dat'
        s, dnds = np.loadtxt(lagache_number_counts_file, unpack = True)

        #20221101 - perform interpolation
        s_ip, dnds_ip = interpolate_dn_ds(s, dnds, increasing_spacing_by = dnds_ip_factor)
        s, dnds = s_ip, dnds_ip

        s150_flux = np.copy(s)
        #perform scaling
        s = (band/band0)**spec_index_radio * s
        lns = np.log(s)
        dlns =  (lns[1]-lns[0])
        dndlns = dnds*s
        nsources = dndlns * dlns  # number of sources obtained

    elif which_dnds == 'agora':
        assert band0 == 150
        mdpl2_number_counts_searchstr = 'publish/data/mdpl2_radio/mdpl2_radiocat_len_universemachine_trinity_95_150_220ghz_randflux_datta2018_truncgauss_*'
        ##print(mdpl2_number_counts_searchstr); sys.exit()
        mdpl2_number_counts_flist = sorted( glob.glob( mdpl2_number_counts_searchstr ) )
        band_ind_dic = {90: 3, 95: 3, 150: 6, 220: 9} #total LK + HK sources
        colind = band_ind_dic[band0]
        s_arr, dnds_s2p5_arr = [], []
        for mdpl2_number_counts_fname in mdpl2_number_counts_flist:
            s, dnds_s2p5 = np.loadtxt(mdpl2_number_counts_fname, usecols = [0, colind], unpack = True) #Flux S is in Jy.
            s_arr.append(s)
            dnds_s2p5_arr.append(dnds_s2p5)
        s_arr = np.asarray( s_arr )
        dnds_s2p5_arr = np.asarray( dnds_s2p5_arr )
        s, dnds_s2p5 = np.mean(s_arr, axis = 0), np.mean(dnds_s2p5_arr, axis = 0)
        dnds = dnds_s2p5/s**2.5

        s_ip, dnds_ip = interpolate_dn_ds(s, dnds, increasing_spacing_by = dnds_ip_factor)
        s, dnds = s_ip, dnds_ip
        s150_flux = np.copy(s)
        s = (band/band0)**spec_index_radio * s

        #get dn/ds and number of sources
        lns = np.log(s)
        dlns = lns[1] - lns[0]
        dndlns = dnds * s
        nsources = dndlns * dlns  # number of sources obtained

    #pick the required flux indices
    sinds = np.where((s150_flux>min_flux_mJy) & (s150_flux<=max_flux_mJy))
    s, nsources, dndlns = s[sinds], nsources[sinds], dndlns[sinds]

    return s, nsources, dndlns

def get_poisson_power_with_spectra_index_scatter(band1, band2 = None, min_flux_mJy_band0 = -1., max_flux_mJy_band0 = 6e-3, band0 = 150., which_dnds = 'dezotti', spec_index_radio = -0.76, spec_index_radio_scatter = 0.6, units = 'uk', dnds_ip_factor = 1, max_radio_band = 230.):

    flux_arr_band0, counts_arr_band0, dndlns_arr_band0 = get_poisson_source_counts(band0, min_flux_mJy = min_flux_mJy_band0, max_flux_mJy = max_flux_mJy_band0, band0 = band0, which_dnds = which_dnds, spec_index_radio = spec_index_radio, dnds_ip_factor = dnds_ip_factor)
    flux_arr, counts_arr, dndlns_arr = flux_arr_band0, counts_arr_band0, dndlns_arr_band0

    if band2 is None: band2 = band1

    s_arr = flux_arr
    #print(s_arr); sys.exit()
    lns = np.log(s_arr)
    dlns =  (lns[1]-lns[0])

    ##print(dlns, s_arr, min_flux_mJy_band0, max_flux_mJy_band0, s_arr**2. * dndlns_arr, 'hihihi'); sys.exit()
    if spec_index_radio_scatter == 0.:
        min_alpha, max_alpha = spec_index_radio - 5. * 0.6, spec_index_radio + 5. * 0.6
        alpha_arr_for_int = np.arange(min_alpha, max_alpha, 0.01)
    else:
        min_alpha, max_alpha = spec_index_radio - 5. * spec_index_radio_scatter, spec_index_radio + 5. * spec_index_radio_scatter
        alpha_arr_for_int = np.arange(min_alpha, max_alpha, 0.01)

    int_arr = []
    for s, dndlns in zip(s_arr, dndlns_arr):
        def integrand_flux_dalpha(s, dndlns, dlns, alpha, alpha_mean, alpha_sigma):                
            return s**2.*(dndlns) * dlns * (band1 * band2 / band0**2.)**alpha * stats.norm.pdf(alpha, alpha_mean, alpha_sigma)
            
        int_val = integrate.simps( integrand_flux_dalpha(s, dndlns, dlns, alpha_arr_for_int, spec_index_radio, spec_index_radio_scatter), x=alpha_arr_for_int )
        ##int_val = s**2.*(dndlns) * dlns
        int_arr.append( int_val )

    integral = np.sum(int_arr)

    Jy_to_K_conv = Jy_Sr_K(band1) * Jy_Sr_K(band2)
    #print(Jy_Sr_K(band1)*1e6, Jy_Sr_K(band2)*1e6); sys.exit()
    cl_poisson = integral * Jy_to_K_conv

    if units == 'uk': cl_poisson *= 1e12

    if band1 > max_radio_band or band2 > max_radio_band: cl_poisson *= 0.

    return cl_poisson

def get_radio_ilc_residuals(els, freqarr, wl1, wl2, which_dnds_arr, min_flux_mJy=0.1e-3, max_flux_mJy=4e-3, spec_index_radio = -0.76, spec_index_radio_scatter_arr = [0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01], dnds_ip_factor = 10, quiet = False):
    radio_cl_dic = {}
    for which_dnds in which_dnds_arr:
        if not quiet: print('\tdN/ds = %s' %(which_dnds))
        radio_cl_dic[which_dnds] = {}
        for spec_index_radio_scatter in spec_index_radio_scatter_arr:
            if not quiet: print('\t\talpha = %.3f; scatter = %.3f' %(spec_index_radio, spec_index_radio_scatter))
            radio_cl_dic[which_dnds][spec_index_radio_scatter]={}
            radio_cl_dic[which_dnds][spec_index_radio_scatter]['cl_dic'] = {}
            for freq1 in freqarr:
                for freq2 in freqarr:
                    cl_analytic = get_poisson_power_with_spectra_index_scatter(freq1, band2 = freq2, min_flux_mJy_band0 = min_flux_mJy, max_flux_mJy_band0 = max_flux_mJy, which_dnds = which_dnds, spec_index_radio = spec_index_radio, spec_index_radio_scatter = spec_index_radio_scatter, dnds_ip_factor = dnds_ip_factor)
                    cl_analytic = np.tile(cl_analytic, len(els) )
                    radio_cl_dic[which_dnds][spec_index_radio_scatter]['cl_dic'][(freq1, freq2)] = cl_analytic
            #print(radio_cl_dic[which_dnds][spec_index_radio_scatter].keys())


            tmp_dic_dic = {}
            tmp_dic_dic['TT'] = radio_cl_dic[which_dnds][spec_index_radio_scatter]['cl_dic']
            res_ilc = ilc.get_ilc_residual_using_weights(tmp_dic_dic, wl1, freqarr, wl2 = wl2, el = els)
            radio_cl_dic[which_dnds][spec_index_radio_scatter]['res_ilc'] = res_ilc

    return radio_cl_dic

##############
def fn_format_axis(ax,fx,fy):
    for label in ax.get_xticklabels(): label.set_fontsize(fx)
    for label in ax.get_yticklabels(): label.set_fontsize(fy)
    if (0):
        import matplotlib.ticker as ticker
        ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
        ax.xaxis.set_minor_locator(ticker.MaNLocator(4))
        ax.yaxis.set_major_locator(ticker.MaNLocator(4))
        ax.yaxis.set_minor_locator(ticker.MaNLocator(4))
    return ax    
