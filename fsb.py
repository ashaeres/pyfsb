# author: Lea Harscouet lea.harscouet@physics.ox.ac.uk
# license: ?

'''Pymaster extension template'''

######################################################


# a function that takes the original map and outputs the filtered squared nmt fields
# a function for the FSB - flag for binned or unbinned (depending on usage: cov or no)
# a function for the gaussian covariance of the FSB
# a function for n222 and n32
    # n32 will be a bit tricky bcs it takes really long 
    # -- C equivalent? maybe for later 
# do we need classes? or just use the Namaster objects


import pymaster as nmt 
import numpy as np
import healpy as hp 


def reduce2(bigm):
    minus1 = np.hstack(tuple(l for l in bigm)); print(minus1.shape)
    minus2 = np.hstack(tuple(l for l in minus1)); print(minus2.shape)
    return minus2

def get_filters(nbands, nside):
    start=2
    nell = (3*nside-1) - start
    dell = nell//nbands
    fls = np.zeros([nbands, 3*nside]) 
    for i in range(nbands):
        fls[i, i*dell+start:(i+1)*dell+start] = 1.0
    return fls

def filtered_sq_fields(map1, filters, mask1, rmask, niter=3, subtractmean=True):
    # maybe change nbands to an array of filters (more flexibility for user)

    npix = len(map1) 
    nside = hp.npix2nside(npix)

    if mask1 is None:
        mask1 = np.ones(npix)
    if rmask is None:
        rmask = mask1 
    
    mask1_bin = mask1>0
    map1 = map1*mask1_bin       

    alm1 = hp.map2alm(map1, lmax=3*nside-1, iter=niter)
    
    mp_filt_sq = np.array([hp.alm2map(hp.almxfl(alm1, fl), nside, lmax=3*nside-1)**2 for fl in filters])   
    if subtractmean == True:
        mp_filt_sq = np.array([(m-np.mean(m[mask1_bin==1]))*mask1_bin for m in mp_filt_sq])   

    f1sq = [nmt.NmtField(rmask, [m], masked_on_input=False, n_iter=niter) for m in mp_filt_sq] 
    
    return np.array(f1sq)


def get_cls_field(field1, mask1=None, field2=None, mask2=None, wksp=None, ells_per_bin=10, cov=False, nside=256):

    if field2 is None:
        field2 = field1
    if mask2 is None:
        mask2 = mask1

    if field1.shape[0]>1: # several fields as input
        print(field1.shape, len(field1.shape))

        if cov is False: # auto power spectra, binned
            if wksp is None:
                wksp = return_wkspace(mask1, mask2, ells_per_bin=ells_per_bin)
            claa = np.array([wksp.decouple_cell(nmt.compute_coupled_cell(fi, fi))[0] for fi in field1])

        else: # cross power spectra, unbinned
            fsky = np.mean(mask1*mask2)
            claa = np.zeros((len(field1), len(field1), 3*nside)) # ?

            if field2 is None:
                for n in range(len(field1)):
                    for m in range(n, len(field2)): # if field2!=field1, should not start from n
                        cross = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
                        claa[n, m] = cross; claa[m, n] = cross
            else:
                for n in range(len(field1)):
                    for m in range(len(field2)): # if field2!=field1, should not start from n
                        claa[n, m] = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
                                
        return claa
    
    else: # one field as input

        if cov is False:
            if wksp is None:
                wksp = return_wkspace(mask1, mask2, ells_per_bin=ells_per_bin)
            clbb = wksp.decouple_cell(nmt.compute_coupled_cell(field1, field2))[0]
        else: 
            fsky = np.mean(mask1*mask2)
            clbb = nmt.compute_coupled_cell(field1, field2)[0] / fsky

        return clbb


def return_wkspace(mask1, mask2=None, ells_per_bin=10): # , mxm=True

    npix = len(mask1); nside = hp.npix2nside(npix)

    if mask2 is None:
        mask2 = mask1 
    
    fmask1 = nmt.NmtField(mask1, None, spin=0)
    fmask2 = nmt.NmtField(mask2, None, spin=0)

    if ells_per_bin==1:
        ells = np.arange(3*nside)
        b = nmt.NmtBin(bpws=ells, ells=ells, weights=np.ones(len(ells)))
    else:
        b = nmt.NmtBin.from_lmax_linear(3*nside-1, ells_per_bin)
    w12 = nmt.NmtWorkspace()
    w12.compute_coupling_matrix(fmask1, fmask2, b)

    return w12 


def get_fsb(map1, filters, mask1, rmask=None, ells_per_bin=10, map2=None, mask2=None, wksp12=None, cov=False, subtractmean=True, niter=3):

    if rmask is None:
        rmask = mask1
    if map2 is None:
        map2 = map1
    if mask2 is None:
        mask2 = mask1
    
    f1s = filtered_sq_fields(map1, filters, mask1, rmask, subtractmean=subtractmean, niter=niter)
    f2 = nmt.NmtField(mask2, [map2], n_iter=niter)

    if cov is False:
        if wksp12 is None:
            wksp12 = return_wkspace(rmask, mask2, ells_per_bin=ells_per_bin)
        fsb = np.array([wksp12.decouple_cell(nmt.compute_coupled_cell(f1, f2))[0] for f1 in f1s])

    else:
        fsky = np.mean(rmask*mask2)
        fsb = np.array([nmt.compute_coupled_cell(f1, f2) for f1 in f1s])[:, 0] / fsky

    return fsb


def get_gauss_cov(fsbs, cls, cls_fsq, mask, rmask, ells_per_bin=10):
    """
    returns a FSB+Cl covariance matrix.
    """

    npix = len(mask); nside = hp.npix2nside(npix)

    bb = nmt.NmtBin.from_lmax_linear(3*nside-1, ells_per_bin)
    b = len(bb.get_effective_ells())

    fmask = nmt.NmtField(mask, None, spin=0)
    fmask_r = nmt.NmtField(rmask, None, spin=0)

    w_fsb = return_wkspace(rmask, mask2=mask, ells_per_bin=ells_per_bin)
    w_cls = return_wkspace(mask, ells_per_bin=ells_per_bin)

    gauss_cov = np.zeros((len(fsbs)+1, len(fsbs)+1, b, b)) 

    fsbs_cls = np.zeros((len(fsbs)+1, fsbs.shape[1]))
    fsbs_cls[:len(fsbs)] = fsbs
    fsbs_cls[-1] = cls

    # FSB-FSB

    for n in range(len(fsbs_cls)):
        for m in range(n, len(fsbs_cls)):

            clad = fsbs_cls[n]
            clbc = fsbs_cls[m]
            clbd = cls
            
            if m==len(fsbs):
                clac = clad 
                clbc = clbd
            else:
                clac = cls_fsq[n, m]

            if n==len(fsbs):
                clac = cls 
                clad = cls
                clbc = cls
                
            cw = nmt.NmtCovarianceWorkspace()
            cw.compute_coupling_coefficients(fmask, fmask_r, fmask, fmask_r)
            covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [clac], [clad], [clbc], [clbd], w_fsb)
            gauss_cov[n, m] = covij_fsb
            gauss_cov[m, n] = covij_fsb
    
    return reduce2(gauss_cov)


def _get_n222_term(cls, filter, mask1, mask2=None):

    if mask2 is None:
        mask2=mask1
    fsky = np.mean(mask1*mask2)

    cls_f = cls*filter

    lmax = len(cls_f)-1
    beam = np.ones_like(cls_f)
    ls = np.arange(lmax+1)
    bin = nmt.NmtBin.from_lmax_linear(lmax, 1)
    w = nmt.nmtlib.comp_coupling_matrix(0, 0, lmax, lmax, 0, 0, 0, 0, beam, beam, cls_f, bin.bin, 0, -1, -1, -1)
    sum_l1 = nmt.nmtlib.get_mcm(w, (lmax+1)**2).reshape([lmax+1, lmax+1])
    sum_l1 /= (2*ls+1)[None, :]
    cov = sum_l1 * 4 * np.outer(cls_f, cls_f)
    nmt.nmtlib.workspace_free(w)

    return cov / fsky

def get_n222_cov(cls, filters, ells_per_bin, mask1, mask2=None):

    bb = nmt.NmtBin.from_lmax_linear(len(cls)-1, ells_per_bin)
    b = len(bb.get_effective_ells())
    n222 = np.zeros((len(filters)+1, len(filters)+1, b, b))

    for nb in range(len(filters)):
        n222ub = _get_n222_term(cls, filters[nb], mask1, mask2=mask2)
        n222_bin1 = np.array([bb.bin_cell(row) for row in n222ub])
        n222_final = np.array([bb.bin_cell(col) for col in n222_bin1.T]).T
        n222[nb, nb] += n222_final 

    return reduce2(n222)


def _get_general_fsb(map1, filters1, filters2, mask1, map2=None, mask2=None):

    npix = len(map1); nside = hp.npix2nside(npix)

    if map2 is None:
        map2=map1
    if mask2 is None:
        mask2=mask1
    fsky = np.mean(mask1*mask2)

    alms = hp.map2alm(map1)
    maps_F = np.array([hp.alm2map(hp.almxfl(alms, fl), nside) for fl in filters1])
    maps_B = np.array([hp.alm2map(hp.almxfl(alms, bl), nside) for bl in filters2])

    cls_FBxOG = np.zeros((len(filters1), len(filters2), len(filters1[0])))

    for f in range(len(filters1)):
        print(f)
        for b in range(len(filters2)):
            print('\t', b)
            map_FB = maps_F[f]*maps_B[b]
            cls_FBxOG[f, b] = hp.anafast(map_FB, map2) / fsky

    return cls_FBxOG


def get_n32_cov(map1, filters1, filters2, ells_per_bin, mask1, map2=None, mask2=None):

    genfsb = _get_general_fsb(map1, filters1, filters2, mask1, map2=map2, mask2=mask2)

    if mask2 is None:
        mask2=mask1
    fsky = np.mean(mask1*mask2)

    cls = hp.anafast(map1) / fsky # not sure if this should be map1 or map2

    bb = nmt.NmtBin.from_lmax_linear(len(cls)-1, ells_per_bin)
    nbins = len(bb.get_effective_ells())
    assert nbins == len(filters2)
    nbands = len(filters1) 

    # bin general fsb
    binned_FSB_Ll1l2 = np.zeros((nbands, nbins, nbins)) # bin the last column of cls_FBxOG_m
    for n in range(nbands):
        for b in range(nbins):
            binned_FSB_Ll1l2[n, b] = bb.bin_cell(genfsb[n, b])

    # bin filtered cls
    cls_f = np.array([cls*fl for fl in filters1])
    un = np.array([bb.bin_cell(e) for e in cls_f]) 
    deux = np.array([np.repeat([e], nbins, axis=0).T for e in un])

    fsb_gen = np.zeros((nbands+1, nbands+1, nbins, nbins))
    cl_filter = np.zeros((nbands+1, nbands+1, nbins, nbins))

    fsb_gen[nbands, :nbands] = binned_FSB_Ll1l2
    cl_filter[nbands, :nbands] = deux

    n32 = cl_filter*fsb_gen / ( np.array([(2*bb.get_effective_ells()+1)]).T * np.pi * fsky)

    return reduce2(n32)+reduce2(n32).T

