# author: Lea Harscouet lea.harscouet@physics.ox.ac.uk

'''Pymaster extension for Filtered Square Bispectra (FSB)

NB: this is the version where 
- we take the bispectrum of fields aab ie the field showing up in the 
power spectrum, is different from the fields that were filtered and squared.
- we also do Not enforce binary masks anymore. nor do we enforce mask intersections.


Functions
---------

_reduce2
get_filters


Classes
---------

FSB
    Methods
    ---------
    __init__
    return_wkspace
    f1s
    fsb_binned
    cls_12_binned
    datavector
    fsb_unbinned
    cls_11_unbinned
    cls_22_unbinned
    cls_12_unbinned 
    cls_1F1Bx2 
    gauss_cov
    cls_1sq1sq_unbinned 
    fsb_unbinned_pure 
    fsb_binned_pure 
    cls_11_binned
    master_datavector
    cls_datavector
    return_wkspace
    filtered_sq_fields
    get_cls_field
    get_fsb
    get_gauss_cov
    get_gauss_cov_cls
    get_gauss_cov_pure
    get_master_gauss_cov
    _get_n222_term
    ined earlier
    get_n222_cov
    _get_general_fsb
    twonickels
    get_n32_cov
    get_full_cov
'''

######################################################

from functools import cached_property
import pymaster as nmt 
import numpy as np
import healpy as hp 


def _reduce2(bigm):
    """
    Reduces the dimensionality of an array by 2 along its first two axes.

    Arguments
    ----------
    bigm : array with shape (i, j, k, ...)

    Returns
    ----------
    An array of shape (k, ...)

    """
    minus1 = np.hstack(tuple(l for l in bigm)); # print(minus1.shape)
    minus2 = np.hstack(tuple(l for l in minus1)); # print(minus2.shape)
    return minus2


def get_filters(nbands, nside):
    """
    Linearly divides the ell range (2; 3*nside-1)
    into nbands filters of equal size.

    Arguments
    ----------
    nbands : int
        number of filters
    nside : int
        resolution of the healpy map used

    Returns
    ----------
    A binary array of shape (nbands, 3*nside)
    
    """
    start=2
    nell = (3*nside-1) - start
    dell = nell//nbands
    fls = np.zeros([nbands, 3*nside]) 
    for i in range(nbands):
        fls[i, i*dell+start:(i+1)*dell+start] = 1.0
    return fls




class FSB():

    def __init__(self, map1, mask1, filters, map2=None, mask2=None, ells_per_bin=10, niter=3): # , rmask=None

        self.niter = niter

        # maps and masks
        self.map1 = map1

        self.npix = len(self.map1)
        self.nside = hp.npix2nside(self.npix)

        if mask1 is None:
            self.mask1 = np.ones(self.npix)
        else:
            self.mask1 = mask1

        if map2 is None:
            self.map2 = map1
            self.twofields = False
        else:
            self.map2 = map2
            self.twofields = True
        if mask2 is None:
            self.mask2 = mask1
        else:
            self.mask2 = mask2
        
        self.rmask = self.mask1 > 0 # binary version of mask1
        self.rmask2 = self.mask2 > 0 # binary version of mask2 (used in generalised FSB)
        # remasking fields appropriately + need to make sure fields within new mask is 0
        self.map1 = (self.map1-np.mean(self.map1[self.mask1==1]))
        self.map2 = (self.map2-np.mean(self.map2[self.mask2==1]))

        # filters
        self.filters = filters
        self.nbands = len(self.filters)

        # binning # TODO: relabel those to make them more explicit
        self.ells_per_bin = ells_per_bin
        self.lmax = 3*self.nside-1
        self.bb = nmt.NmtBin.from_lmax_linear(self.lmax, self.ells_per_bin)
        self.b = len(self.bb.get_effective_ells())
        self.bins = get_filters(self.b, self.nside) # filters corresponding to bins (for generalized fsb)
        # this is specifically for n222 calculation
        self._ls = np.arange(self.lmax+1)
        self._bin_n222 = nmt.NmtBin.from_lmax_linear(self.lmax, 1)
        self._beam = np.ones(self.lmax+1)


        # workspaces and fsky
        self.w_fsb = self.return_wkspace(self.rmask, self.mask2, self.ells_per_bin)
        self.fsky_fsb = np.mean(self.rmask*self.mask2)

        self.w_fsb_pure = self.return_wkspace(self.rmask, self.mask1, self.ells_per_bin)
        self.fsky_fsb_pure = np.mean(self.rmask*self.mask1)
        
        self.w_cls_11 = self.return_wkspace(self.mask1, self.mask1, self.ells_per_bin)
        self.fsky_cls_11 = np.mean(self.mask1*self.mask1)

        self.w_cls_12 = self.return_wkspace(self.mask1, self.mask2, self.ells_per_bin)
        self.fsky_cls_12 = np.mean(self.mask1*self.mask2)

        self.w_cls_22 = self.return_wkspace(self.mask2, self.mask2, self.ells_per_bin)
        self.fsky_cls_22 = np.mean(self.mask2*self.mask2)

        # self.w_cls_rr = self.return_wkspace(self.rmask, self.rmask, self.ells_per_bin)
        self.fsky_cls_rr = np.mean(self.rmask*self.rmask)

        # fields
        self.field1 = nmt.NmtField(self.mask1, [self.map1], masked_on_input=False, n_iter=self.niter)
        if map2 is None:
            self.field2 = self.field1
        else:
            self.field2 = nmt.NmtField(self.mask2, [self.map2], masked_on_input=False, n_iter=self.niter)

        # self.cls_1F1Bx2 = None # TODO: make default usage with self.filters and self.binfilters
        self._genfsbs = {}

    
    @cached_property
    def f1s(self):
        return self.filtered_sq_fields()

    @cached_property
    def cls_1sq1sq_unbinned(self): 
        return self.get_cls_field(self.f1s) / self.fsky_cls_rr
    
    @cached_property
    def fsb_unbinned(self):
        return self.get_fsb() / self.fsky_fsb
    
    @cached_property
    def fsb_binned(self):
        return self.get_fsb(wksp=self.w_fsb)

    @cached_property
    def fsb_unbinned_pure(self): 
        return self.get_cls_field(self.f1s, field2=np.array([self.field1])) / self.fsky_fsb_pure
    
    @cached_property
    def fsb_binned_pure(self): 
        return self.get_cls_field(self.f1s, field2=np.array([self.field1]), wksp=self.w_fsb_pure) 
    
    @cached_property
    def cls_11_unbinned(self):
        return self.get_cls_field(np.array([self.field1])) / self.fsky_cls_11
    
    @cached_property
    def cls_11_binned(self):
        return self.get_cls_field(np.array([self.field1]), wksp=self.w_cls_11)
    
    @cached_property
    def cls_12_unbinned(self): 
        return self.get_cls_field(np.array([self.field1]), field2=np.array([self.field2])) / self.fsky_cls_12
    
    @cached_property
    def cls_12_binned(self):
        return self.get_cls_field(np.array([self.field1]), field2=np.array([self.field2]), wksp=self.w_cls_12)

    @cached_property
    def cls_22_unbinned(self):
        return self.get_cls_field(np.array([self.field2])) / self.fsky_cls_22
    
    @cached_property
    def datavector(self):
        return np.concatenate((self.fsb_binned.flatten(), self.cls_12_binned))

    @cached_property
    def cls_datavector(self):
        return np.concatenate((self.cls_11_binned, self.cls_12_binned))

    @cached_property
    def master_datavector(self):
        return np.concatenate((self.fsb_binned_pure.flatten(), self.cls_11_binned, self.fsb_binned.flatten(), self.cls_12_binned))
    
    @cached_property
    def cov_cls(self):
        return self.get_cov_cls()
    
    @cached_property
    def cov_auto_gauss(self):
        return self.get_cov_auto(n222=False, n32=False)
    
    @cached_property
    def cov_cross_gauss(self):
        return self.get_cov_cross(n222=False, n32=False)
    

    def return_wkspace(self, mask1, mask2, lpb): 

        """
        Creates an `NmtWorkspace` object from mask(s) and
        a binning scheme.

        Arguments
        ----------
        mask1 : array
            healpy map of size (12*nside**2)
        mask2 : array
            a healpy map of the same size as mask1
        ells_per_bin : int (default = 10)
            the number of ells in each bin in the binning scheme.

        Returns
        ----------
        A `NmtWorkspace` object corresponding to the
        given mask(s) and binning scheme.
        
        """
            
        fmask1 = nmt.NmtField(mask1, None, spin=0)
        fmask2 = nmt.NmtField(mask2, None, spin=0)

        b = nmt.NmtBin.from_lmax_linear(self.lmax, lpb)

        w12 = nmt.NmtWorkspace()
        w12.compute_coupling_matrix(fmask1, fmask2, b)

        return w12 



    def filtered_sq_fields(self):

        """
        Creates NmtField objects for each filtered-squared map.

        Returns
        ----------
        An array containing as many NmtField objects as 
        there are filters.
        
        """ 
        # # already done up there i believe
        # mask1_bin = self.effmask>0   
        # map1 = self.map1*mask1_bin  

        alm1 = hp.map2alm(self.map1, iter=self.niter)
        
        mp_filt_sq = np.array([hp.alm2map(hp.almxfl(alm1, fl), self.nside, lmax=self.lmax)**2 for fl in self.filters])  

        # # FIXME: sara's fix (potentially affects covariance)
        # # print(np.mean(mp_filt_sq, axis=1))
        # mp_filt_sq = np.array([(mp_filt_sq[i]-np.average(mp_filt_sq[i], weights=self.rmask))*self.rmask for i in range(len(mp_filt_sq))]) 
        
        f1sq = [nmt.NmtField(self.rmask, [m], masked_on_input=False, n_iter=self.niter) for m in mp_filt_sq] 
        
        return np.array(f1sq)
    
    

    def get_cls_field(self, field1, field2=None, wksp=None):

        """
        Computes the power spectra of the given field(s).
        Handles the following cases:
            - if field1 is a single field, will return 
            its auto-power spectrum
            - if field1 is several fields, will return 
            the corresponding cross-power spectra
            - if field1 and field2 are both single fields,
            will return their cross-power spectrum.
            - if both field1 and field2 are several fields,
            will return their cross-power spectra. (NOT SURE THIS WORKS)

        Arguments
        ----------
        field1 : `NmtField` object or array of `NmtField` objects
        field2 : `NmtField` object or array of `NmtField` objects,
        optional
        wksp : `NmtWorkspace` object, optional
            a workspace to compute the mode coupling and binning of
            the power spectra. If not given, it will be assumed that
            no binning is required, and the fsky correction is applied
            to the unbinned cls.

        Returns
        ----------
        An array of the cls.
            
        """

        if field2 is None:
            field2 = field1
            same = True
        else:
            same = False

        if field1.shape[0]>1: # several fields as input

            if wksp is None: # cross power spectra, unbinned

                claa = np.zeros((len(field1), len(field2), 3*self.nside)) 

                if same:
                    for n in range(len(field1)):
                        for m in range(n, len(field2)):
                            cross = nmt.compute_coupled_cell(field1[n], field2[m])[0] 
                            claa[n, m] = cross; claa[m, n] = cross
                else:
                    for n in range(len(field1)):
                        for m in range(len(field2)): # if field2!=field1, should not start from n?
                            claa[n, m] = nmt.compute_coupled_cell(field1[n], field2[m])[0] 
                            # TODO: are we using this bit in the covariance? why /fsky?
                                      
            else: 

                if same: # auto power spectra, binned
                    claa = np.array([wksp.decouple_cell(nmt.compute_coupled_cell(fi, fi))[0] for fi in field1])

                else: # cross power spectra, binned
                    claa = np.zeros((len(field1), len(field2), self.b)) 
                    for n in range(len(field1)):
                        for m in range(len(field2)): # if field2!=field1, should not start from n?
                            claa[n, m] = wksp.decouple_cell(nmt.compute_coupled_cell(field1[n], field2[m]))[0]
                    
            return claa.squeeze() 
        
        else: # one field as input, inside a np.array

            if wksp is None: # auto power spectra, unbinned
                clbb = nmt.compute_coupled_cell(field1[0], field2[0])[0]

            else: # auto power spectra, binned
                clbb = wksp.decouple_cell(nmt.compute_coupled_cell(field1[0], field2[0]))[0]

            return clbb



    def get_fsb(self, wksp=None):

        """
        Computes the FSB of the class field(s)
        for a set of filters. if map1 and map2 are different, 
        only map1 will be filtered.
        
        Arguments
        ----------
        
        wksp12 : `NmtWorkspace` object, optional
            a workspace to compute the mode coupling and binning of
            the power spectra. if not given, it will be assumed that
            no binning is required, and the fsky correction is applied
            to the unbinned FSB.
        
        Returns
        ----------
        An array of the FSBs.
            
        """

        return self.get_cls_field(self.f1s, field2=np.array([self.field2]), wksp=wksp) 
        

    def get_cov_cross(self, n222=True, n32=False, insquares=True):

        # if hasattr(self, 'cov_cross_gauss') is False:
        
        fmask_r = nmt.NmtField(self.rmask, None, spin=0)
        fmask_1 = nmt.NmtField(self.mask1, None, spin=0)
        fmask_2 = nmt.NmtField(self.mask2, None, spin=0)
        # TODO: maybe make them attributes?
        cw_fsbfsb = nmt.NmtCovarianceWorkspace(fmask_r, fmask_2, fmask_r, fmask_2)
        cw_fsbcls = nmt.NmtCovarianceWorkspace(fmask_r, fmask_2, fmask_1, fmask_2)
        cw_clscls = nmt.NmtCovarianceWorkspace(fmask_1, fmask_2, fmask_1, fmask_2)

        gauss_cov = np.zeros((self.nbands+1, self.nbands+1, self.b, self.b)) 

        for n in range(self.nbands):

            # this is for the fsb-cls cross covariance
            cla2b1 = self.cls_12_unbinned #/ self.fsky_cls_12
            cla2b2 = self.cls_22_unbinned #/ self.fsky_cls_22

            if self.nbands==1:
                cla1b1 = self.fsb_unbinned_pure #/ self.fsky_fsb_pure
                cla1b2 = self.fsb_unbinned #/ self.fsky_fsb 
            else:
                cla1b1 = self.fsb_unbinned_pure[n] #/ self.fsky_fsb_pure
                cla1b2 = self.fsb_unbinned[n] #/ self.fsky_fsb 
            
            covij_fsb = cw_fsbcls.gaussian_covariance(
                # cw_fsbcls, 0, 0, 0, 0, 
                [cla1b1], [cla1b2], [cla2b1], [cla2b2], 
                self.w_fsb, self.w_cls_12)
            gauss_cov[n, -1] = covij_fsb
            gauss_cov[-1, n] = covij_fsb.T # TODO: changed

            for m in range(n, self.nbands):

                # this is for the fsb-fsb covariance

                if self.nbands==1:
                    cla1b1 = self.cls_1sq1sq_unbinned #/ self.fsky_cls_rr
                    cla1b2 = self.fsb_unbinned #/ self.fsky_fsb 
                else:
                    cla1b1 = self.cls_1sq1sq_unbinned[n,m] #/ self.fsky_cls_rr
                    cla2b1 = self.fsb_unbinned[m] #/ self.fsky_fsb 

                covij_fsb = cw_fsbfsb.gaussian_covariance(
                    # cw_fsbfsb, 0, 0, 0, 0, 
                    [cla1b1], [cla1b2], [cla2b1], [cla2b2], 
                    self.w_fsb)
                gauss_cov[n, m] = covij_fsb
                gauss_cov[m, n] = covij_fsb.T # TODO: new transpose

        gauss_cov[-1, -1] = cw_clscls.gaussian_covariance(
                                                    # cw_clscls, 0, 0, 0, 0, 
                                                    [self.cls_11_unbinned], # /self.fsky_cls_11
                                                    [self.cls_12_unbinned], [self.cls_12_unbinned], # /self.fsky_cls_12
                                                    [self.cls_22_unbinned], 
                                                    self.w_cls_12) # /self.fsky_cls_22
        self.cov_cross = gauss_cov
        
        # else:
        #     self.cov_cross = self.cov_cross_gauss

        if n222:
            self.cov_cross += self.get_n222_cov(self.cls_12_unbinned, self.cls_12_unbinned, self.fsky_fsb)
        if n32: 
            temp = self.get_n32_cov(self.cls_11_unbinned, '122', self.filters, self.bins, self.fsky_fsb, self.cls_12_unbinned, '112')
            # print(temp.shape, self.cov_cross.shape)
            self.cov_cross += temp
        if insquares==False:
            return _reduce2(self.cov_cross)
        else:
            return self.cov_cross
        
    
    def get_cov_cls(self, insquares=True):

        fmask_1 = nmt.NmtField(self.mask1, None, spin=0)
        fmask_2 = nmt.NmtField(self.mask2, None, spin=0)

        cw_1212 = nmt.NmtCovarianceWorkspace(fmask_1, fmask_2, fmask_1, fmask_2) # gk, gk
        cw_1211 = nmt.NmtCovarianceWorkspace(fmask_1, fmask_1, fmask_1, fmask_2) # gg, gk
        cw_1111 = nmt.NmtCovarianceWorkspace(fmask_1, fmask_1, fmask_1, fmask_1) # gg, gg

        gauss_cov_cls = np.zeros((2, 2, self.b, self.b)) 

        cls = np.zeros((3, self.cls_11_unbinned.shape[0]))
        cls[0] = self.cls_11_unbinned #/ self.fsky_cls_11 # gg
        cls[1] = self.cls_12_unbinned #/ self.fsky_cls_12 # gk
        cls[2] = self.cls_22_unbinned #/ self.fsky_cls_22 # kk

        # TODO: careful: below is not the same ordering as in the master cov :/
        gauss_cov_cls[0,0] = cw_1111.gaussian_covariance([cls[0]], [cls[0]], [cls[0]], [cls[0]], self.w_cls_11) # gg, gg # cw_1111, 0, 0, 0, 0, 
        gauss_cov_cls[1,1] = cw_1212.gaussian_covariance([cls[0]], [cls[1]], [cls[1]], [cls[2]], self.w_cls_12) # gk, gk # cw_1212, 0, 0, 0, 0, 
        gauss_cov_cls[0,1] = cw_1211.gaussian_covariance([cls[0]], [cls[1]], [cls[0]], [cls[1]], self.w_cls_11, self.w_cls_12) # gg, gk # cw_1211, 0, 0, 0, 0, 
        gauss_cov_cls[1,0] = gauss_cov_cls[0,1].T # symmetric
        
        self.cov_cls = gauss_cov_cls

        if insquares==False:
            return _reduce2(self.cov_cls)
        else:
            return self.cov_cls
        

    def get_cov_auto(self, n222=True, n32=False, insquares=True):

        # if hasattr(self, 'cov_auto_gauss') is False:
    
        fmask_r = nmt.NmtField(self.rmask, None, spin=0)
        fmask_1 = nmt.NmtField(self.mask1, None, spin=0)

        cw_fsbfsb = nmt.NmtCovarianceWorkspace(fmask_r, fmask_1, fmask_r, fmask_1)
        cw_fsbcls = nmt.NmtCovarianceWorkspace(fmask_1, fmask_r, fmask_1, fmask_1)
        cw_clscls = nmt.NmtCovarianceWorkspace(fmask_1, fmask_1, fmask_1, fmask_1)

        gauss_cov = np.zeros((self.nbands+1, self.nbands+1, self.b, self.b)) 

        for n in range(self.nbands):

            cla2b1 = self.cls_11_unbinned #/ self.fsky_cls_11
            cla2b2 = self.cls_11_unbinned #/ self.fsky_cls_11

            if self.nbands==1:
                cla1b1 = self.fsb_unbinned_pure #/ self.fsky_fsb_pure
                cla1b2 = self.fsb_unbinned_pure #/ self.fsky_fsb_pure
            else:
                cla1b1 = self.fsb_unbinned_pure[n] #/ self.fsky_fsb_pure
                cla1b2 = self.fsb_unbinned_pure[n] #/ self.fsky_fsb_pure
            
            covij_fsb = cw_fsbcls.gaussian_covariance([cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb_pure, self.w_cls_11) # cw_fsbcls, 0, 0, 0, 0, 
            gauss_cov[n, -1] = covij_fsb
            gauss_cov[-1, n] = covij_fsb.T # TODO: changed

            for m in range(n, self.nbands):

                if self.nbands==1:
                    cla1b1 = self.cls_1sq1sq_unbinned #/ self.fsky_cls_rr 
                    cla2b1 = self.fsb_unbinned_pure #/ self.fsky_fsb_pure
                else:
                    cla1b1 = self.cls_1sq1sq_unbinned[n,m] #/ self.fsky_cls_rr 
                    cla2b1 = self.fsb_unbinned_pure[m] #/ self.fsky_fsb_pure

                covij_fsb = cw_fsbfsb.gaussian_covariance([cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb_pure) # cw_fsbfsb, 0, 0, 0, 0, 
                gauss_cov[n, m] = covij_fsb
                gauss_cov[m, n] = covij_fsb.T # TODO: new .T, let's see if it makes things better

        gauss_cov[-1, -1] = cw_clscls.gaussian_covariance([self.cls_11_unbinned], # /self.fsky_cls_11
                                                    [self.cls_11_unbinned], [self.cls_11_unbinned], 
                                                    [self.cls_11_unbinned], self.w_cls_11) # cw_clscls, 0, 0, 0, 0, 
        
        self.cov_auto = gauss_cov

        # else:
        #     self.cov_auto = self.cov_auto_gauss

        if n222:
            self.cov_auto += self.get_n222_cov(self.cls_11_unbinned, self.cls_11_unbinned, self.fsky_fsb_pure)
        if n32:
            self.cov_auto += self.get_n32_cov(self.cls_11_unbinned, '111', self.filters, self.bins, self.fsky_fsb_pure) 
        
        if insquares==False:
            return _reduce2(self.cov_auto)
        else:
            return self.cov_auto
        


    def get_cov_all(self, n222=True, n32=False, insquares=True):

        """TODO: add description of shape of matrix"""

        fmask_r = nmt.NmtField(self.rmask, None, spin=0)
        fmask_1 = nmt.NmtField(self.mask1, None, spin=0)
        fmask_2 = nmt.NmtField(self.mask2, None, spin=0)

        cw_ggkggg = nmt.NmtCovarianceWorkspace(fmask_r, fmask_2, fmask_r, fmask_1)
        cw_ggkgg = nmt.NmtCovarianceWorkspace(fmask_r, fmask_2, fmask_1, fmask_1)
        cw_ggggk = nmt.NmtCovarianceWorkspace(fmask_r, fmask_1, fmask_1, fmask_2)

        # structure of data vector if fsb_ggk + cl_gk + fsb_ggg + cl_gg
        master_cov = np.zeros((2*(self.nbands+1), 2*(self.nbands+1), self.b, self.b)) 
        
        oldggg = self.get_cov_auto(n222=n222, n32=n32) 
        oldggk = self.get_cov_cross(n222=n222, n32=n32) 
        oldcls = self.cov_cls
        master_cov[self.nbands+1:, self.nbands+1:] = oldggk
        master_cov[:self.nbands+1, :self.nbands+1] = oldggg
        master_cov[self.nbands, -1] = oldcls[0,1] # gg, gk
        master_cov[-1, self.nbands] = oldcls[1,0] # gk, gg # one of them will be overwritten by transpose bit, but not sure which one

        # now we need master_cov[self.nbands+1:, :self.nbands+1]
        fsbs_cls_mixed = np.zeros((self.nbands+1, len(self.cls_11_unbinned)))
        fsbs_cls_mixed[:self.nbands] = self.fsb_unbinned #/ self.fsky_fsb
        fsbs_cls_mixed[-1] = self.cls_12_unbinned #/ self.fsky_cls_12
        fsbs_cls_pure = np.zeros_like(fsbs_cls_mixed)
        fsbs_cls_pure[:self.nbands] = self.fsb_unbinned_pure #/ self.fsky_fsb_pure
        fsbs_cls_pure[-1] = self.cls_11_unbinned #/ self.fsky_cls_11

        for n in range(self.nbands):

            if self.nbands==1:
                cla1b1 = fsbs_cls_pure[0]
                cla1b2 = fsbs_cls_pure[-1]
                cla2b1 = fsbs_cls_mixed[0]
                cla2b2 = fsbs_cls_mixed[-1]
            else:
                cla1b1 = fsbs_cls_pure[n] # for both cases
                cla1b2 = fsbs_cls_pure[-1]
                cla2b1 = fsbs_cls_mixed[n]
                cla2b2 = fsbs_cls_mixed[-1]

            # ggg, gk (along axis 1, horizontal)
            covij_fsb = cw_ggggk.gaussian_covariance([cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb_pure, self.w_cls_12) # cw_ggggk, 0, 0, 0, 0, 
            master_cov[2*(self.nbands+1)-1, n] = covij_fsb
            # ggk, gg (along axis 0, vertical)
            covij_fsb = cw_ggkgg.gaussian_covariance([cla1b1], [cla1b1], [cla2b2], [cla2b2], self.w_fsb, self.w_cls_11) # cw_ggkgg, 0, 0, 0, 0, 
            master_cov[(self.nbands+1) + n, (self.nbands)] = covij_fsb

            for m in range(self.nbands):

                if self.nbands==1:
                    cla1b1 = self.cls_1sq1sq_unbinned #/ self.fsky_cls_rr
                    cla1b2 = fsbs_cls_pure[0]
                    cla2b1 = fsbs_cls_mixed[0]
                else:
                    cla1b1 = self.cls_1sq1sq_unbinned[n,m] #/ self.fsky_cls_rr
                    cla1b2 = fsbs_cls_pure[n]
                    cla2b1 = fsbs_cls_mixed[m]
                    # cla2b2 = fsbs_cls_mixed[-1] # same as above
    
                covij_fsb = cw_ggkggg.gaussian_covariance([cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb) # cw_ggkggg, 0, 0, 0, 0, 
                master_cov[(self.nbands+1) + n, m] = covij_fsb
                # master_cov[m, (self.nbands+1) + n] = covij_fsb # NO! not symmetric bcs Phi_ggg and Phi_ggk

        if n222:
            master_cov[self.nbands+1:, :self.nbands+1] += self.get_n222_cov(self.cls_11_unbinned, self.cls_12_unbinned, self.fsky_cls_11)
        if n32:
            # compute ggg x gk first
            master_cov[self.nbands+1:, :self.nbands+1] += self.get_n32_cov(self.cls_11_unbinned, '112', self.filters, self.bins, self.fsky_fsb_pure, self.cls_12_unbinned, '111', symmetric=False)
            # then ggk x gg
            temp = self.get_n32_cov(self.cls_11_unbinned, '112', self.filters, self.bins, self.fsky_fsb_pure, symmetric=False) # TODO: can make this faster by saving cl-genfsb combinations?
            # and transpose that one
            temp_t = np.transpose(temp, (1, 0, 3, 2))
            master_cov[self.nbands+1:, :self.nbands+1] += temp_t
                    
        transposed = np.transpose(master_cov[self.nbands+1:, :self.nbands+1], (1, 0, 3, 2)) # used to be (1, 0, 3, 2)
        master_cov[:self.nbands+1, self.nbands+1:] = transposed 

        self.master_cov = master_cov

        if insquares==False:
            return _reduce2(self.master_cov)
        else:
            return self.master_cov



    def _get_n222_term(self, filter, cl_out_one, cl_out_two):

        """
        Computes the N222 term for each FSB.

        Arguments
        ----------
        cls : array 
            an array containing the power spectrum of the
            original map (before filtering/squaring)
        filter : array
            a single binary array of shape (3*nside, )

        Returns
        ----------
        The not yet mask-corrected N222 term in a
        block matrix of size (len(cls), len(cls)).
        """

        # the cls within the sum are always the same for ggg x ggg, ggg x ggk, ggk x ggk
        # the cls outside varies for each of the three cases listed above
        cls_in_f = self.cls_11_unbinned*filter #/ self.fsky_cls_11
        cls_out1_f = cl_out_one*filter
        cls_out2_f = cl_out_two*filter

        w = nmt.nmtlib.comp_coupling_matrix(0, 0, self.lmax, self.lmax, 0, 0, 0, 0, 0,
                                            0.0, self._beam, self._beam, cls_in_f,
                                            self._bin_n222.bin, 0, -1, -1, -1)
        sum_l1 = nmt.nmtlib.get_mcm(w, (self.lmax+1)**2).reshape([self.lmax+1, self.lmax+1])
        sum_l1 /= (2*self._ls+1)[None, :]
        cov = sum_l1 * 4 * np.outer(cls_out1_f, cls_out2_f)
        nmt.nmtlib.workspace_free(w)

        return cov

    def get_n222_cov(self, cl_out_one, cl_out_two, fskycorrection):

        """
        Computes the mask-corrected, binned N222
        terms for all FSBs.

        Returns
        ----------
        The binned mask-corrected N222 term in an
        array of dimensions (nbands+1, nbands+1, nbins).
        """

        n222 = np.zeros((len(self.filters)+1, len(self.filters)+1, self.b, self.b))

        for nb in range(self.nbands):
            n222ub = self._get_n222_term(self.filters[nb], cl_out_one, cl_out_two)
            n222_bin1 = np.array([self.bb.bin_cell(row) for row in n222ub])
            n222_final = np.array([self.bb.bin_cell(col) for col in n222_bin1.T]).T
            n222[nb, nb] += n222_final 

        return n222 / fskycorrection


    def _get_general_fsb(self, id, m1, m2, m3, filters1, filters2): 

        """
        Computes the generalized FSB for 
        two sets of filters.

        Arguments
        ----------
        id : string
            the map combination of the general fsb: 
            '111', '112', '122' is all one should need for 
            covariance calculations
        m1, m2, m3: NmtField objects
            field corresponding to the three fields
            indicated in the id argument: if id is 
            '112', pass m1, m2, m3 = (self.field1, 
            self.field1, self.field2)
            # TODO: could just get m1, m2, m3 from id
        filters1 : array 
            a binary array of shape (nbands1, 3*nside)
        filters2 : array 
            a binary array of shape (nbands2, 3*nside)
        
        Returns
        ----------
        The binned generalized FSBs in an array
        of shape (nbands1, nbands2, 3*nside).
        """

        # yk what let's pass it the fields, rather than the maps only (at least we got alms and masks for free)

        print(f'computing {id} hopefully for the first and only time')

        localfsky = np.mean(m1.mask * m2.mask * m3.mask) # TODO: hopefully correct?

        maps_F = np.array([hp.alm2map(hp.almxfl(m1.alm[0], fl), self.nside)*m1.mask for fl in filters1])
        maps_B = np.array([hp.alm2map(hp.almxfl(m2.alm[0], bl), self.nside)*m2.mask for bl in filters2])

        genfsb = np.zeros((len(filters1), len(filters2), len(filters1[0])))

        for f in range(len(filters1)):
            print('\t', f)
            for b in range(len(filters2)):
                print('\t\t', b)
                map_FB = maps_F[f]*maps_B[b]
                genfsb[f, b] = hp.anafast(map_FB, m3.maps[0]) / (localfsky * self.ells_per_bin)

        self._genfsbs[id] = genfsb



    def twonickels(self, cls, genfsb, filters1, filters2):
        """binning and multiplying the arguments
        if i had a nickel everytime this function was used, 
        i'd have two nickels - which is not a lot but 
        it's weird it's happened twice"""

        # bin general fsb
        binned_genFSB = np.zeros((len(filters1), len(filters2), self.b)) 
        for n in range(len(filters1)):
            for b in range(len(filters2)):
                binned_genFSB[n, b] = self.bb.bin_cell(genfsb[n, b])
        # bin filtered cls
        cls_f = np.array([cls*fl for fl in filters1])
        un = np.array([self.bb.bin_cell(e) for e in cls_f]) 
        deux = np.array([np.repeat([e], len(filters2), axis=0).T for e in un])

        fsb_gen = np.zeros((len(filters1)+1, len(filters1)+1, len(filters2), self.b))
        cl_filter = np.zeros((len(filters1)+1, len(filters1)+1, len(filters2), self.b))

        fsb_gen[len(filters1), :len(filters1)] = binned_genFSB
        cl_filter[len(filters1), :len(filters1)] = deux

        return cl_filter*fsb_gen



    def get_n32_cov(self, cls1, id1, filters1, filters2, fskycorrection, cls2=None, id2=None, symmetric=True):

        """
        Computes the mask-corrected, binned N32
        terms for all FSBs x Cls.

        Arguments
        ----------
        cls1: array
        id1: string
        cls2: array
        id2: string
        filters1 : array
            a binary array of shape (nbands1, 3*nside)
        filters2 : array 
            a binary array of shape (nbands2, 3*nside)
        fskycorrection: 


        Returns
        ----------
        The binned mask-corrected N32 term in an
        array of dimensions (nbands+1, nbands+1, nbins).
        (if we assume filters2 to be the equivalent of 
        the binning scheme.)
        """
        
        if self._genfsbs.__contains__(id1) is False:
            ms = []
            for t in id1:
                if t=='1':
                    ms.append(self.field1)
                else:
                    ms.append(self.field2)
            self._get_general_fsb(id1, *ms, filters1, filters2)
            
            
        n32_term = self.twonickels(cls1, self._genfsbs[id1], filters1, filters2)

        if id2 is not None:
            if self._genfsbs.__contains__(id2) is False:
                ms = []
                for t in id2:
                    if t=='1':
                        ms.append(self.field1)
                    else:
                        ms.append(self.field2)
                self._get_general_fsb(id2, *ms, filters1, filters2)
            n32_term += self.twonickels(cls2, self._genfsbs[id2], filters1, filters2)
        else:
            n32_term += n32_term # so we can just do 2* (instead of either 2* or 4*)
            
        n32 = 2* n32_term / ( np.array([(2*self.bb.get_effective_ells()+1)]).T )

        if symmetric:
            for i in range(self.nbands): # make it symmetric
                n32[i, self.nbands] = n32[self.nbands, i].T

        return n32 / fskycorrection
    
    

    # def get_full_cov(self, insquares=False, n32=True):

    #     """
    #     Adds the mask-corrected, binned gaussian-
    #     limit covariance and its main mask-corrected, 
    #     binned non-gaussian contributions. 

    #     Returns
    #     ----------
    #     An array of dimensions ((nbands+1)*nbins, (nbands+1)*nbins),
    #     or if insquares set to True, an array of dimensions 
    #     ((nbands+1), (nbands+1), nbins, nbins).
    #     """
    #     # if self.gauss_cov is None:
    #     #     self.gauss_cov = self.get_gauss_cov()

    #     cl_out = self.cls_12_unbinned #/self.fsky_cls_12

    #     if n32:
    #         self.full_cov_large = self.gauss_cov + self.get_n222_cov(cl_out, cl_out, self.fsky_cls_12) + self.get_n32_cov(self.filters, self.bins)
    #     else: # bypass n32 altogether if not needed
    #         self.full_cov_large = self.gauss_cov + self.get_n222_cov(cl_out, cl_out, self.fsky_cls_12)
        
    #     if insquares:
    #         return self.full_cov_large
    #     else:
    #         return _reduce2(self.full_cov_large)





