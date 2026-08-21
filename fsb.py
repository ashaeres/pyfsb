# author: Lea Harscouet lea.harscouet@physics.ox.ac.uk

'''Pymaster extension for Filtered Square Bispectra (FSB)

version in which we accomodate for an array input for map1

we also need a version in which map2 can be the array:(
(thinking of gaussian planck realisations to compute ggk cov)

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
    filtered_sq_fields
    get_cls_field
    get_fsb
    get_gauss_cov
    _get_n222_term
    get_n222_cov
    _get_general_fsb
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








# class FSBmap1():

#     def __init__(self, map1, mask1, filters, map2=None, mask2=None, ells_per_bin=10, niter=3): # , rmask=None

#         self.map1 = map1

#         self.npix = len(self.map1[0])
#         self.nside = hp.npix2nside(self.npix)

#         if mask1 is None:
#             self.mask1 = np.ones(self.npix)
#         else:
#             self.mask1 = mask1

#         if map2 is None:
#             self.map2 = map1
#             self.twofields = False
#         else:
#             self.map2 = map2
#             self.twofields = True
#         if mask2 is None:
#             self.mask2 = mask1
#         else:
#             self.mask2 = mask2
        

#         # set the effective mask to the intersection of the 2 masks
#         self.effmask = (self.mask1*self.mask2) > 0 # make binary
#         self.rmask = self.effmask
#         # remasking fields appropriately + need to make sure fields within new mask is 0
#         self.map1 = np.array([m-np.mean(m[self.effmask==1]) for m in self.map1])*self.effmask # FIXME:
#         self.map2 = np.array([m-np.mean(m[self.effmask==1]) for m in self.map2])*self.effmask # FIXME:

#         self.filters = filters
#         self.ells_per_bin = ells_per_bin
#         self.niter = niter

#         # binning
#         self.bb = nmt.NmtBin.from_lmax_linear(3*self.nside-1, self.ells_per_bin)
#         self.b = len(self.bb.get_effective_ells())
#         self.bins = get_filters(self.b, self.nside) # filters corresponding to bins (for generalized fsb)
        
        
#         self.w_fsb = self.return_wkspace(self.rmask, self.effmask, self.ells_per_bin)
#         self.w_cls = self.return_wkspace(self.effmask, self.effmask, self.ells_per_bin)

#         self.fsky_fsb = np.mean(self.rmask*self.effmask)
#         self.fsky_cls = np.mean(self.effmask*self.effmask)

#         self.nbands = len(self.filters)

#         self.field1 = np.array([nmt.NmtField(self.effmask, [m], masked_on_input=False, n_iter=self.niter) for m in self.map1]) # FIXME:
#         if map2 is None:
#             self.field2 = self.field1
#         else:
#             self.field2 = np.array([nmt.NmtField(self.effmask, [m], masked_on_input=False, n_iter=self.niter) for m in self.map2]) # FIXME:

#         self.cls_1F1Bx2 = None # TODO: make default usage with self.filters and self.binfilters

    
    
#     @cached_property
#     def f1s(self):
#         print('computed f1s for the 1st (and hopefully only) time')
#         return self.filtered_sq_fields()
    
#     @cached_property
#     def fsb_binned(self): # average 
#         print('computed fsb_binned for the 1st (and hopefully only) time')
#         return self.get_fsb(wksp=self.w_fsb)
    
#     @cached_property
#     def cls_12_binned(self): # average
#         print('computed cls_12_binned for the 1st (and hopefully only) time')
#         # return self.get_cls_field(np.array([self.field1]), self.effmask, field2=np.array([self.field2]), mask2=self.effmask, wksp=self.w_cls)

#         cls_12_b_sims = [] # FIXME:
#         for n in range(len(self.map1)): # FIXME:
#             cls_12_b_sims.append( self.get_cls_field(np.array([self.field1[n]]), self.rmask, field2=np.array([self.field2[n]]), mask2=self.mask2, wksp=self.w_cls) ) # FIXME:
#         cls_12_b_sims = np.array(cls_12_b_sims) # FIXME:

#         return np.mean(cls_12_b_sims, axis=0) # FIXME:


#     @cached_property
#     def datavector(self):
#         return np.concatenate(self.fsb_binned.flatten(), self.cls_12_binned)
    
#     @cached_property
#     def fsb_unbinned(self): # average
#         print('computed fsb_unbinned for the 1st (and hopefully only) time')
#         return self.get_fsb() 
    
#     @cached_property
#     def cls_11_unbinned(self): # average
#         print('computed cls_11_unbinned for the 1st (and hopefully only) time')

#         cls_11_ub_sims = [] # FIXME:
#         for n in range(len(self.map1)): # FIXME:
#             cls_11_ub_sims.append( self.get_cls_field(np.array([self.field1[n]]), self.effmask) ) # FIXME:
#         cls_11_ub_sims = np.array(cls_11_ub_sims) # FIXME:

#         return np.mean(cls_11_ub_sims, axis=0) # FIXME:
    
#     @cached_property
#     def cls_22_unbinned(self): # average
#         print('computed cls_22_unbinned for the 1st (and hopefully only) time')

#         cls_22_ub_sims = [] # FIXME:
#         for n in range(len(self.map1)): # FIXME:
#             cls_22_ub_sims.append( self.get_cls_field(np.array([self.field2[n]]), self.effmask) ) # FIXME:
#         cls_22_ub_sims = np.array(cls_22_ub_sims) # FIXME:

#         return np.mean(cls_22_ub_sims, axis=0) # FIXME:


#     @cached_property
#     def cls_12_unbinned(self): # average
#         print('computed cls_12_unbinned for the 1st (and hopefully only) time')

#         cls_12_ub_sims = [] # FIXME:
#         for n in range(len(self.map1)): # FIXME:
#             cls_12_ub_sims.append( self.get_cls_field(np.array([self.field1[n]]), self.effmask, field2=np.array([self.field2[n]]), mask2=self.effmask) ) # FIXME:
#         cls_12_ub_sims = np.array(cls_12_ub_sims) # FIXME:

#         return np.mean(cls_12_ub_sims, axis=0) # FIXME:
    
#     # @cached_property
#     # def cls_1F1Bx2(self): 
#     #     return self._get_general_fsb(filters1, filters2)

#     @cached_property
#     def gauss_cov(self):
#         print('computed gauss_cov for the 1st (and hopefully only) time')
#         return self.get_gauss_cov()
    



#     def return_wkspace(self, mask1, mask2, lpb): 

#         """
#         Creates an `NmtWorkspace` object from mask(s) and
#         a binning scheme.

#         Arguments
#         ----------
#         mask1 : array
#             healpy map of size (12*nside**2)
#         mask2 : array
#             a healpy map of the same size as mask1
#         ells_per_bin : int (default = 10)
#             the number of ells in each bin in the binning scheme.

#         Returns
#         ----------
#         A `NmtWorkspace` object corresponding to the
#         given mask(s) and binning scheme.
        
#         """
            
#         fmask1 = nmt.NmtField(mask1, None, spin=0)
#         fmask2 = nmt.NmtField(mask2, None, spin=0)

#         b = nmt.NmtBin.from_lmax_linear(3*self.nside-1, lpb)

#         w12 = nmt.NmtWorkspace()
#         w12.compute_coupling_matrix(fmask1, fmask2, b)

#         return w12 



#     def filtered_sq_fields(self):

#         """
#         Creates NmtField objects for each filtered-squared map.

#         Returns
#         ----------
#         An array containing as many NmtField objects as 
#         there are filters.
        
#         """ 
#         # # already done up there i believe
#         # mask1_bin = self.effmask>0   
#         # map1 = self.map1*mask1_bin

#         f1sq = [] # np.empty((len(self.map1), self.nbands)) # FIXME: 

#         for n in range(len(self.map1)): # FIXME: 

#             alm1 = hp.map2alm(self.map1[n], iter=self.niter) # FIXME: 
            
#             mp_filt_sq = np.array([hp.alm2map(hp.almxfl(alm1, fl), self.nside, lmax=3*self.nside-1)**2 for fl in self.filters]) # FIXME: 
#             step1 = [nmt.NmtField(self.rmask, [m], masked_on_input=False, n_iter=self.niter) for m in mp_filt_sq] 
#             f1sq.append(step1) # FIXME:

#         return np.array(f1sq) # FIXME: 
    
    

#     def get_cls_field(self, field1, mask1, field2=None, mask2=None, wksp=None):

#         """
#         Computes the power spectra of the given field(s).
#         Handles the following cases:
#             - if field1 is a single field, will return 
#             its auto-power spectrum
#             - if field1 is several fields, will return 
#             the corresponding cross-power spectra
#             - if field1 and field2 are both single fields,
#             will return their cross-power spectrum.
#             - if both field1 and field2 are several fields,
#             will return their cross-power spectra. (NOT SURE THIS WORKS)

#         Arguments
#         ----------
#         field1 : `NmtField` object or array of `NmtField` objects
#         mask1 : array
#             a mask of shape corresponding to the individual field1
#         field2 : `NmtField` object or array of `NmtField` objects,
#         optional
#         mask2 : array, optional
#             a mask of shape corresponding to the individual field2
#         wksp : `NmtWorkspace` object, optional
#             a workspace to compute the mode coupling and binning of
#             the power spectra. If not given, it will be assumed that
#             no binning is required, and the fsky correction is applied
#             to the unbinned cls.

#         Returns
#         ----------
#         An array of the cls.
            
#         """

#         if field2 is None:
#             field2 = field1
#             same = True
#         else:
#             same = False
#         if mask2 is None:
#             mask2 = mask1


#         fsky = np.mean(mask1*mask2) # TODO: sus

#         # fsb = self.get_cls_field(self.f1s, self.rmask, field2=np.array([f2]), mask2=self.mask1, wksp=wksp) 


#         if field1.shape[0]>1: # several fields as input

#             if wksp is None: # cross power spectra, unbinned

#                 claa = np.zeros((len(field1), len(field2), 3*self.nside)) # ?

#                 if same is True:
#                     for n in range(len(field1)):
#                         for m in range(n, len(field2)):
#                             cross = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
#                             claa[n, m] = cross; claa[m, n] = cross
#                 else:
#                     for n in range(len(field1)):
#                         for m in range(len(field2)): # if field2!=field1, should not start from n?
#                             claa[n, m] = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
#                             # TODO: are we using this bit in the covariance? why /fsky?
                                    
#             else: 

#                 if same is True: # auto power spectra, binned
#                     claa = np.array([wksp.decouple_cell(nmt.compute_coupled_cell(fi, fi))[0] for fi in field1])

#                 else: # cross power spectra, binned
#                     claa = np.zeros((len(field1), len(field2), self.b)) 
#                     for n in range(len(field1)):
#                         for m in range(len(field2)): # if field2!=field1, should not start from n?
#                             claa[n, m] = wksp.decouple_cell(nmt.compute_coupled_cell(field1[n], field2[m]))[0]
                    
#             return claa.squeeze() 
        
#         else: # one field as input, inside a np.array

#             if wksp is None: # auto power spectra, unbinned
#                 clbb = nmt.compute_coupled_cell(field1[0], field2[0])[0] / fsky

#             else: # auto power spectra, binned
#                 clbb = wksp.decouple_cell(nmt.compute_coupled_cell(field1[0], field2[0]))[0]

#             return clbb



#     def get_fsb(self, wksp=None):

#         """
#         Computes the FSB of the class field(s)
#         for a set of filters. if map1 and map2 are different, 
#         only map1 will be filtered.
        
#         Arguments
#         ----------
        
#         wksp12 : `NmtWorkspace` object, optional
#             a workspace to compute the mode coupling and binning of
#             the power spectra. if not given, it will be assumed that
#             no binning is required, and the fsky correction is applied
#             to the unbinned FSB.
        
#         Returns
#         ----------
#         An array of the FSBs.
            
#         """

#         fsb_sims = [] # FIXME:
#         for n in range(len(self.f1s)): # FIXME:
#             fsb_sims.append( self.get_cls_field(self.f1s[n], self.rmask, field2=np.array([self.field2[n]]), mask2=self.mask2, wksp=wksp) ) # FIXME:
#         fsb_sims = np.array(fsb_sims) # FIXME:

#         return np.mean(fsb_sims, axis=0) # FIXME:
        


#     def get_gauss_cov(self, insquares=True):

#         """
#         Computes the gaussian-limit approximation
#         of the FSB+Cl covariance. 
        
#         Returns
#         ----------
#         A FSB+Cl covariance matrix.
#         Should be of size ((nbands+1)*ndatapoints, (nbands+1)*ndatapoints)
#         if the same binning is used for both FSBs and Cls.
            
#         """

#         # self.cls_1sq1sq_unbinned = self.get_cls_field(self.f1s, self.effmask) # FIXME: changed to below
#         cls_1sq1sq_unbinned = []
#         for n in range(len(self.map1)):
#             step1 = self.get_cls_field(self.f1s[n], self.effmask)
#             cls_1sq1sq_unbinned.append(step1)
#         self.cls_1sq1sq_unbinned = np.mean(np.array(cls_1sq1sq_unbinned), axis=0)

#         fmask_r = nmt.NmtField(self.rmask, None, spin=0)
#         fmask = nmt.NmtField(self.effmask, None, spin=0)
#         cw = nmt.NmtCovarianceWorkspace()
#         cw.compute_coupling_coefficients(fmask, fmask_r, fmask, fmask_r)

#         gauss_cov = np.zeros((self.nbands+1, self.nbands+1, self.b, self.b)) 

#         fsbs_cls = np.zeros((self.nbands+1, self.fsb_unbinned.shape[1]))
#         fsbs_cls[:self.nbands] = self.fsb_unbinned
#         fsbs_cls[-1] = self.cls_22_unbinned # last col is for cls


#         for n in range(len(fsbs_cls)):
#             for m in range(n, len(fsbs_cls)):

#                 clad = fsbs_cls[n]
#                 clbc = fsbs_cls[m]
#                 clbd = fsbs_cls[-1] # cls
                
#                 if m==self.nbands: # when computing cross FSB-Cl cov
#                     clac = clad 
#                     clbc = clbd
#                 else: # when computing cross FSB-FSB cov
#                     clac = self.cls_1sq1sq_unbinned[n, m] 

#                 if n==self.nbands: # when computing auto Cl cov
#                     clac = fsbs_cls[-1] # cls 
#                     clad = fsbs_cls[-1] # cls
#                     clbc = fsbs_cls[-1] # cls
                    
#                 covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [clac], [clad], [clbc], [clbd], self.w_fsb)
#                 gauss_cov[n, m] = covij_fsb
#                 gauss_cov[m, n] = covij_fsb
        
#         self.gauss_cov = gauss_cov

#         if insquares==False:
#             return _reduce2(self.gauss_cov)
#         else:
#             return self.gauss_cov



#     def _get_n222_term(self, filter):

#         """
#         Computes the N222 term for each FSB.

#         Arguments
#         ----------
#         cls : array 
#             an array containing the power spectrum of the
#             original map (before filtering/squaring)
#         filter : array
#             a single binary array of shape (3*nside, )

#         Returns
#         ----------
#         The not yet mask-corrected N222 term in a
#         block matrix of size (len(cls), len(cls)).
#         """

#         # if self.cls_12_unbinned is None:
#         #     self.cls_12_unbinned = self.get_cls_field(np.array([self.field2]), self.effmask) # obvs wrong
#         # if self.cls_11_unbinned is None:
#         #     self.cls_11_unbinned = self.get_cls_field(np.array([self.field2]), self.effmask) # obvs wrong

#         cls_11_f = self.cls_11_unbinned*filter 
#         cls_12_f = self.cls_12_unbinned*filter

#         # wondering if this can be defined earlier
#         lmax = len(cls_11_f)-1 
#         beam = np.ones_like(cls_11_f)
#         ls = np.arange(lmax+1)
#         bin = nmt.NmtBin.from_lmax_linear(lmax, 1)

#         w = nmt.nmtlib.comp_coupling_matrix(0, 0, lmax, lmax, 0, 0, 0, 0, beam, beam, cls_11_f, bin.bin, 0, -1, -1, -1)
#         sum_l1 = nmt.nmtlib.get_mcm(w, (lmax+1)**2).reshape([lmax+1, lmax+1])
#         sum_l1 /= (2*ls+1)[None, :]
#         cov = sum_l1 * 4 * np.outer(cls_12_f, cls_12_f)
#         nmt.nmtlib.workspace_free(w)

#         return cov

#     def get_n222_cov(self):

#         """
#         Computes the mask-corrected, binned N222
#         terms for all FSBs.

#         Returns
#         ----------
#         The binned mask-corrected N222 term in an
#         array of dimensions (nbands+1, nbands+1, nbins).
#         """

#         n222 = np.zeros((len(self.filters)+1, len(self.filters)+1, self.b, self.b))

#         for nb in range(self.nbands):
#             n222ub = self._get_n222_term(self.filters[nb])
#             n222_bin1 = np.array([self.bb.bin_cell(row) for row in n222ub])
#             n222_final = np.array([self.bb.bin_cell(col) for col in n222_bin1.T]).T
#             n222[nb, nb] += n222_final 

#         return n222/self.fsky_cls


#     def _get_general_fsb(self, filters1, filters2): # maybe use fsq_fields

#         """
#         Computes the generalized FSB for 
#         two sets of filters.

#         Arguments
#         ----------
#         filters1 : array
#             a binary array of shape (nbands1, 3*nside)
#         filters2 : array 
#             a binary array of shape (nbands2, 3*nside)
        
#         Returns
#         ----------
#         The binned generalized FSBs in an array
#         of shape (nbands1, nbands2, 3*nside).
#         """

#         alms1 = hp.map2alm(self.map1) # already computed earlier i believe


#         maps_F = np.array([hp.alm2map(hp.almxfl(alms1, fl), self.nside) for fl in filters1])
#         maps_B = np.array([hp.alm2map(hp.almxfl(alms1, bl), self.nside) for bl in filters2])

#         self.cls_1F1Bx2 = np.zeros((len(filters1), len(filters2), len(filters1[0])))

#         for f in range(len(filters1)):
#             # print(f)
#             for b in range(len(filters2)):
#                 # print('\t', b)
#                 map_FB = maps_F[f]*maps_B[b]
#                 self.cls_1F1Bx2[f, b] = hp.anafast(map_FB, self.map2) / (self.fsky_fsb * self.ells_per_bin) # could use get_cls_fields? or too much work?
#                 # last bit above: must divide by normalisation factor lpb

#         # return self.cls_mFBxm
        
#         if self.twofields is True:
#             alms2 = hp.map2alm(self.map2)
#             maps_B = np.array([hp.alm2map(hp.almxfl(alms2, bl), self.nside) for bl in filters2])

#             self.cls_1F2Bx2 = np.zeros((len(filters1), len(filters2), len(filters1[0])))

#             for f in range(len(filters1)):
#                 # print(f)
#                 for b in range(len(filters2)):
#                     # print('\t', b)
#                     map_FB = maps_F[f]*maps_B[b]
#                     self.cls_1F2Bx2[f, b] = hp.anafast(map_FB, self.map2) / (self.fsky_fsb * self.ells_per_bin)
#                     # must divide by normalisation factor lpb

#             # return self.cls_mFBxm



#     def twonickels(self, cls, genfsb, filters1, filters2):
#         """binning and multiplying the arguments"""

#         # bin general fsb
#         binned_genFSB = np.zeros((len(filters1), len(filters2), self.b)) 
#         for n in range(len(filters1)):
#             for b in range(len(filters2)):
#                 binned_genFSB[n, b] = self.bb.bin_cell(genfsb[n, b])
#         # bin filtered cls
#         cls_f = np.array([cls*fl for fl in filters1])
#         un = np.array([self.bb.bin_cell(e) for e in cls_f]) 
#         deux = np.array([np.repeat([e], len(filters2), axis=0).T for e in un])

#         fsb_gen = np.zeros((len(filters1)+1, len(filters1)+1, len(filters2), self.b))
#         cl_filter = np.zeros((len(filters1)+1, len(filters1)+1, len(filters2), self.b))

#         fsb_gen[len(filters1), :len(filters1)] = binned_genFSB
#         cl_filter[len(filters1), :len(filters1)] = deux

#         return cl_filter*fsb_gen



#     def get_n32_cov(self, filters1, filters2):

#         """
#         Computes the mask-corrected, binned N32
#         terms for all FSBs x Cls.

#         Arguments
#         ----------
#         filters1 : array
#             a binary array of shape (nbands1, 3*nside)
#         filters2 : array 
#             a binary array of shape (nbands2, 3*nside)

#         Returns
#         ----------
#         The binned mask-corrected N32 term in an
#         array of dimensions (nbands+1, nbands+1, nbins).
#         (if we assume filters2 to be the equivalent of 
#         the binning scheme.)
#         """

#         if self.cls_1F1Bx2 is None:
#             self._get_general_fsb(filters1, filters2)

#         n32_term1 = self.twonickels(self.cls_11_unbinned, self.cls_1F1Bx2, filters1, filters2)

#         if self.twofields is False:
#             n32 = 4*n32_term1 / ( np.array([(2*self.bb.get_effective_ells()+1)]).T )
#         else:
#             n32_term2 = self.twonickels(self.cls_12_unbinned, self.cls_1F2Bx2, filters1, filters2)
#             n32 = 2*(n32_term1 + n32_term2) / ( np.array([(2*self.bb.get_effective_ells()+1)]).T )

#         for i in range(self.nbands): # make it symmetric
#             n32[i, self.nbands] = n32[self.nbands, i].T

#         return n32 / self.fsky_cls
    
#     def get_full_cov(self, insquares=False):

#         """
#         Adds the mask-corrected, binned gaussian-
#         limit covariance and its main mask-corrected, 
#         binned non-gaussian contributions. 

#         Returns
#         ----------
#         An array of dimensions ((nbands+1)*nbins, (nbands+1)*nbins),
#         or if insquares set to True, an array of dimensions 
#         ((nbands+1), (nbands+1), nbins, nbins).
#         """
#         # if self.gauss_cov is None:
#         #     self.gauss_cov = self.get_gauss_cov()
#         self.full_cov_large = self.gauss_cov + self.get_n222_cov() # + self.get_n32_cov(self.filters, self.bins) # FIXME: n32 does not matter for studying linterm + a pain to compute over sims
#         self.full_cov = _reduce2(self.full_cov_large)
        
#         if insquares is True:
#             return self.full_cov_large
#         else:
#             return self.full_cov
    






class FSBmap2():

    def __init__(self, map1, mask1, filters, map2=None, mask2=None, ells_per_bin=10, niter=3): # , rmask=None

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
        

        # set the effective mask to the intersection of the 2 masks
        self.effmask = (self.mask1*self.mask2) > 0 # make binary
        self.rmask = self.effmask # TODO: construct options for remasking
        # remasking fields appropriately + need to make sure fields within new mask is 0
        self.map1 = (self.map1-np.mean(self.map1[self.effmask==1]))*self.effmask
        self.map2 = np.array([(m-np.mean(m[self.effmask==1]))*self.effmask for m in self.map2]) # FIXME: changed


        self.filters = filters
        self.ells_per_bin = ells_per_bin
        self.niter = niter

        # binning
        self.bb = nmt.NmtBin.from_lmax_linear(3*self.nside-1, self.ells_per_bin)
        self.b = len(self.bb.get_effective_ells())
        self.bins = get_filters(self.b, self.nside) # filters corresponding to bins (for generalized fsb)
        
        
        self.w_fsb = self.return_wkspace(self.rmask, self.effmask, self.ells_per_bin)
        self.w_cls = self.return_wkspace(self.effmask, self.effmask, self.ells_per_bin)

        self.fsky_fsb = np.mean(self.rmask*self.effmask)
        self.fsky_cls = np.mean(self.effmask*self.effmask)

        self.nbands = len(self.filters)

        self.field1 = nmt.NmtField(self.effmask, [self.map1], masked_on_input=False, n_iter=self.niter)
        if map2 is None:
            self.field2 = self.field1
        else:
            self.field2 = np.array([nmt.NmtField(self.effmask, [m], masked_on_input=False, n_iter=self.niter) for m in self.map2]) # FIXME: changed

        # self.cls_11_binned = self.get_cls_field(np.array([self.field1]), self.effmask, wksp=self.w_cls)
        # self.cls_22_binned = self.get_cls_field(np.array([self.field2]), self.effmask, wksp=self.w_cls)
        
        # for None statements
        # self.gauss_cov = None
        # self.cls_11_unbinned = None # by 1, we do NOT mean the filter squared field, simply the original f1
        # self.cls_22_unbinned = None
        # self.cls_12_unbinned = None
        self.cls_1F1Bx2 = None # TODO: make default usage with self.filters and self.binfilters

    
    
    @cached_property
    def f1s(self):
        print('computed f1s for the 1st (and hopefully only) time')
        return self.filtered_sq_fields()
    
    @cached_property
    def fsb_binned(self): # TODO: average
        print('computed fsb_binned for the 1st (and hopefully only) time')
        return self.get_fsb(wksp=self.w_fsb)
    
    @cached_property
    def cls_12_binned(self): # TODO: average
        print('computed cls_12_binned for the 1st (and hopefully only) time')
        # return self.get_cls_field(np.array([self.field1]), self.effmask, field2=np.array([self.field2]), mask2=self.effmask, wksp=self.w_cls)
        cls_12_b_sims = [] # FIXME:
        for n in range(len(self.map2)): # FIXME:
            cls_12_b_sims.append( self.get_cls_field(np.array([self.field1]), self.rmask, field2=np.array([self.field2[n]]), mask2=self.mask2, wksp=self.w_cls) ) # FIXME:
        cls_12_b_sims = np.array(cls_12_b_sims) # FIXME:
        return np.mean(cls_12_b_sims, axis=0) # FIXME:

    @cached_property
    def datavector(self):
        return np.concatenate((self.fsb_binned.flatten(), self.cls_12_binned))
    
    @cached_property
    def fsb_unbinned(self): # TODO: average
        print('computed fsb_unbinned for the 1st (and hopefully only) time')
        return self.get_fsb() 
    
    @cached_property
    def cls_11_unbinned(self): 
        print('computed cls_11_unbinned for the 1st (and hopefully only) time')
        return self.get_cls_field(np.array([self.field1]), self.effmask)

    @cached_property
    def cls_22_unbinned(self): # TODO: average
        print('computed cls_22_unbinned for the 1st (and hopefully only) time')
        # return self.get_cls_field(np.array([self.field2]), self.effmask)

        cls_22_ub_sims = [] # FIXME:
        for n in range(len(self.map2)): # FIXME:
            cls_22_ub_sims.append( self.get_cls_field(np.array([self.field2[n]]), self.effmask) ) # FIXME:
        cls_22_ub_sims = np.array(cls_22_ub_sims) # FIXME:

        return np.mean(cls_22_ub_sims, axis=0) # FIXME:
    
    @cached_property
    def cls_12_unbinned(self): # TODO: average
        print('computed cls_12_unbinned for the 1st (and hopefully only) time')
        # return self.get_cls_field(np.array([self.field1]), self.effmask, field2=np.array([self.field2]), mask2=self.effmask)

        cls_12_ub_sims = [] # FIXME:
        for n in range(len(self.map2)): # FIXME:
            cls_12_ub_sims.append( self.get_cls_field(np.array([self.field1]), self.effmask, field2=np.array([self.field2[n]]), mask2=self.effmask) ) # FIXME:
        cls_12_ub_sims = np.array(cls_12_ub_sims) # FIXME:

        return np.mean(cls_12_ub_sims, axis=0) # FIXME:
    
    # @cached_property
    # def cls_1F1Bx2(self): 
    #     return self._get_general_fsb(filters1, filters2)

    @cached_property
    def gauss_cov(self):
        print('computed gauss_cov for the 1st (and hopefully only) time')
        return self.get_gauss_cov()
    

    # NEW!
    # ----------------------------------------------------------------------------------------------------
    @cached_property
    def cls_1sq1sq_unbinned(self): 
        print('computed cls_1sq1sq_unbinned for the 1st (and hopefully only) time')
        return self.get_cls_field(self.f1s, self.effmask)

    @cached_property
    def fsb_unbinned_pure(self): 
        print('computed fsb_unbinned_pure for the 1st (and hopefully only) time')
        return self.get_cls_field(self.f1s, self.effmask, field2=np.array([self.field1]), mask2=self.effmask) 
    
    @cached_property
    def fsb_binned_pure(self): 
        print('computed fsb_binned_pure for the 1st (and hopefully only) time')
        return self.get_cls_field(self.f1s, self.effmask, field2=np.array([self.field1]), mask2=self.effmask, wksp=self.w_fsb) 
    
    @cached_property
    def cls_11_binned(self): 
        print('computed cls_11_binned for the 1st (and hopefully only) time')
        return self.get_cls_field(np.array([self.field1]), self.effmask, wksp=self.w_cls)
    
    @cached_property
    def master_datavector(self):
        return np.concatenate((self.fsb_binned_pure.flatten(), self.cls_11_binned, self.fsb_binned.flatten(), self.cls_12_binned))
    
    @cached_property
    def cls_datavector(self):
        return np.concatenate((self.cls_11_binned, self.cls_12_binned))
    # ----------------------------------------------------------------------------------------------------


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

        b = nmt.NmtBin.from_lmax_linear(3*self.nside-1, lpb)

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
        
        mp_filt_sq = np.array([hp.alm2map(hp.almxfl(alm1, fl), self.nside, lmax=3*self.nside-1)**2 for fl in self.filters])   
        
        f1sq = [nmt.NmtField(self.rmask, [m], masked_on_input=False, n_iter=self.niter) for m in mp_filt_sq] 
        
        return np.array(f1sq)
    
    

    def get_cls_field(self, field1, mask1, field2=None, mask2=None, wksp=None):

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

            TODO: must add an option equivalent to option 2, 
            but where field1 is one single field and field2 an array

        Arguments
        ----------
        field1 : `NmtField` object or array of `NmtField` objects
        mask1 : array
            a mask of shape corresponding to the individual field1
        field2 : `NmtField` object or array of `NmtField` objects,
        optional
        mask2 : array, optional
            a mask of shape corresponding to the individual field2
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
        if mask2 is None:
            mask2 = mask1


        fsky = np.mean(mask1*mask2) # TODO: sus

        if field1.shape[0]>1: # several fields as input

            if wksp is None: # cross power spectra, unbinned

                claa = np.zeros((len(field1), len(field2), 3*self.nside)) # ?

                if same is True:
                    for n in range(len(field1)):
                        for m in range(n, len(field2)):
                            cross = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
                            claa[n, m] = cross; claa[m, n] = cross
                else:
                    for n in range(len(field1)):
                        for m in range(len(field2)): # if field2!=field1, should not start from n?
                            claa[n, m] = nmt.compute_coupled_cell(field1[n], field2[m])[0] / fsky
                            # TODO: are we using this bit in the covariance? why /fsky?
                                      
            else: 

                if same is True: # auto power spectra, binned
                    claa = np.array([wksp.decouple_cell(nmt.compute_coupled_cell(fi, fi))[0] for fi in field1])

                else: # cross power spectra, binned
                    claa = np.zeros((len(field1), len(field2), self.b)) 
                    for n in range(len(field1)):
                        for m in range(len(field2)): # if field2!=field1, should not start from n?
                            claa[n, m] = wksp.decouple_cell(nmt.compute_coupled_cell(field1[n], field2[m]))[0]
                    
            return claa.squeeze() 
        
        else: # one field as input, inside a np.array

            if wksp is None: # auto power spectra, unbinned
                clbb = nmt.compute_coupled_cell(field1[0], field2[0])[0] / fsky

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

        fsb_sims = [] # FIXME:
        for n in range(len(self.field2)): # FIXME:
            fsb_sims.append( self.get_cls_field(self.f1s, self.rmask, field2=np.array([self.field2[n]]), mask2=self.mask2, wksp=wksp) ) # FIXME:
        fsb_sims = np.array(fsb_sims) # FIXME:

        return np.mean(fsb_sims, axis=0) # FIXME:

        # return self.get_cls_field(self.f1s, self.rmask, field2=np.array([self.field2]), mask2=self.mask2, wksp=wksp) 


    def get_gauss_cov(self, insquares=True):
        
        fmask_r = nmt.NmtField(self.rmask, None, spin=0)
        fmask = nmt.NmtField(self.effmask, None, spin=0)
        cw = nmt.NmtCovarianceWorkspace()
        cw.compute_coupling_coefficients(fmask, fmask_r, fmask, fmask_r)

        gauss_cov = np.zeros((self.nbands+1, self.nbands+1, self.b, self.b)) 

        for n in range(self.nbands):

            # this is for the fsb-cls cross covariance
            cla1b1 = self.fsb_unbinned_pure[n]
            cla1b2 = self.fsb_unbinned[n]
            cla2b1 = self.cls_12_unbinned
            cla2b2 = self.cls_22_unbinned

            covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
            gauss_cov[n, -1] = covij_fsb
            gauss_cov[-1, n] = covij_fsb

            for m in range(n, self.nbands):

                # this is for the fsb-fsb covariance
                cla1b1 = self.cls_1sq1sq_unbinned[n,m]
                cla2b1 = self.fsb_unbinned[m]

                covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
                gauss_cov[n, m] = covij_fsb
                gauss_cov[m, n] = covij_fsb

        gauss_cov[-1, -1] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [self.cls_11_unbinned], [self.cls_12_unbinned], [self.cls_12_unbinned], [self.cls_22_unbinned], self.w_cls)
        
        self.gauss_cov = gauss_cov

        if insquares==False:
            return _reduce2(self.gauss_cov)
        else:
            return self.gauss_cov
        
    

    def get_gauss_cov_cls(self, insquares=True):

        fmask = nmt.NmtField(self.effmask, None, spin=0)
        cw = nmt.NmtCovarianceWorkspace()
        cw.compute_coupling_coefficients(fmask, fmask, fmask, fmask)

        gauss_cov_cls = np.zeros((2, 2, self.b, self.b)) 

        cls = np.zeros((3, self.cls_11_unbinned.shape[0]))
        cls[0] = self.cls_11_unbinned # gg
        cls[1] = self.cls_12_unbinned # gk
        cls[2] = self.cls_22_unbinned # kk

        gauss_cov_cls[0,0] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cls[0]], [cls[0]], [cls[0]], [cls[0]], self.w_cls) # gg, gg
        gauss_cov_cls[1,1] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cls[0]], [cls[1]], [cls[1]], [cls[2]], self.w_cls) # gk, gk
        gauss_cov_cls[0,1] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cls[0]], [cls[1]], [cls[0]], [cls[1]], self.w_cls) # gg, gk
        gauss_cov_cls[1,0] = gauss_cov_cls[0,1] # symmetric
        
        self.cov_cls = gauss_cov_cls

        if insquares==False:
            return _reduce2(self.cov_cls)
        else:
            return self.cov_cls
        


    def get_gauss_cov_pure(self, insquares=True):
    
        fmask_r = nmt.NmtField(self.rmask, None, spin=0)
        fmask = nmt.NmtField(self.effmask, None, spin=0)
        cw = nmt.NmtCovarianceWorkspace()
        cw.compute_coupling_coefficients(fmask, fmask_r, fmask, fmask_r)

        gauss_cov = np.zeros((self.nbands+1, self.nbands+1, self.b, self.b)) 

        for n in range(self.nbands):

            cla1b1 = self.fsb_unbinned_pure[n]
            cla1b2 = self.fsb_unbinned_pure[n]
            cla2b1 = self.cls_11_unbinned
            cla2b2 = self.cls_11_unbinned

            covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
            gauss_cov[n, -1] = covij_fsb
            gauss_cov[-1, n] = covij_fsb

            for m in range(n, self.nbands):

                cla1b1 = self.cls_1sq1sq_unbinned[n,m]
                cla2b1 = self.fsb_unbinned_pure[m]

                covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
                gauss_cov[n, m] = covij_fsb
                gauss_cov[m, n] = covij_fsb

        gauss_cov[-1, -1] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [self.cls_11_unbinned], [self.cls_11_unbinned], [self.cls_11_unbinned], [self.cls_11_unbinned], self.w_cls)
        
        self.gauss_cov_pure = gauss_cov

        if insquares==False:
            return _reduce2(self.gauss_cov_pure)
        else:
            return self.gauss_cov_pure
        


    def get_master_gauss_cov(self, insquares=True):

            # self.cls_1sq1sq_unbinned

            fmask_r = nmt.NmtField(self.rmask, None, spin=0)
            fmask = nmt.NmtField(self.effmask, None, spin=0)
            cw = nmt.NmtCovarianceWorkspace()
            cw.compute_coupling_coefficients(fmask, fmask_r, fmask, fmask_r)

            master_cov = np.zeros((2*(self.nbands+1), 2*(self.nbands+1), self.b, self.b)) 
            # let's not use this for now (we check if indices for mixed cov are correct first)
            oldggg = self.get_gauss_cov_pure() # call gauss cov in squares
            oldggk = self.gauss_cov # call gauss cov in squares
            master_cov[:self.nbands+1, :self.nbands+1] = oldggg
            master_cov[self.nbands+1:, self.nbands+1:] = oldggk
            # now we need master_cov[self.nbands+1:, :self.nbands+1]

            fsbs_cls_mixed = np.zeros((self.nbands+1, self.fsb_unbinned.shape[1]))
            fsbs_cls_mixed[:self.nbands] = self.fsb_unbinned
            fsbs_cls_mixed[-1] = self.cls_12_unbinned
            fsbs_cls_pure = np.zeros_like(fsbs_cls_mixed)
            fsbs_cls_pure[:self.nbands] = self.fsb_unbinned_pure
            fsbs_cls_pure[-1] = self.cls_11_unbinned


            for n in range(self.nbands):

                cla1b1 = fsbs_cls_pure[n] # for both cases
                cla1b2 = fsbs_cls_pure[-1]
                cla2b1 = fsbs_cls_mixed[n]
                cla2b2 = fsbs_cls_mixed[-1]

                # ggg, gk (along axis 1, horizontal)
                covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
                master_cov[2*(self.nbands+1)-1, n] = covij_fsb
                # ggk, gg (along axis 0, vertical)
                covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b1], [cla2b2], [cla2b2], self.w_fsb)
                master_cov[(self.nbands+1) + n, (self.nbands)] = covij_fsb

                for m in range(self.nbands):

                    cla1b1 = self.cls_1sq1sq_unbinned[n,m]
                    cla1b2 = fsbs_cls_pure[n]
                    cla2b1 = fsbs_cls_mixed[m]
                    # cla2b2 = fsbs_cls_mixed[-1] # same as above
        
                    covij_fsb = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [cla1b1], [cla1b2], [cla2b1], [cla2b2], self.w_fsb)
                    master_cov[(self.nbands+1) + n, m] = covij_fsb
                    # master_cov[m, (self.nbands+1) + n] = covij_fsb # not symmetric bcs Phi_ggg and Phi_ggk
                        
            master_cov[2*(n+1)+1, n+1] = nmt.gaussian_covariance(cw, 0, 0, 0, 0, [self.cls_11_unbinned], [self.cls_12_unbinned], [self.cls_11_unbinned], [self.cls_12_unbinned], self.w_cls)
            # add symmetry in upper right corner // there's probably an easier way to do this, using ::-1?
            transposed = np.transpose(master_cov[self.nbands+1:, :self.nbands+1], (1, 0, 2, 3))
            # transposed = np.zeros_like(master_cov[self.nbands+1:, :self.nbands+1])
            # for i in range(self.nbands+1):
            #     for j in range(self.nbands+1):
            #         transposed[i,j] = master_cov[self.nbands+1:, :self.nbands+1][j,i]
            master_cov[:self.nbands+1, self.nbands+1:] = transposed # master_cov[self.nbands+1:, :self.nbands+1]


            self.master_cov = master_cov

            if insquares==False:
                return _reduce2(self.master_cov)
            else:
                return self.master_cov



    def _get_n222_term(self, filter):

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

        # if self.cls_12_unbinned is None:
        #     self.cls_12_unbinned = self.get_cls_field(np.array([self.field2]), self.effmask) # obvs wrong
        # if self.cls_11_unbinned is None:
        #     self.cls_11_unbinned = self.get_cls_field(np.array([self.field2]), self.effmask) # obvs wrong

        cls_11_f = self.cls_11_unbinned*filter 
        cls_12_f = self.cls_12_unbinned*filter

        # wondering if this can be defined earlier
        lmax = len(cls_11_f)-1 
        beam = np.ones_like(cls_11_f)
        ls = np.arange(lmax+1)
        bin = nmt.NmtBin.from_lmax_linear(lmax, 1)

        w = nmt.nmtlib.comp_coupling_matrix(0, 0, lmax, lmax, 0, 0, 0, 0, 0, 0.0,
                                            beam, beam, cls_11_f, bin.bin, 0, -1, -1, -1)
        sum_l1 = nmt.nmtlib.get_mcm(w, (lmax+1)**2).reshape([lmax+1, lmax+1])
        sum_l1 /= (2*ls+1)[None, :]
        cov = sum_l1 * 4 * np.outer(cls_12_f, cls_12_f)
        nmt.nmtlib.workspace_free(w)

        return cov

    def get_n222_cov(self):

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
            n222ub = self._get_n222_term(self.filters[nb])
            n222_bin1 = np.array([self.bb.bin_cell(row) for row in n222ub])
            n222_final = np.array([self.bb.bin_cell(col) for col in n222_bin1.T]).T
            n222[nb, nb] += n222_final 

        return n222/self.fsky_cls


    # def _get_general_fsb(self, filters1, filters2): # maybe use fsq_fields

    #     """
    #     Computes the generalized FSB for 
    #     two sets of filters.

    #     Arguments
    #     ----------
    #     filters1 : array
    #         a binary array of shape (nbands1, 3*nside)
    #     filters2 : array 
    #         a binary array of shape (nbands2, 3*nside)
        
    #     Returns
    #     ----------
    #     The binned generalized FSBs in an array
    #     of shape (nbands1, nbands2, 3*nside).
    #     """

    #     alms1 = hp.map2alm(self.map1) # already computed earlier i believe


    #     maps_F = np.array([hp.alm2map(hp.almxfl(alms1, fl), self.nside) for fl in filters1])
    #     maps_B = np.array([hp.alm2map(hp.almxfl(alms1, bl), self.nside) for bl in filters2])

    #     self.cls_1F1Bx2 = np.zeros((len(filters1), len(filters2), len(filters1[0])))

    #     for f in range(len(filters1)):
    #         # print(f)
    #         for b in range(len(filters2)):
    #             # print('\t', b)
    #             map_FB = maps_F[f]*maps_B[b]
    #             self.cls_1F1Bx2[f, b] = hp.anafast(map_FB, self.map2) / (self.fsky_fsb * self.ells_per_bin) # could use get_cls_fields? or too much work?
    #             # last bit above: must divide by normalisation factor lpb

    #     # return self.cls_mFBxm
        
    #     if self.twofields is True:
    #         alms2 = hp.map2alm(self.map2)
    #         maps_B = np.array([hp.alm2map(hp.almxfl(alms2, bl), self.nside) for bl in filters2])

    #         self.cls_1F2Bx2 = np.zeros((len(filters1), len(filters2), len(filters1[0])))

    #         for f in range(len(filters1)):
    #             # print(f)
    #             for b in range(len(filters2)):
    #                 # print('\t', b)
    #                 map_FB = maps_F[f]*maps_B[b]
    #                 self.cls_1F2Bx2[f, b] = hp.anafast(map_FB, self.map2) / (self.fsky_fsb * self.ells_per_bin)
    #                 # must divide by normalisation factor lpb

    #         # return self.cls_mFBxm



    def twonickels(self, cls, genfsb, filters1, filters2):
        """binning and multiplying the arguments"""

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



    # def get_n32_cov(self, filters1, filters2):

    #     """
    #     Computes the mask-corrected, binned N32
    #     terms for all FSBs x Cls.

    #     Arguments
    #     ----------
    #     filters1 : array
    #         a binary array of shape (nbands1, 3*nside)
    #     filters2 : array 
    #         a binary array of shape (nbands2, 3*nside)

    #     Returns
    #     ----------
    #     The binned mask-corrected N32 term in an
    #     array of dimensions (nbands+1, nbands+1, nbins).
    #     (if we assume filters2 to be the equivalent of 
    #     the binning scheme.)
    #     """

    #     if self.cls_1F1Bx2 is None:
    #         self._get_general_fsb(filters1, filters2)

    #     n32_term1 = self.twonickels(self.cls_11_unbinned, self.cls_1F1Bx2, filters1, filters2)

    #     if self.twofields is False:
    #         n32 = 4*n32_term1 / ( np.array([(2*self.bb.get_effective_ells()+1)]).T )
    #     else:
    #         n32_term2 = self.twonickels(self.cls_12_unbinned, self.cls_1F2Bx2, filters1, filters2)
    #         n32 = 2*(n32_term1 + n32_term2) / ( np.array([(2*self.bb.get_effective_ells()+1)]).T )

    #     for i in range(self.nbands): # make it symmetric
    #         n32[i, self.nbands] = n32[self.nbands, i].T

    #     return n32 / self.fsky_cls
    
    
    def get_full_cov(self, insquares=False, n32=True):

        """
        Adds the mask-corrected, binned gaussian-
        limit covariance and its main mask-corrected, 
        binned non-gaussian contributions. 

        Returns
        ----------
        An array of dimensions ((nbands+1)*nbins, (nbands+1)*nbins),
        or if insquares set to True, an array of dimensions 
        ((nbands+1), (nbands+1), nbins, nbins).
        """
        # if self.gauss_cov is None:
        #     self.gauss_cov = self.get_gauss_cov()

        if n32 is True:
            self.full_cov_large = self.gauss_cov + self.get_n222_cov() # + self.get_n32_cov(self.filters, self.bins) # TODO: add n32 back
        else: # bypass n32 altogether if not needed
            self.full_cov_large = self.gauss_cov + self.get_n222_cov()
        
        if insquares is True:
            return self.full_cov_large
        else:
            return _reduce2(self.full_cov_large)












