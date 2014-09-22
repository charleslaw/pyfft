# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 00:09:26 2010

@author: Charles Law
"""

import math

def fft(fin, inverse):
    nfft = len(fin)
    twiddles, factors = fft_alloc(nfft, inverse)
    
    fout = []
    for i in xrange(nfft):
        fout.append((0, 0))
        
    fout_ind_start = 0
    fin_ind_start = 0
    in_stride = 1
    fft_work(fout, fout_ind_start, fin, fin_ind_start, 1, in_stride, factors,
            twiddles, nfft)

    return fout




def fft_work(fout, fout_ind, f, f_ind, fstride, in_stride, factors,
            twiddles, nfft_orig):
    p = factors[0][0] # the radix
    m = factors[0][1] # stage's fft length/p
    factors = factors[1:]
    
    fout_beg = fout_ind
    fout_end = fout_ind + p*m

    if m == 1 :
        dowhile_if = 1
        while ( dowhile_if ):
            fout[fout_ind] = f[f_ind]
            f_ind = f_ind + fstride*in_stride
            fout_ind = fout_ind + 1
            if fout_ind == fout_end:
                dowhile_if = 0
    else:
        dowhile_if = 1
        while ( dowhile_if ):
            # recursive call:
            # DFT of size m*p performed by doing
            # p instances of smaller DFTs of size m, 
            # each one takes a decimated version of the input
            fft_work(fout, fout_ind , f, f_ind, fstride*p, in_stride,
                     factors, twiddles, nfft_orig)
            f_ind = f_ind + fstride*in_stride
            
            #}while( (fout += m) != fout_end )
            fout_ind = fout_ind + m
            if ( fout_ind == fout_end ):
                dowhile_if = 0

    fout_ind = fout_beg
    
    # recombine the p smaller DFTs 
    if p == 2:
        fft_bfly2(fout, fout_ind, fstride, twiddles, m)
    elif p == 3:
        fft_bfly3(fout, fout_ind, fstride, twiddles, m)
    else:
        fft_bfly_generic(fout, fout_ind, fstride, twiddles, m, p, nfft_orig)
    
    return fout




def fft_bfly2(fout, fout_ind, fstride, twiddles, m):
    tw1_ind = 0
    fout2_ind = fout_ind + m
    
    dowhile_if = 1
    while(dowhile_if):

        t = _mult ( fout[fout2_ind], twiddles[tw1_ind] )
        tw1_ind = tw1_ind + fstride

        fout[fout2_ind] = _sub( fout[fout_ind], t )
        fout[fout_ind] = _addto( fout[fout_ind], t )
        
        fout2_ind = fout2_ind + 1
        fout_ind = fout_ind + 1
        
        m -= 1
        if not(m):
            dowhile_if = 0
        
    return fout
 
 
def fft_bfly3(fout, fout_ind, fstride, twiddles, m):
    k = m
    m2 = 2*m
    scratch = [(0, 0), (0, 0), (0, 0), (0, 0)]
    epi3_i = twiddles[fstride*m][1]

    tw1_ind = 0
    tw2_ind = tw1_ind

    dowhile_if = 1
    while (dowhile_if):
        scratch[1] = _mult( fout[fout_ind+m],  twiddles[tw1_ind] )
        scratch[2] = _mult( fout[fout_ind+m2], twiddles[tw2_ind] )

        scratch[3] = _add( scratch[1], scratch[2] )
        scratch[0] = _sub( scratch[1], scratch[2] )
        tw1_ind = tw1_ind + fstride
        tw2_ind = tw2_ind + fstride*2

        fout[fout_ind+m] = ( fout[fout_ind][0] - (scratch[3][0])/2, \
            fout[fout_ind][1] - (scratch[3][1])/2 )

        scratch[0] = _mult_by_scalar( scratch[0], epi3_i )

        fout[fout_ind] = _addto( fout[fout_ind], scratch[3] )

        fout[fout_ind+m2] = ( fout[fout_ind+m][0] + scratch[0][1], \
            fout[fout_ind+m][1] - scratch[0][0] )

        fout[fout_ind+m] = ( fout[fout_ind+m][0] - scratch[0][1], \
            fout[fout_ind+m][1] + scratch[0][0] )

        fout_ind = fout_ind + 1

        k -= 1
        if not(k):
            dowhile_if = 0

    return fout


def fft_bfly_generic(fout, fout_ind, fstride, twiddles, m, p, nfft_orig):

    n_orig = nfft_orig

    # initialize scratch
    scratch = []
    for q1 in xrange(p): #( q1=0 ; q1<p ; ++q1 )
        scratch.append(0)

    for u in xrange(m): #( u=0; u<m; ++u )
        k = u
        for q1 in xrange(p): #( q1=0 ; q1<p ; ++q1 )
            scratch[q1] = fout[fout_ind+k]
            k = k + m

        k = u
        for q1 in xrange(p):
            twidx = 0
            fout[fout_ind+k] = scratch[0]
            for q in xrange(1, p):
                twidx = twidx + fstride * k
                if (twidx >= n_orig):
                    twidx = twidx - nfft_orig
 
                t = _mult( scratch[q],  twiddles[twidx] )
                fout[fout_ind+k] = _addto( fout[fout_ind+k], t )
            k = k + m
    return fout




def fft_alloc(nfft, inverse):
        
    twiddles = []
    
    for i in xrange(nfft):
        phase = -2*math.pi*float(i) / float(nfft)
        
        if (inverse):
            phase = phase * float(-1)
        
        twiddles.append(fft_cexp(phase))
    
    factors = fft_factor(nfft)
    
    return twiddles, factors


def fft_cexp(phase):
    x = (math.cos(phase), math.sin(phase))
    return x


def fft_factor(n):
    facbuf = []

    p = 4
    floor_sqrt = math.floor( math.sqrt( float(n) ) )

    # factor out powers of 4, powers of 2, then any remaining primes
    dowhile_test = 1
    while (dowhile_test):
        while n % p:
            if p == 4:
                p = 2
            elif p == 2:
                p = 3
            else:
                p = p + 2

            if (p > floor_sqrt):
                p = n # no more factors, skip to end
            
        n = n / p
        facbuf.append((p, n))
        
        if not(n > 1):
            dowhile_test = 0

    return facbuf


def _mult( a, b ):
    return ( a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0] )
        
def _sub( a, b ):
    return ( a[0]-b[0], a[1]-b[1] )

def _add( a, b ):
    return ( a[0] + b[0], a[1] + b[1] )

def _addto( res , a):
    return ( res[0] + a[0], res[1] + a[1] )

def _mult_by_scalar( c, s ):
    return ( c[0] * s, c[1] * s)
   

def main():
    fin = [(0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (0, 0)]
    inverse = 0
    print fft(fin, inverse)

if __name__ == '__main__':
    main()
