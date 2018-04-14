import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import BSpline
from scipy import interpolate

#---------------------------CORRECTOR---------------------------#


def rescale_centroids(col,row):
    col = col - np.mean(col)
    row = row - np.mean(row)
    return col, row

def get_eigen_vectors(centroid_col, centroid_row):
    '''get the eigenvalues and eigenvectors given centroid x, y positions'''
    centroids = np.array([centroid_col, centroid_row])
    eig_val, eig_vec = np.linalg.eigh(np.cov(centroids))
    return eig_val, eig_vec

def rotate(eig_vec, centroid_col, centroid_row):
    """Rotate the coordinate frame of the (col, row) centroids to a new (x,y)
    frame in which the dominant motion of the spacecraft is aligned with
    the x axis.  This makes it easier to fit a characteristic polynomial
    that describes the motion."""
    centroids = np.array([centroid_col, centroid_row])
    return np.dot(eig_vec, centroids)

def fit_curve(rot_rowp, rot_colp, order):
    z = np.polyfit(rot_rowp, rot_colp, order)
    p = np.poly1d(z)
    p_deriv = p.deriv()
    return p,p_deriv

def arclength(x_prime, x_dense, p_deriv):
    """
    Compute the arclength of the polynomial used to fit the centroid
    measurements.

    Parameters
    ----------
    x_prime : float
        Upper limit of the integration domain.
    x_dense : ndarray
        Domain at which the arclength integrand is defined.

    Returns
    -------
    arclength : float
        Result of the arclength integral from x[0] to x1.
    """

    s = []
    for i in x_prime:
        gi = x_dense < i
        s_integrand = np.sqrt(1 + p_deriv(x_dense[gi]) ** 2)
        s.append(np.trapz(s_integrand, x=x_dense[gi]))
    return np.array(s)

def find_thruster_events(time,data,Xc,Yc):
    '''
    Find events when the spacecruft thruster are fired.
    Usually no useful data points are gathered when this happens
    '''

    diff_centroid = np.diff(Xc)**2 + np.diff(Yc)**2

    thruster_mask = diff_centroid < (1.5*np.mean(diff_centroid) + 0.*np.std(diff_centroid))

    # this little trick helps us remove 2 data points each time instead of just 1
    thruster_mask1 = np.insert(thruster_mask,0, False)
    thruster_mask2 = np.append(thruster_mask,False)
    thruster_mask = thruster_mask1*thruster_mask2

    time_thruster = time[ thruster_mask ]
    diff_centroid_thruster = diff_centroid[ thruster_mask[1:] ]

#     Xc_clipped = Xc[:][thruster_mask]
#     Yc_clipped = Yc[:][thruster_mask]
#     time_clipped = time[:][thruster_mask]
#     data_clipped = data[:][thruster_mask]

    return ~thruster_mask



def apply_sff(times, raw_fluxes, centroid_col, centroid_row,
        polyorder=5, bins=15, showfig=False, return_trend=False):
    """Returns a systematics-corrected LightCurve.

    Note that it is assumed that time and flux do not contain NaNs.

    Parameters
    ----------
    times : array-like
        Time measurements
    raw_fluxes : array-like
        Data flux for every time point
    centroid_col, centroid_row : array-like, array-like
        Centroid column and row coordinates as a function of time
    polyorder : int
        Degree of the polynomial which will be used to fit one
        centroid as a function of the other.
    bins : int
        Number of bins to be used in step (6) to create the
        piece-wise interpolation of arclength vs flux correction.
    windows : int
        Number of windows to subdivide the data.  The SFF algorithm
        is ran independently in each window.
    sigma_1, sigma_2 : float, float
        Sigma values which will be used to reject outliers
        in steps (6) and (2), respectivelly.

    Returns
    -------
    corrected_lightcurve : LightCurve object
        Returns a corrected lightcurve object.
    """

    col,row = rescale_centroids(centroid_col,centroid_row)
    eig_val, eig_vec = get_eigen_vectors(col, row)

    platescale = 4.0 # The Kepler plate scale; has units of arcseconds / pixel
    if showfig:
        v1, v2 = eig_vec

        pl.figure()#figsize=(5, 6))
        pl.plot(col * platescale, row * platescale, 'ko', ms=4)
        pl.plot(col * platescale, row * platescale, 'ro', ms=1)
        pl.xticks([-2, -1,0, 1, 2])
        pl.yticks([-2, -1,0, 1, 2])
        pl.xlabel('X position [arcseconds]')
        pl.ylabel('Y position [arcseconds]')
        #pl.xlim(-2, 2)
        #pl.ylim(-2, 2)
        pl.plot([0, v1[0]], [0, v1[1]], color='blue', lw=3)
        pl.plot([0, v2[0]], [0, v2[1]], color='blue', lw=3);

    #(1) Rotate the centroid measurements onto the subspace spanned by the
    #eigenvectors of the centroid covariance matrix
    rot_colp, rot_rowp = rotate(eig_vec, col, row) #units in pixels

    if showfig:
        pl.figure()#figsize=(5, 6))
        pl.plot(rot_rowp * platescale, rot_colp * platescale, 'ko', ms=4)
        pl.plot(rot_rowp * platescale, rot_colp * platescale, 'ro', ms=1)
        pl.xticks([-2, -1,0, 1, 2])
        pl.yticks([-2, -1,0, 1, 2])
        pl.xlabel("X' position [arcseconds]")
        pl.ylabel("Y' position [arcseconds]")
        #pl.xlim(-2, 2)
        #pl.ylim(-2, 2)
        pl.plot([0, 1], [0, 0], color='blue')
        pl.plot([0, 0], [0, 1], color='blue');

    #(2) Fit a polynomial to the rotated centroids
    p, p_deriv = fit_curve(rot_rowp, rot_colp, order=polyorder)
    x_dense = np.linspace(np.min(rot_rowp),  np.max(rot_rowp), 2000)

    if showfig:
        pl.figure()
        pl.plot(rot_rowp, rot_colp, '.')
        pl.plot(x_dense, p(x_dense))
        pl.ylabel('Position along minor axis (pixels)')
        pl.xlabel('Position along major axis (pixels)')
        pl.title('Performance of polynomial regression')
#         pl.ylim(-0.1, 0.1);

    #(3) Compute the arclength of such polynomial
    N = 1.5 #thruster firing?
    interior_knots = np.arange(times[0]+N, times[0]+6, N)
    #(4) Fit a BSpline of the raw flux as a function of time

    t,c,k = interpolate.splrep(times, raw_fluxes, s=0, task=-1, t=interior_knots)
    bspl = BSpline(t,c,k)

    if showfig:
        pl.figure()
        pl.plot(times, raw_fluxes, '.', label='raw flux')
        pl.plot(times, bspl(times), label='trend')
        pl.xlabel('$t$ (days)')
        pl.ylabel('Raw Flux');

    #Plot the normalized flux versus arclength to see the position-dependent flux
    #(5) Normalize the raw flux by the fitted BSpline computed in step (4)
    fluxes = raw_fluxes/bspl(times)

    #Mask the data by keeping only the good samples
    mask = find_thruster_events(times,fluxes,rot_colp,rot_rowp);

    clean_fluxes = fluxes[~mask]

    al = arclength(rot_rowp[~mask], x_dense, p_deriv) * platescale

    sorted_inds = np.argsort(al)

    #(6) Bin and interpolate the normalized flux as a function of the arclength
    #interpolate flux versus arclength position in 15 bins of means
    #piecewise linear fit

    bins = 15
    knots = np.array([np.min(al)]+
                 [np.median(splt) for splt in np.array_split(al[sorted_inds], bins)]+
                 [np.max(al)])


    bin_means = np.array([clean_fluxes[sorted_inds][0]]+
                 [np.mean(splt) for splt in np.array_split(clean_fluxes[sorted_inds], 15)]+
                 [clean_fluxes[sorted_inds][-1]])

    zz = np.polyfit(al, clean_fluxes, polyorder)
    linfit = np.poly1d(zz)
    #al_dense = np.linspace(0, 4, 1000)
    interp_func = interpolate.interp1d(knots, bin_means)

    if showfig:
        pl.figure()#figsize=(5, 6))
        pl.plot(arclength(rot_rowp,x_dense,p_deriv)*platescale, fluxes, 'ko', ms=4)
        pl.plot(arclength(rot_rowp,x_dense,p_deriv)*platescale, fluxes, 'o', color='#3498db', ms=3)
        pl.plot(arclength(rot_rowp[mask],x_dense,p_deriv)*platescale, fluxes[mask], 'o', color='r', ms=3)
        pl.plot(np.sort(al), interp_func(np.sort(al)), '-', color='#e67e22')
        #pl.xticks([0, 1, 2, 3, 4])

        pl.minorticks_on()
        pl.xlabel('Arclength [arcseconds]')
        pl.ylabel('Relative Brightness')
        #pl.title('EPIC 60021426, Kp =10.3')
#         pl.xlim(0,4)
#         pl.ylim(0.997, 1.002);

    #(7) Divide the raw flux by the piecewise linear interpolation done in step (6)
    corr_fluxes = clean_fluxes / interp_func(al)

    #(8) Set raw flux as the flux computed in step (7) and repeat

    if showfig:
        pl.figure(figsize=(10,6))
        dy = 0.01
        pl.plot(times, raw_fluxes+dy, 'ko', ms=4)
        pl.plot(times, raw_fluxes+dy, 'o', color='#3498db', ms=3)
        pl.plot(times[mask], raw_fluxes[mask]+dy, 'o', color='r', ms=3)
        pl.plot(times[~mask], corr_fluxes*bspl(times[~mask]), 'o', color='k', ms = 4)
        pl.plot(times[~mask], corr_fluxes*bspl(times[~mask]), 'o', color='#e67e22', ms = 3)

        pl.xlabel('BJD - 2454833')
        pl.ylabel('Relative Brightness')

#         pl.xlim(1862, 1870)
#         pl.ylim(0.994, 1.008);

    if return_trend:
        return corr_fluxes, bspl(times[~mask])
    #(9) Multiply back the fitted BSpline
    else:
        return corr_fluxes*bspl(times[~mask])
