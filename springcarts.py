#!/usr/bin/env python
#------------------------------------------------------------------------------
# MCEN40024 SOLID MECHANICS (2012)
#  DEPARTMENT OF MECHANICAL ENGINEERING
#  THE UNIVERSITY OF MELBOURNE
#     ___________ 
#    / __/ __/ _ |
#   / _// _// __ |
#  /_/ /___/_/ |_|  - Finite Element Analysis
#
# This code may be used freely under the Creative Commons Attribution Unported
# (CC BY 3.0) license:  http://creativecommons.org/licenses/by/3.0/
# Cite the attribution as: 
#  "Portions of this code written by Dr Jonathan Merritt"
#
# Original Author: Jonathan Merritt <j.s.merritt@gmail.com>
#                                   <merritt@unimelb.edu.au>
#------------------------------------------------------------------------------

"""
==== ANALYSIS OF A SPRING-CART SYSTEM ====
This code performs three different types of analysis on the spring cart system
developed during Lecture 1 of the FEA component of the course:
  1. A transient dynamic analysis.
  2. A static analysis.
  3. An undamped, unforced modal frequency analysis.
For the development of the system equations, please see the lecture slides.
The code demonstrates an application of the Direct Stiffness Method to a
particular example system.

== BOUNDARY CONDITIONS ==
= TRANSIENT DYNAMICS AND STATIC ANALYSIS =
For these analyses, the displacement of cart 1 is fixed at zero.  A force is
applied in the positive direction to cart 4.
= MODAL FREQUENCY ANALYSIS =
In this analysis, no boundary conditions or forces are applied.

== UNITS ==
This example uses the following units:
  mm - length
  s  - time
  N  - force
  T  - mass

== PYTHON NOTES ==
You need relatively recent versions of scipy and matplotlib installed to run 
this code.
"""


# Attempt to import scipy and matplotlib; report an error if they're not
#  present
try:
    import numpy
    import scipy
    import matplotlib
    from numpy import array, asarray, bmat, diag, dot, eye, pi, sort, zeros
    from numpy.linalg import eig, inv
    from scipy.integrate import ode
    from matplotlib.pyplot import figure, legend, plot, show, xlabel, xlim, \
        ylabel, ylim
    from matplotlib import rc
except ImportError, e:
    print "You need to have NumPy / SciPy and matplotlib installed."


# Try to import matplotlibtotikz for figure export.  If it doesn't exist,
#  don't worry, we'll skip using it later.  This is used to generate the
#  figures in a pgfplots format, suitable for import into the lecture notes.
HAVE_TIKZ_EXPORT = False
try:
    from matplotlib2tikz import save as tikz_save
    HAVE_TIKZ_EXPORT = True
except:
    pass


def systemMatrices(k1, k2, k3, k4, k5, c6, m1, m2, m3, m4):

    """
    Computes the system matrices.

    This function simply assembles the system matrices in the form provided in
    the lecture.

    Args:
        k1, k2, k3, k4, k5: spring stiffnesses
        c6: viscous damping coefficient
        m1, m2, m3, m4: cart masses

    Returns:
        (K, C, M) a tuple of numpy arrays containing the system matrices,
        where K is the global stiffness matrix, C is the global damping matrix,
        and M is the global mass matrix.
    """

    # global system stiffness matrix
    K = array([[ k1,         -k1,        0,     0],
               [-k1, k1+k2+k3+k4,   -k2-k3,   -k4],
               [  0,      -k2-k3, k2+k3+k5,   -k5],
               [  0,         -k4,      -k5, k4+k5]], dtype='float64')

    # global system damping matrix
    Cz = zeros((4, 4), dtype='float64')
    Cz[0:2, 0:2] = array([[c6, -c6], [-c6, c6]], dtype='float64')
    C = Cz

    # global system mass matrix
    M = diag([m1, m2, m3, m4])

    return (K, C, M)


def transientDynamics(K, C, M, R, tstep, tf, d0, fixedDOF):

    """
    Analysis of transient dynamics.

    This function performs forward integration of the system equations over a
    given time interval.  Fixed (time-invariant) boundary conditions are
    supported by enforcing zero velocity and acceleration on any non-moving
    degrees of freedom.

    Args:
        K, C, M: system stiffness, damping and mass matrices (NxN)
        R: applied force matrix (Nx1)
        tstep: time-step to use
        tf: final time of integration
        d0: initial conditions (Nx1)
        fixedDOF: column vector containing 1.0 in any degrees of freedom that
          are fixed (Nx1)
    Returns:
        (T, d): a tuple of time samples and their corresponding displacement
          vectors
    """

    # find the size of our system (NB: I am assuming C and M are the correct
    #  sizes; not doing any checks on them)
    nDOF, n = K.shape
    assert (nDOF == n), 'K must be a square matrix'

    # zero and identity matrices of appropriate sizes for assembly
    Z = zeros((nDOF, nDOF))
    I = eye(nDOF)

    # Mbar and its inverse
    Mbar = bmat([[I, Z], [Z, M]])
    MbarInv = inv(Mbar)

    # Kbar and Rbar
    Kbar = bmat([[Z, I], [-K, -C]])
    Rbar = bmat([[zeros((nDOF, 1))], [R]])

    # Initial conditions (displacements are specified, velocities are zero)
    q0 = asarray(bmat([zeros((1, nDOF)), d0.transpose()])).flatten()

    # Matrix to set the derivatives of any fixed DOF to zero.  Mobility is a
    #  matrix which multiplies the output derivatives of all "active" DOF
    #  by 1, and all "fixed" or "inactive" DOF by 0.  If the first and
    #  second derivatives of a DOF are zero, then it will remain fixed during
    #  the forward integration of the system equations, thus enforcing the
    #  fixed boundary condition.
    MobilityD = (asarray(fixedDOF != 0, dtype='float64') * -1.0) + 1.0
    MobilityD = MobilityD.transpose()
    Mobility = diag(asarray(bmat([MobilityD, MobilityD])).flatten())
    
    # First-order ODE function for integration
    #  Mobility * MbarInv * (Kbar * q + Rbar)
    def qdot(t, q):
        # q is passed in as a row vector, so we need to convert it to a column
        #  vector
        qcol = asarray([q]).transpose()
        # evaluate the matrix equation
        qdotmat = dot(Mobility, dot(MbarInv, dot(Kbar, qcol) + Rbar))
        # flatten qdotmat back out to an array before returning
        return qdotmat.flatten()

    # Solve the ODE
    # time increment: 30 frames-per-second, 16 subframes-per-frame
    t0 = 0.0
    integrator = ode(qdot).set_integrator('dopri5', 
        atol=1E-4, rtol=1E-4, first_step=0.1, max_step=0.1)
    integrator.set_initial_value(q0, t0)
    T = []
    d = []
    while integrator.successful() and integrator.t < tf:
        Tnext = integrator.t + tstep
        integrator.integrate(Tnext)
        T.append(Tnext)
        d.append(integrator.y[0:nDOF])  # grab only the d components of q

    return (asarray(T),asarray(d))


if __name__ == '__main__':

    K, C, M = systemMatrices(
        # spring constants, k1 ... k6 (N/mm)
        5, 0.1, 0.1, 3, 0.2,
        # viscous damping coefficient (N/mm/s)
        1,
        # cart masses, m1 ... m4 (T)
        5E-3, 21E-3, 12E-3, 22E-3)


    #--------------------#
    # Transient dynamics #
    #--------------------#
    # external forces along each DOF (N)
    R = array([[0, 0, 0, 100]], dtype='float64').transpose()
    # initial conditions for all DOF (zero initial displacement)
    d0 = array([[0, 0, 0, 0]], dtype='float64').transpose()
    # fix the first DOF
    fixedDOF = array([[1, 0, 0, 0]], dtype='float64').transpose()
    # time step to use
    tStep = 0.01
    # final time of integration (s)
    tFinal = 10.0
    # perform the system integration
    T, d = transientDynamics(K, C, M, R, tStep, tFinal, d0, fixedDOF)


    #-----------------#
    # Static analysis #
    #-----------------#
    # essential boundary conditions / known displacements
    d_E = asarray([[0]])
    # force vector for the free nodes; f_F = [ R2, R3, R4 ] = [ 0, 0, R4 ]
    f_F = R[1:]
    # partition K into [ KE KEF ; KEF' KF].  we're taking a special case here,
    #  since we know that we have one essential boundary condition
    #  corresponding to the first DOF.  this would need to be more general
    #  in the case of other BCs.
    KE  = asarray([[K[1,1]]])
    KEF = asarray([K[0,1:]])  # row 0, column 1 to end
    KF  = K[1:,1:]            # rows 1->end, 2->end
    # solve for unknown static displacements
    d_F = dot(inv(KF), (f_F - dot(KEF.transpose(), d_E)))
    # solve for the unknown reactions at node 1
    rE_1 = dot(KE, d_E) + dot(KEF, d_F)


    #-------------------#
    # Modal frequencies #
    #-------------------#
    # solve the eigenvalue problem
    lambd, d_Modal = eig(dot(inv(M), K))  # "lambda" is a Python reserved word
    lambd = lambd[abs(lambd) > 1E-4]  # pick only non-zero frequencies
    omega = sort(lambd ** (0.5))
    freq = omega / (2*pi)
    print("Un-damped, un-forced modal frequencies (Hz):")
    print(freq)


    #--------------------------------#
    # Plot dynamic vs static results #
    #--------------------------------#
    # set plotting options
    lw = 2
    color1 = '#5FEA7C'
    color2 = '#F9E65C'
    color3 = '#DA61D2'
    color4 = '#4CABFF'
    fig = figure()
    xlabel('Time (s)')
    ylabel('Displacement (mm)')
    # plot the dynamic response curves
    plot(T, d[:,0], color=color1, lw=lw)
    plot(T, d[:,1], color=color2, lw=lw)
    plot(T, d[:,2], color=color3, lw=lw)
    plot(T, d[:,3], color=color4, lw=lw)
    # plot the static response horizontal lines
    plot([0, tFinal], [d_F[0], d_F[0]], color=color2, lw=lw, ls='--')
    plot([0, tFinal], [d_F[1], d_F[1]], color=color3, lw=lw, ls='--')
    plot([0, tFinal], [d_F[2], d_F[2]], color=color4, lw=lw, ls='--')
    xlim(0, tFinal)
    ylim(-10, 80)
    legend([r'$u_1$', r'$u_2$', r'$u_3$', r'$u_4$'])
    if HAVE_TIKZ_EXPORT:
        tikz_save('springcarts-plot.tex', 
            figureheight='7cm', figurewidth='10cm')
    show()
