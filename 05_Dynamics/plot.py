#!/usr/bin/env python3

##########

from IPython.display import IFrame, display
from Slides import math_helper as mh
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib
from Slides import presentation_helper as ph
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
from sympy import *

################################################################
matplotlib.rcParams['figure.figsize'] = 9, 6
matplotlib.rcParams.update({'font.size': 22})
matplotlib.rcParams.update({'legend.fontsize': 22})
matplotlib.rcParams.update({'lines.linewidth': 4})
matplotlib.rcParams.update({'lines.markersize': 10})
matplotlib.rcParams.update({'axes.labelsize': 30})
################################################################


def plot_matrix(matrix, matrix_name):
    mh.print_latex(matrix_name+"= {0}", Matrix(matrix))

################################################################


def spring_animation(nodes, disp, xlim=None, ylim=(-0.5, 0.5)):
    nsteps = disp.shape[1]

    if xlim is None:
        xlim = (nodes.min(), nodes.max())
    # La figure, les axes et les éléments à tracer sont définis.
    fig = plt.figure()
    ax = plt.axes(xlim=xlim, ylim=ylim)
    line, = ax.plot([], [], '-')

    # Quelques fonctions pour la visualisation.
    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        line.set_data(nodes, disp[i, :])
        return line,

    # Animator (blit=True means only re-draw the parts that have changed)
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=int(nsteps), interval=20, blit=True)
    # Visualisation
    return ph.display_animation(anim)

################################################################


def computeEnergyEvolution(displacements, velocities, computeEnergy):
    epot = []
    ekin = []
    etot = []

    nsteps = displacements.shape[0]
    for i in range(nsteps):
        _epot, _ekin = computeEnergy(displacements[i], velocities[i])
        _etot = _epot+_ekin

        epot.append(_epot)
        ekin.append(_ekin)
        etot.append(_etot)

    epot = np.array(epot)
    ekin = np.array(ekin)
    etot = np.array(etot)

    return(epot, ekin, etot)

################################################################


def plotEnergyEvolution(displacements, velocities, computeEnergy, loglog=False,
                        semilog=False, **kwargs):
    nsteps = displacements.shape[0]

    epot, ekin, etot = computeEnergyEvolution(
        displacements, velocities, computeEnergy)

    # Représentation graphique des énergies
    T = range(nsteps)
    if semilog:
        plt.semilogy(T, etot, label='$E^{Tot}$', **kwargs)
    elif loglog:
        plt.loglog(T, etot, label='$E^{Tot}$', **kwargs)
    else:
        plt.plot(T, epot, label='$E^{pot}$', **kwargs)
        plt.plot(T, ekin, label='$E^{cin}$', **kwargs)
        plt.plot(T, etot, label='$E^{Tot}$', **kwargs)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', **kwargs)
    plt.ylabel('Energy [J]', **kwargs)
    plt.xlabel('Time frame [itérations]', **kwargs)

################################################################


def modesAnimation(nodes, xlim=None, ylim=(-0.5, 0.5), eigs=None,
                   marker=None):

    eigen_values = eigs[0].real
    eigen_vectors = eigs[1]

    if nodes.shape[0] != eigen_vectors[:, 0].shape[0]:
        nodes = nodes[1:-1]

    eigen_values[eigen_values < 0] = 0
    omegas = np.sqrt(eigen_values)
    T = 2.*np.pi/omegas[0]
    if np.abs(omegas[0]) < 1e-12:
        T = 2.*np.pi/omegas[1]

    n_modes = eigen_values.shape[0]

    nsteps = 200*n_modes
    time_factor = 2*T/nsteps

    if xlim is None:
        xlim = (nodes.min(), nodes.max())

    fig = plt.figure()
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.15)
    ax = plt.axes(xlim=xlim, ylim=ylim)
    ax.set_xlabel('Position $x$')
    ax.set_ylabel('Deplacement $d(t)$')
    lines = [ax.plot(nodes, eigen_vectors[:, 0].real, marker=marker, lw=2)[0]
             for i in range(0, n_modes)]

    def init():
        for i in range(0, n_modes):
            lines[i].set_data([], [])
        return lines

    def animate(ts):
        nn_modes = 1 + int(ts/200)
        if nn_modes > n_modes:
            nn_modes = n_modes

        for i in range(0, nn_modes):
            lines[i].set_data(nodes, eigen_vectors[:, i].real *
                              np.cos(omegas[i]*ts*time_factor))

        return lines

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=nsteps, interval=10, blit=True)

    return anim

################################################################


def NewMarkIntegrationBeta025Gamma05(U, V, A, M, K, dt):
    Fprime = M@(U + dt*V+dt**2/4.*A)
    Mprime = M + dt**2/4*K
    Mprime = scipy.sparse.csr_matrix(Mprime)
    U[:] = scipy.sparse.linalg.spsolve(Mprime, Fprime)
    V += dt/2*A
    Fprime = -K@U
    Mprime = scipy.sparse.csr_matrix(M)
    new_A = scipy.sparse.linalg.spsolve(Mprime, Fprime)
    A[:] = new_A[:]
    V += dt/2*A
################################################################


def makeEvolution(nodes, nsteps, dt=1.,
                  U=None, K=None, M=None,
                  time_integration=NewMarkIntegrationBeta025Gamma05):

    displacements = []
    forces = []
    velocities = []

    V = np.zeros_like(nodes)
    A = np.zeros_like(nodes)
    if callable(U):
        U = U(nodes)
    if callable(K):
        K = K(nodes.shape[0])
    if callable(M):
        M = M(nodes.shape[0])

    K = scipy.sparse.csr_matrix(K)
    M = scipy.sparse.csr_matrix(M)

    for s in range(0, nsteps):
        displacements.append(U.copy())
        velocities.append(V.copy())
        forces.append((K@U).copy())
        if s % 100 == 0:
            print(s)
        time_integration(U, V, A, M, K, dt)

    displacements = np.array(displacements)
    velocities = np.array(velocities)
    forces = np.array(forces)
    return displacements, velocities, forces
################################################################


def votre_opinion_compte(name):
    display(IFrame('https://www.surveymonkey.com/r/NOTOSURVEY?notebook_set=CIVIL-321&notebook_id=CIVIL-321'+name, 600, 1000))
