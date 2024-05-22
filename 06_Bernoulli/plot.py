#!/usr/bin/env python3

##########

import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import numpy as np
from matplotlib import collections as mc
from sympy import Matrix, simplify
from Slides import math_helper as mh
from IPython.display import Image
from sympy import *
from IPython.display import IFrame, display


def plot_matrix(matrix, matrix_name):
    mh.print_latex(matrix_name+"= {0}", Matrix(matrix))


def plot_matrix_product(matrix1, matrix2, matrix3, matrix_name):
    mh.print_latex(matrix_name + '= {0}{1}{2}', matrix1, matrix2, matrix3)


def compute_minmax(positions):
    return np.array(positions.min(axis=0), dtype=float), \
        np.array(positions.max(axis=0), dtype=float)


def compute_range(positions):
    _min, _max = compute_minmax(positions)
    return _max - _min


def compute_node_size(positions):
    _range = compute_range(positions)
    _range = _range.max()
    return _range*.08


def packed_eqn(node_idx, dof_per_node, number_of_nodes):
    return np.arange(dof_per_node)+np.ones(dof_per_node)*node_idx*dof_per_node


def eqn_number_node(f_node_eqn, number_of_nodes):
    eqn_number = []
    for n in range(0, number_of_nodes):
        eqn = f_node_eqn(n, 2, number_of_nodes)
        eqn_number.append(np.array(eqn).flatten())
    return np.array(eqn_number, dtype=int)


def eqn_number_elem(eqn_number_node, conn):
    eqn_number = []
    for e in conn:
        eqn = []
        eqn.append(eqn_number_node[e[0], :])
        eqn.append(eqn_number_node[e[1], :])
        eqn_number.append(np.array(eqn).flatten())

    return np.array(eqn_number, dtype=int)


def convert_colors(cols, default_color='#1f77b4'):
    if cols is None:
        return default_color

    if isinstance(cols, str):
        if cols == '':
            return default_color
        return cols

    for i in range(len(cols)):
        if cols[i] == '':
            cols[i] = default_color

    cols = [pltcol.to_rgb(c) for c in cols]
    return cols


def plot_nodes(ax, positions, eqn_num_node=None, node_colors=None,
               node_size_px=None, **kwargs):

    _min, _max = compute_minmax(positions)
    _range = compute_range(positions)
    ax.set_xlim((_min[0]-_range.max()*0.2, _max[0]+_range.max()*0.2))
    ax.set_ylim((_min[1]-_range.max()*0.2, _max[1]+_range.max()*0.2))

    if node_size_px is None:
        disp_positions = ax.transData.transform(positions)
        node_size_px = compute_node_size(disp_positions)

    node_colors = convert_colors(node_colors)
    ax.scatter(positions[:, 0], positions[:, 1], s=node_size_px**2,
               c=node_colors, edgecolors='k')

    for i in range(2):
        if _range[i] == 0:
            _range[i] = .5

    # center_gravity = np.average(positions, axis=0)
    center_gravity = (_max+_min)/2
    # print(center_gravity)

    for i, p in enumerate(positions):
        _n = p - center_gravity
        norm = np.linalg.norm(_n)
        if norm < 1e-5:
            _n = np.array([1, 1])
        # print(norm, p, center_gravity,  _n)
        norm = _range.max()*.08/np.linalg.norm(_n)
        _n *= norm
        # print(p, _n)
        pos = p  # + _n
        # print('aaa', _range, p, pos, np.linalg.norm(_n))
        ax.text(pos[0], pos[1], str(i), horizontalalignment='center',
                verticalalignment='center')

        if eqn_num_node is not None:
            eqns = eqn_num_node[i, :]
            ax.text(pos[0]+_n[0]*1.5, pos[1]+_n[1]*1.5,
                    "[" + ",".join([str(int(e)) for e in eqns]) + "]",
                    horizontalalignment='center',
                    verticalalignment='center')

    return node_size_px


def N1(L, xi):
    return 1/L**3*(2*xi**3-3*xi**2*L+L**3)


def N2(L, xi):
    return 1/L**3*(xi**3*L-2*xi**2*L**2+xi*L**3)


def N3(L, xi):
    return 1/L**3*(-2*xi**3+3*xi**2*L)


def N4(L, xi):
    return 1/L**3*(xi**3*L-xi**2*L**2)


def create_element_lines(p1, p2, u1=None, u2=None, n=20, interpolate=None):
    L = np.linalg.norm(p2-p1)
    e1 = np.zeros(3)
    e3 = np.zeros(3)
    e3[2] = 1.
    e1[:2] = (p2-p1)/L
    e2 = np.cross(e3, e1)
    xi = np.arange(n+1)/n*L
    X = np.einsum('i,j->ij', N1(L, xi), p1)
    X += np.einsum('i,j->ij', N3(L, xi), p2)
    if (u1 is not None) and (u2 is not None):
        X += np.einsum('i,j->ij', N1(L, xi), u1[:2])
        X += np.einsum('i,j->ij', N3(L, xi), u2[:2])
        X += np.einsum('i,j->ij', N2(L, xi), e2[:2]*u1[2])
        X += np.einsum('i,j->ij', N4(L, xi), e2[:2]*u2[2])
    res = [x for x in zip(X[:-1], X[1:])]
    return res


def plot_elements(ax, positions, conn, displacement=None,
                  linestyle=None, elem_colors=None, n=20,
                  show_number=True, **kwargs):
    lines = []

    for i, e in enumerate(conn):
        p1 = positions[e[0]]
        p2 = positions[e[1]]
        u1 = None
        u2 = None
        if displacement is not None:
            u1 = displacement[e[0]]
            u2 = displacement[e[1]]
        lines += create_element_lines(p1, p2, n=n, u1=u1, u2=u2)

    elem_colors = convert_colors(elem_colors)

    # if elem_colors is not None:
    #    elem_colors = [pltcol.to_rgb(c) for c in elem_colors]

    if linestyle is None:
        linestyle = '-'

    lc = mc.LineCollection(lines, linewidths=2,
                           linestyles=linestyle, colors=elem_colors)
    ax.add_collection(lc)

    node_size = compute_node_size(positions)
    # print(node_size)

    if not show_number:
        return

    for i, e in enumerate(conn):
        p1 = positions[e[0]]
        p2 = positions[e[1]]
        center = (p1 + p2)/2
        _l = np.zeros(3)
        _l[:2] = p2 - p1
        _l /= np.linalg.norm(_l)
        _n = np.cross(_l, [0, 0, 1])
        # print(_n, node_size)
        center += _n[:2]*node_size*1.5
        ax.text(center[0], center[1], f"({str(i)})",
                horizontalalignment='center',
                verticalalignment='center')


def plot_structure(positions, conn,
                   plot_eqn=None,
                   displacement=None,
                   **kwargs):

    positions = np.array(positions)
    conn = np.array(conn)
    eqn_num_node = None

    if type(plot_eqn) == str:
        if plot_eqn == "packed":
            plot_eqn = packed_eqn
        else:
            raise RuntimeError("not known eqn packing strategy")
        eqn_num_node = eqn_number_node(plot_eqn, positions.shape[0])

    elif type(plot_eqn) == np.ndarray:
        eqn_num_node = plot_eqn
    else:
        if plot_eqn is not None:
            raise RuntimeError("could not get eqn " + str(type(plot_eqn)))

    if plot_eqn is not None:
        eqn_num_elem = eqn_number_elem(eqn_num_node, conn)

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    if displacement is not None:
        displacement = np.array(displacement, dtype=float)
        plot_elements(ax, positions, conn, show_number=False,
                      elem_colors=pltcol.CSS4_COLORS['silver'],
                      linestyle='--', n=2, **kwargs)

        plot_elements(ax, positions, conn, displacement=displacement, **kwargs)

        node_size_px = plot_nodes(ax, positions,
                                  node_colors=pltcol.CSS4_COLORS['silver'],
                                  **kwargs)

        plot_nodes(ax, positions + displacement[:, :2],
                   node_size_px=node_size_px, **kwargs)
    else:
        plot_nodes(ax, positions, **kwargs)
        plot_elements(ax, positions, conn, **kwargs)

    if plot_eqn is None:
        return

    ret = {}
    ret['eqn_node'] = np.array(eqn_num_node, dtype=int)
    ret['eqn_elem'] = np.array(eqn_num_elem, dtype=int)
    return ret


################################################################

def calculerMatriceRotation(p1, p2):
    # Calcul de la matrice de rotation du repère local vers le repère global
    # La fonction retourne également L pour le calcul de la rigidité de
    # l'élément e global coordonnees connectivites

    barre = Matrix(p2 - p1)
    L = barre.norm()

    R = 1/L * np.array([[barre[0], -barre[1]],
                        [barre[1], barre[0]]])

    return simplify(Matrix(R)), L

################################################################


def calculerMatriceRigiditeLocale(k):

    Kl = k*Matrix([[1, 0, -1, 0],
                   [0, 0, 0, 0],
                   [- 1, 0, 1, 0],
                   [0, 0, 0, 0]])
    return simplify(Kl)

################################################################


def assemblerMatriceRigidite(coordonnees, connectivites, materiau,
                             elem_colors=None, elem_to_assemble=None,
                             **kwargs):

    eqn_num_node = eqn_number_node(packed_eqn, coordonnees.shape[0])
    equations_num = eqn_number_elem(eqn_num_node, connectivites)

    nb_noeuds = coordonnees.shape[0]
    nb_elem = connectivites.shape[0]
    nb_noeuds_p_elem = connectivites.shape[1]
    nb_ddl_p_noeuds = 2

    K = np.zeros((nb_noeuds*nb_ddl_p_noeuds, nb_noeuds *
                  nb_ddl_p_noeuds), dtype=object)

    if elem_to_assemble is None:
        elem_to_assemble = range(0, nb_elem)

    for e in elem_to_assemble:
        n_1 = coordonnees[connectivites[e, 0], :]
        n_2 = coordonnees[connectivites[e, 1], :]
        R, L = calculerMatriceRotation(n_1, n_2, **kwargs)
        T = np.zeros((nb_noeuds_p_elem*nb_ddl_p_noeuds,
                      nb_noeuds_p_elem*nb_ddl_p_noeuds), dtype=object)
        T[: 2, :2] = R
        T[2:, 2:] = R
        # Construction de la matrice de rigidité locale
        Kl = calculerMatriceRigiditeLocale(materiau[e]/L)
        # Rotation de la matrice dans le système global de coordonnée
        Kg = T@Kl@T.T
        # Assemblage dans la matrice de rigidité du système à l'aide des
        idx = equations_num[e, :]
        for i, gi in enumerate(idx):
            for j, gj in enumerate(idx):
                K[gi, gj] += Kg[i, j]
    K = simplify(Matrix(K))
    if elem_colors is None:
        return K

    # elem_colors = np.array(convert_colors(elem_colors))
    # print(elem_colors)

    Kcolor = np.zeros_like(K, dtype='S10')
    for e in elem_to_assemble:
        idx = equations_num[e, :]
        for i, gi in enumerate(idx):
            for j, gj in enumerate(idx):
                c1 = Kcolor[gi, gj].decode()
                c2 = elem_colors[e]

                if c1 != c2:
                    c3 = c1+c2
                else:
                    c3 = c1
                Kcolor[gi, gj] = c3.encode()

    K = mh.ColoredMatrix(K)
    K.colors = Kcolor
    return K

################################################################

def votre_opinion_compte(name):
    display(IFrame('https://www.surveymonkey.com/r/NOTOSURVEY?notebook_set=CIVIL-321&notebook_id=CIVIL-321'+name, 600, 1000))

def plot_bloc_stiffness(positions, conn, E, A, nb_barre):

    materiau = E*A*np.ones(conn.shape[0], dtype=object)
    elem_colors = np.array(['b', 'g', 'r', 'y'])

    for i in range(nb_barre):
        colors = np.ones(conn.shape[0], dtype=str)
        colors[:] = 'k'
        colors[i] = elem_colors[i]
        ret = plot_structure(positions, conn, elem_colors=colors)

        K = assemblerMatriceRigidite(
            positions, conn, materiau, elem_colors=elem_colors, elem_to_assemble=range(i+1))
        K.evalf(5)
        display(K.profile(remove_zeros=True))
        # display(K.profile())
