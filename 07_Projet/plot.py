#!/usr/bin/env python3

##########
import matplotlib.pyplot as plt
import meshio
import matplotlib.tri as tri
import subprocess
import numpy as np
import os
from IPython.display import IFrame, display
from Slides import math_helper as mh
from sympy import *
import matplotlib

################################################################
matplotlib.rcParams['figure.figsize'] = 9, 6
################################################################


def plot_matrix(matrix, matrix_name):
    mh.print_latex(matrix_name+"= {0}", Matrix(matrix))


def readMesh(filename, element_type='triangle'):
    mesh = meshio.read(filename)
    return mesh.points[:, :2], mesh.cells_dict[element_type]

def generate_scalar_field(field, axis=None):
    if axis is None:
        field = np.linalg.norm(field, axis=1)
    else:
        field = field[:, axis]
    return field

def plotMesh(coords, connectivity, nodal_field=None, elemental_field=None, **kwargs):
    triangles = tri.Triangulation(coords[:, 0], coords[:, 1], connectivity)
    plt.axes().set_aspect('equal')
    if nodal_field is not None:
        nodal_field = nodal_field.reshape(coords.shape[0], nodal_field.size//coords.shape[0])
        nodal_field = generate_scalar_field(nodal_field, **kwargs)
        contour = plt.tricontourf(triangles, nodal_field)
        plt.colorbar(contour)
    if elemental_field is not None:
        elemental_field = elemental_field.reshape(connectivity.shape[0], elemental_field.size//connectivity.shape[0])
        elemental_field = generate_scalar_field(elemental_field, **kwargs)
        contour = plt.tripcolor(triangles, elemental_field)
        plt.colorbar(contour)

    t = plt.triplot(triangles, '--', lw=.8)


def meshGeo(filename, dim=2, order=1):
    ret = subprocess.run(f"gmsh -2 -order 1 -o tmp.msh {filename}", shell=True)
    if ret.returncode:
        print("Beware, gmsh could not run: mesh is not generated")
    else:
        print("Mesh generated")
        mesh = readMesh('tmp.msh')
        os.remove('tmp.msh')
        return mesh
    return None


def votre_opinion_compte(name):
    display(IFrame('https://www.surveymonkey.com/r/NOTOSURVEY?notebook_set=CIVIL-321&notebook_id=CIVIL-321'+name, 600, 1000))
