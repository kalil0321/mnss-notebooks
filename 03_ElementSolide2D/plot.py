#!/usr/bin/env python3

from IPython.display import IFrame, display
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy as np
import meshio
from sympy import *
import os
import subprocess
##########


def readMesh(filename):
    mesh = meshio.read(filename)
    return mesh.points[:, :2], np.array([i for i in mesh.cells[0].data])


def plotMesh(coords, connectivity, field=None):
    triangles = tri.Triangulation(coords[:, 0], coords[:, 1], connectivity)
    plt.axes().set_aspect('equal')
    if field is not None:
        plt.tricontourf(triangles, np.linalg.norm(
            field.reshape(field.shape[0]//2, 2), axis=1))
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

##########
