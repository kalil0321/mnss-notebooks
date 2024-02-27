#!/usr/bin/env python3

##########

from IPython.display import IFrame, display
from sympy import Matrix
from Slides import math_helper as mh


def plot_matrix(matrix, matrix_name):
    mh.print_latex(matrix_name+"= {0}", Matrix(matrix))


def plot_matrix_product(matrix1, matrix2, matrix3, matrix_name):
    mh.print_latex(matrix_name + '= {0}{1}{2}', matrix1, matrix2, matrix3)
    
def votre_opinion_compte():
    display(IFrame('https://docs.google.com/forms/d/e/1FAIpQLSfYU5qWjevZSh2jLFoN05yXHcjpLAdAaTmLlpVrKX-1gfuf_A/viewform?usp=sf_link', 600, 1000))