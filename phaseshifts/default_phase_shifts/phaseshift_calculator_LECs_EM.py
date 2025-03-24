#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 16:05:07 2022

@author: tom & yannick
"""

# import of necessary modules
import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt, sqrt, zeros, array, exp
from numpy import linalg
import scipy as sc
from scipy.integrate import odeint
import random
import os





operator_path = '/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/operators'

potential_path = '/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/EM500/N3LO'


# physical constants
hbarc = 197.326
M = (938.272 + 939.565) / 2.0  # averaged neutron/proton mass in MeV
units_factor = hbarc * hbarc / M
Nrows = 100  # 100: accurate, SRG takes some time! 50: less accurate, but faster

# arrays for uncoupled interaction channels
V = zeros([Nrows, Nrows], float)
Vmat = zeros([Nrows, Nrows], float)

# arrays for coupled interaction channels
Vcoupled = zeros([2 * Nrows, 2 * Nrows], float)
Vmatcoupled = zeros([2 * Nrows, 2 * Nrows], float)

# matrices for kinetic energy
Tkin = zeros([Nrows, Nrows], float)
Tkincoupled = zeros([2 * Nrows, 2 * Nrows], float)


# read in nuclear interaction matrix elements for given uncoupled partial wave channels







def read_Vchiral(file, Nrows):

    mesh = np.genfromtxt(file, dtype=(float, float), skip_header=0, max_rows=Nrows)
    mesh_weights = mesh[:, 0]
    mesh_points = mesh[:, 1]



    Vread = np.genfromtxt(file, dtype=(float, float, float), skip_header=Nrows, max_rows=Nrows * Nrows)
    V = zeros([Nrows, Nrows], float)
    Vmat = zeros([Nrows, Nrows], float)
    Tkin = zeros([Nrows, Nrows], float)

    for i in range(Nrows):
        Tkin[i, i] = mesh_points[i] * mesh_points[i]
        for j in range(Nrows):
            V[i, j] = Vread[i * Nrows + j][2]
            Vmat[i, j] = 2.0 / np.pi * sqrt(mesh_weights[i]) * sqrt(mesh_weights[j]) * mesh_points[i] * mesh_points[j] * \
                         V[i, j]
    return [V, Vmat, Tkin, mesh_points, mesh_weights]


# routines for computing phase shifts by solving the Lippmann Schwinger equation

# offset for subtracting the pole
eps = 1e-3


# computation of phase shifts in uncoupled channels
def counterterm(pmax, pole):
    return np.arctanh(pole / pmax) / pole


def compute_phase_shifts(K, i, mesh_points):
    return 180.0 / np.pi * np.arctan(-mesh_points[i] * K[i, i])


# note that ambiguity regarding the brach of arctan exists, make sure that solution is continuous as a function of energy
def compute_phase_shifts_coupled(K, i, mesh_points):
    epsilon = np.arctan(2 * K[i, i + Nrows] / (K[i, i] - K[i + Nrows, i + Nrows])) / 2.0
    r_epsilon = (K[i, i] - K[i + Nrows, i + Nrows]) / (np.cos(2 * epsilon))
    delta_a = - np.arctan(mesh_points[i] * (K[i, i] + K[i + Nrows, i + Nrows] + r_epsilon) / 2.0)
    delta_b = - np.arctan(mesh_points[i] * (K[i, i] + K[i + Nrows, i + Nrows] - r_epsilon) / 2.0)

    epsilonbar = np.arcsin(np.sin(2 * epsilon) * np.sin(delta_a - delta_b)) / 2.0
    delta_1 = 180.0 / np.pi * (delta_a + delta_b + np.arcsin(np.tan(2 * epsilonbar) / (np.tan(2 * epsilon)))) / 2.0
    delta_2 = 180.0 / np.pi * (delta_a + delta_b - np.arcsin(np.tan(2 * epsilonbar) / (np.tan(2 * epsilon)))) / 2.0

    epsilonbar *= -180.0 / np.pi
    return [delta_1, delta_2, epsilonbar]


def delta(i, j):
    if (i == j):
        return 1
    else:
        return 0


def Elab(p):
    return 2 * p ** 2 * hbarc ** 2 / M


def mom(E):
    return np.sqrt(M * E / 2 / hbarc ** 2)


def compute_K_matrix(V, Nrows, pmax, mesh_points, mesh_weights):
    A = zeros([Nrows + 1, Nrows + 1], float)
    K = zeros([Nrows, Nrows], float)

    for x in range(Nrows):
        pole = mesh_points[x] + eps

        for i in range(Nrows):
            for j in range(Nrows):
                A[i, j] = delta(i, j) - 2.0 / np.pi * mesh_weights[j] * V[i, j] * mesh_points[j] ** 2 / (
                            pole ** 2 - mesh_points[j] ** 2)

        sum = 0.0
        for i in range(Nrows):
            sum += mesh_weights[i] / (pole ** 2 - mesh_points[i] ** 2)

        for i in range(Nrows):
            A[Nrows, i] = - 2.0 / np.pi * mesh_weights[i] * mesh_points[i] ** 2 * V[x, i] / (
                        pole ** 2 - mesh_points[i] ** 2)
            A[i, Nrows] = + 2.0 / np.pi * V[i, x] * pole ** 2 * (sum - counterterm(pmax, pole))

        A[Nrows, Nrows] = 1 + 2.0 / np.pi * V[x, x] * pole ** 2 * (sum - counterterm(pmax, pole))

        bvec = zeros([Nrows + 1], float)
        for i in range(Nrows):
            bvec[i] = V[i, x]
        bvec[Nrows] = V[x, x]

        xvec = np.linalg.solve(A, bvec)

        for i in range(Nrows):
            K[i, x] = xvec[i]

    return K


def compute_K_matrix_coupled(V00, V01, V10, V11, Nrows, pmax, mesh_points, mesh_weights):
    A = zeros([2 * Nrows + 2, 2 * Nrows + 2], float)
    K = zeros([2 * Nrows, 2 * Nrows], float)

    for x in range(Nrows):
        pole = mesh_points[x] + eps

        for i in range(Nrows):
            for j in range(Nrows):
                A[i, j] = delta(i, j) - 2.0 / np.pi * mesh_weights[j] * V00[i, j] * mesh_points[j] ** 2 / (
                            pole ** 2 - mesh_points[j] ** 2)
                A[i, j + Nrows] = - 2.0 / np.pi * mesh_weights[j] * V01[i, j] * mesh_points[j] ** 2 / (
                            pole ** 2 - mesh_points[j] ** 2)
                A[i + Nrows, j] = - 2.0 / np.pi * mesh_weights[j] * V10[i, j] * mesh_points[j] ** 2 / (
                            pole ** 2 - mesh_points[j] ** 2)
                A[i + Nrows, j + Nrows] = delta(i, j) - 2.0 / np.pi * mesh_weights[j] * V11[i, j] * mesh_points[
                    j] ** 2 / (pole ** 2 - mesh_points[j] ** 2)

        sum = 0.0
        for i in range(Nrows):
            sum += mesh_weights[i] / (pole ** 2 - mesh_points[i] ** 2)

        for i in range(Nrows):
            A[2 * Nrows, i] = - 2.0 / np.pi * mesh_weights[i] * mesh_points[i] ** 2 * V00[x, i] / (
                        pole ** 2 - mesh_points[i] ** 2)
            A[i, 2 * Nrows] = + 2.0 / np.pi * V00[i, x] * pole ** 2 * (sum - counterterm(pmax, pole))

            A[2 * Nrows, i + Nrows] = - 2.0 / np.pi * mesh_weights[i] * mesh_points[i] ** 2 * V01[x, i] / (
                        pole ** 2 - mesh_points[i] ** 2)
            A[i, 2 * Nrows + 1] = + 2.0 / np.pi * V01[i, x] * pole ** 2 * (sum - counterterm(pmax, pole))

            A[2 * Nrows + 1, i] = - 2.0 / np.pi * mesh_weights[i] * mesh_points[i] ** 2 * V10[x, i] / (
                        pole ** 2 - mesh_points[i] ** 2)
            A[i + Nrows, 2 * Nrows] = + 2.0 / np.pi * V10[i, x] * pole ** 2 * (sum - counterterm(pmax, pole))

            A[2 * Nrows + 1, i + Nrows] = - 2.0 / np.pi * mesh_weights[i] * mesh_points[i] ** 2 * V11[x, i] / (
                        pole ** 2 - mesh_points[i] ** 2)
            A[i + Nrows, 2 * Nrows + 1] = + 2.0 / np.pi * V11[i, x] * pole ** 2 * (sum - counterterm(pmax, pole))

        A[2 * Nrows, 2 * Nrows] = 1 + 2.0 / np.pi * V00[x, x] * pole ** 2 * (sum - counterterm(pmax, pole))
        A[2 * Nrows, 2 * Nrows + 1] = + 2.0 / np.pi * V01[x, x] * pole ** 2 * (sum - counterterm(pmax, pole))
        A[2 * Nrows + 1, 2 * Nrows] = + 2.0 / np.pi * V10[x, x] * pole ** 2 * (sum - counterterm(pmax, pole))
        A[2 * Nrows + 1, 2 * Nrows + 1] = 1 + 2.0 / np.pi * V11[x, x] * pole ** 2 * (sum - counterterm(pmax, pole))

        bvec = zeros([2 * Nrows + 2, 2], float)
        for i in range(Nrows):
            bvec[i, 0] = V00[i, x]
            bvec[i, 1] = V01[i, x]
            bvec[i + Nrows, 0] = V10[i, x]
            bvec[i + Nrows, 1] = V11[i, x]

        bvec[2 * Nrows, 0] = V00[x, x]
        bvec[2 * Nrows, 1] = V01[x, x]
        bvec[2 * Nrows + 1, 0] = V10[x, x]
        bvec[2 * Nrows + 1, 1] = V11[x, x]

        xvec = np.linalg.solve(A, bvec)

        for i in range(Nrows):
            K[i, x] = xvec[i, 0]
            K[i, x + Nrows] = xvec[i, 1]
            K[i + Nrows, x] = xvec[i + Nrows, 0]
            K[i + Nrows, x + Nrows] = xvec[i + Nrows, 1]

    return K


def genenerate_wavename(S, L, J):
    dictio = {1: "P", 0: "S", 2: "D"}
    i = str(2 * S + 1)
    j = dictio[L]
    k = str(J)
    return "delta" + i + j + k + ".txt"


# read in particular uncoupled or coupled matrix elements

title2 = "LO_EM500new"
S =      [0]
L =      [0]
Lprime = [0]            #first index respectively 10010 3s1
J =      [0]
T =      [1]
SVD_rank = 4 #5 operators and 5 singular values



grid_size = 100

partial_wave = "$^1$D$_2$"

Nrows = 100
colors = ['r', 'navy', 'orange', 'teal', 'forestgreen']











#[V_LO, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral_og('N3LO_EM500', S[0], L[0], Lprime[0], J[0], T[0], Nrows)
#[V_NLO, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral_og('NLO_EM500new', S[0], L[0], Lprime[0], J[0], T[0], Nrows)
#[V_N2LO, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral_og('N2LO_EM500new', S[0], L[0], Lprime[0], J[0], T[0], Nrows)
#[V_N3LO, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral_og('N3LO_EM500new', S[0], L[0], Lprime[0], J[0], T[0], Nrows)
#[V_N4LO, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral_og('N4LO_EM500new', S[0], L[0], Lprime[0], J[0], T[0], Nrows)






def phase_shift_correct(phase_shifts):
    p = np.copy(phase_shifts)
    found = False
    l = len(phase_shifts)
    for i in range(l - 2):
        if np.abs(phase_shifts[l - i - 2] - phase_shifts[l - i - 1]) > 80:  # previously < 1
            index = l - i - 1
            found = True
            # print(index)
            break
    if found == True:
        for m in range(index):
            phase_shifts[m] = np.copy(phase_shifts[m]) + 180
    return phase_shifts
    #return phase_shifts
'''
def phase_shift_correct(phase_shifts):
    return phase_shifts
'''

#energy = Elab(mesh_points)


orders = ["LO", "NLO", "N2LO", "N3LO", "N4LO"]


### calculates phase shifts for:
# given partial wave ("10010" : 3S1),
# chiral order (0: LO, 1: NLO...),
# SVD_order (0: 1 singular value and operator, 1: 2 singular values and operators...),
# lambda ("2.00" for 2.00fm-1),
# singular values as an array of singular values,
# example function call: SVD('00001', 3, 4, '2.00', [30, 15, 5, 1, 0.5])
# for an N3LO lambda = 2.00, SVD up to order 5, 1S0 potential
def SVD(partial_wave,  SVD_order,  singular_values, kmax):
        potential_sum = 0

        for o in range(SVD_order + 1):
            sv = singular_values[o]

            operator_file = f'_operator_{o}_VNN_N3LO_EM500_SLLJT_{partial_wave}_lambda_1.80_Np_100_np_nocut.dat'
            file_00 = os.path.join(operator_path, operator_file)


            [V00, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral(file_00 , Nrows)
            V00_new = np.copy(V00) * sv
            potential_sum += V00_new

        K = compute_K_matrix(potential_sum , Nrows, kmax, mesh_points, mesh_weights)
        phase_shifts = np.array([compute_phase_shifts(K, i, mesh_points) for i in range(Nrows)])
        phase_shifts = phase_shift_correct(phase_shifts)
        energy = Elab(mesh_points)
        return phase_shifts, energy


def SVD_coupled(partial_wave_00, partial_wave_01, partial_wave_10, partial_wave_11 , SVD_order, singular_values, kmax):
    potential_sum = 0

    for o in range(SVD_order + 1):
        sv = singular_values[o]

        operator_file_00 = f'_operator_{o}_VNN_N3LO_EM500_SLLJT_{partial_wave_00}_lambda_1.80_Np_100_np_nocut.dat'
        file_00 = os.path.join(operator_path, operator_file_00)



        [V00, Vmat_1, Tkin_1, mesh_points, mesh_weights] = read_Vchiral(file_00, Nrows)
        V00_new = np.copy(V00) * sv
        potential_sum += V00_new

    file_01 = f'{potential_path}/SLLJT_{partial_wave_01}/VNN_N3LO_EM500_SLLJT_{partial_wave_01}_lambda_1.80_Np_100_np_nocut.dat'
    file_10 = f'{potential_path}/SLLJT_{partial_wave_10}/VNN_N3LO_EM500_SLLJT_{partial_wave_10}_lambda_1.80_Np_100_np_nocut.dat'
    file_11 = f'{potential_path}/SLLJT_{partial_wave_11}/VNN_N3LO_EM500_SLLJT_{partial_wave_11}_lambda_1.80_Np_100_np_nocut.dat'

    [V01, _, _, _, _] = read_Vchiral(file_01, Nrows)
    [V10, _, _, _, _] = read_Vchiral(file_10, Nrows)
    [V11, _, _, _, _] = read_Vchiral(file_11, Nrows)

    K = compute_K_matrix_coupled(potential_sum,V01 , V10, V11, Nrows, kmax, mesh_points, mesh_weights)
    #phase_shifts = np.array([compute_phase_shifts(K, i, mesh_points) for i in range(Nrows)])
    phase_shifts = np.array([compute_phase_shifts_coupled(K, i, mesh_points)[0] for i in range(Nrows)])
    phase_shifts = phase_shift_correct(phase_shifts)
    energy = Elab(mesh_points)
    return phase_shifts, energy





