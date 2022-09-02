#!/usr/bin/env python

import scipy.interpolate as interpolate
import numpy as np

def no_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind):
    mu_rho = mu_rho**-1.#inverting the data to help the interp
    interp_func = interpolate.interp1d(photon_MeV, mu_rho, kind='cubic')
    interp = interp_func(E)

    mu_rho = mu_rho**-1.
    interp = interp**-1.

    return interp

def one_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind):
    mu_rho = mu_rho**-1.#inverting the data to help the interp
    edge_1 = interp_edge[1]
    booles_1 = (E<np.max(photon_MeV[0:edge_1]))
    inds_1   = np.where(booles_1)
    E_1      = E[inds_1]
    interp_func_1 = interpolate.interp1d(photon_MeV[0:edge_1], mu_rho[0:edge_1], kind=interp_kind[1])
    interp_1 = interp_func_1(E_1)

    booles_2 = (E>np.min(photon_MeV[edge_1:]))
    inds_2   = np.where(booles_2)
    E_2      = E[inds_2]

    interp_func_2 = interpolate.interp1d(photon_MeV[edge_1:], mu_rho[edge_1:], kind='cubic')
    interp_2 = interp_func_2(E_2)

    interp = np.concatenate((interp_1,interp_2))

    mu_rho = mu_rho**-1.
    interp = interp**-1.

    return interp

def two_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind):
    mu_rho = mu_rho**-1.#inverting the data to help the interp
    edge_1 = interp_edge[1]
    booles_1 = (E<np.max(photon_MeV[0:edge_1]))
    inds_1   = np.where(booles_1)
    E_1      = E[inds_1]
    interp_func_1 = interpolate.interp1d(photon_MeV[0:edge_1], mu_rho[0:edge_1], kind=interp_kind[1])
    interp_1 = interp_func_1(E_1)

    edge_2 = interp_edge[2]
    booles_2 = np.logical_and(E>=np.max(photon_MeV[0:edge_1]), E<np.min(photon_MeV[edge_2:]))
    inds_2   = np.where(booles_2)
    E_2      = E[inds_2]
    interp_func_2 = interpolate.interp1d(photon_MeV[edge_1:edge_2], mu_rho[edge_1:edge_2], kind=interp_kind[2])
    interp_2 = interp_func_2(E_2)

    booles_3 = (E>np.min(photon_MeV[edge_2:]))
    inds_3   = np.where(booles_3)
    E_3      = E[inds_3]

    interp_func_3 = interpolate.interp1d(photon_MeV[edge_2:], mu_rho[edge_2:], kind='cubic')
    interp_3 = interp_func_3(E_3)

    interp = np.concatenate((interp_1,interp_2))
    interp = np.concatenate((interp,interp_3))

    mu_rho = mu_rho**-1.
    interp = interp**-1.

    return interp

def four_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind):
    mu_rho = mu_rho**-1.#inverting the data to help the interp
    edge_1 = interp_edge[1]
    booles_1 = (E<np.max(photon_MeV[0:edge_1]))
    inds_1   = np.where(booles_1)
    E_1      = E[inds_1]
    interp_func_1 = interpolate.interp1d(photon_MeV[0:edge_1], mu_rho[0:edge_1], kind=interp_kind[1])
    interp_1 = interp_func_1(E_1)

    edge_2 = interp_edge[2]
    booles_2 = np.logical_and(E>=np.max(photon_MeV[0:edge_1]), E<np.min(photon_MeV[edge_2:]))
    inds_2   = np.where(booles_2)
    E_2      = E[inds_2]
    interp_func_2 = interpolate.interp1d(photon_MeV[edge_1:edge_2], mu_rho[edge_1:edge_2], kind=interp_kind[2])
    interp_2 = interp_func_2(E_2)

    edge_3 = interp_edge[3]
    booles_3 = np.logical_and(E>=np.max(photon_MeV[edge_1:edge_2]), E<np.min(photon_MeV[edge_3:]))
    inds_3   = np.where(booles_3)
    E_3      = E[inds_3]
    interp_func_3 = interpolate.interp1d(photon_MeV[edge_2:edge_3], mu_rho[edge_2:edge_3], kind=interp_kind[3])
    interp_3 = interp_func_3(E_3)

    edge_4 = interp_edge[4]
    booles_4 = np.logical_and(E>=np.max(photon_MeV[edge_2:edge_3]), E<np.min(photon_MeV[edge_4:]))
    inds_4   = np.where(booles_4)
    E_4      = E[inds_4]
    interp_func_4 = interpolate.interp1d(photon_MeV[edge_3:edge_4], mu_rho[edge_3:edge_4], kind=interp_kind[4])
    interp_4 = interp_func_4(E_4)

    booles_5 = (E>np.min(photon_MeV[edge_4:]))
    inds_5   = np.where(booles_5)
    E_5      = E[inds_5]

    interp_func_5 = interpolate.interp1d(photon_MeV[edge_4:], mu_rho[edge_4:], kind='cubic')
    interp_5 = interp_func_5(E_5)

    interp = np.concatenate((interp_1,interp_2))
    interp = np.concatenate((interp,interp_3))
    interp = np.concatenate((interp,interp_4))
    interp = np.concatenate((interp,interp_5))

    mu_rho = mu_rho**-1.
    interp = interp**-1.

    return interp
