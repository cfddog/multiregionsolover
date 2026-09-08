#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analytical reference solutions for the Beavers-Joseph validation case.

Model solved by OpenCFD-EC porous/low-speed solvers (SI units):
   free fluid (0<y<h)      : mu u'' = dp/dx
   porous bed (-Hp<y<0)    : mu u'' - (eps^2*mu/K) u = dp/dx
with the Ergun permeability K = dp_p^2 eps^3/[150 (1-eps)^2], so the
effective Darcy permeability entering the bed equation is Keff = K/eps^2,
lambda = sqrt(Keff), kappa = 1/lambda, and the deep-bed Darcy velocity is
uD = -(Keff/mu) dp/dx.

Interface at y=0: continuity of u and of mu du/dy; no-slip at y=h (top
wall of the fluid channel) and at y=-Hp (impermeable bed bottom).

The exact finite-bed solution (composite) is closed form; in the deep-bed
limit kappa*Hp >> 1 it reduces to the classical Beavers-Joseph slip
profile with slip coefficient alpha = 1 (based on Keff):
   du/dy|_{0+} = (u_int - uD)/lambda .

Functions:
   params(eps, dp_p, mu)
   solve_composite(G, eps, dp_p, mu, h, Hp)
       -> dict with uD, lambda, kappa, u_int, C, A, B
   u_fluid(y, sol), u_porous(y, sol)
"""
import math


def params(eps=0.8, dp_p=0.5, mu=0.02):
    K = dp_p**2 * eps**3 / (150.0 * (1.0 - eps)**2)
    Keff = K / eps**2
    lam = math.sqrt(Keff)
    return dict(eps=eps, dp_p=dp_p, mu=mu, K=K, Keff=Keff,
                lam=lam, kappa=1.0 / lam)


def solve_composite(G, eps=0.8, dp_p=0.5, mu=0.02, h=1.0, Hp=1.2):
    """G = -dp/dx > 0 (pressure-driven flow toward +x).  Returns dict."""
    p = params(eps, dp_p, mu)
    Keff, lam, kappa = p['Keff'], p['lam'], p['kappa']
    uD = G * Keff / mu
    ch = math.cosh(kappa * Hp)
    sh = math.sinh(kappa * Hp)
    tanh = sh / ch
    sech = 1.0 / ch
    # ui [1 + tanh/(kappa*h)] = uD (1-sech) + G*h*tanh/(2*mu*kappa)
    ui = (uD * (1.0 - sech) + G * h * tanh / (2.0 * mu * kappa)) / \
         (1.0 + tanh / (kappa * h))
    C = (G * h * h / (2.0 * mu) - ui) / h     # du_f/dy(0)
    B = C / kappa                              # bed sinh coefficient
    A = (B * sh - uD) / ch                     # bed cosh coefficient
    sol = dict(p=p, G=G, uD=uD, lam=lam, kappa=kappa,
               u_int=ui, C=C, A=A, B=B, h=h, Hp=Hp, mu=mu)
    return sol


def u_fluid(y, sol):
    """Free-fluid channel profile, 0 <= y <= h."""
    G, mu, ui, C = sol['G'], sol['mu'], sol['u_int'], sol['C']
    return ui + C * y - 0.5 * G * y * y / mu


def u_porous(y, sol):
    """Porous bed profile, -Hp <= y <= 0 (Brinkman layer over Darcy)."""
    uD, A, B, kappa = sol['uD'], sol['A'], sol['B'], sol['kappa']
    return uD + A * math.cosh(kappa * y) + B * math.sinh(kappa * y)


def profile(y_pts, sol):
    """Evaluate the composite analytic u(y) on arbitrary y points."""
    import numpy as np
    y = np.asarray(y_pts, dtype=float)
    u = np.empty_like(y)
    for n, yy in enumerate(y):
        u[n] = u_fluid(yy, sol) if yy >= 0.0 else u_porous(yy, sol)
    return u


def solve_deepbed(G, eps=0.8, dp_p=0.5, mu=0.02, h=1.0):
    """Classical Beavers-Joseph (deep-bed, alpha=1) fluid-gap solution.

    Semi-infinite Darcy bed: du/dy|0 = (u_int-uD)/lambda  =>  slip coefficient
    alpha = 1 (based on Keff).  Returns dict with uD, lam, kappa, u_int, C.
    """
    p = params(eps, dp_p, mu)
    Keff, lam, kappa = p['Keff'], p['lam'], p['kappa']
    uD = G * Keff / mu
    ui = (uD + G * h / (2.0 * mu * kappa)) / (1.0 + 1.0 / (kappa * h))
    C = kappa * (ui - uD)
    return dict(G=G, uD=uD, lam=lam, kappa=kappa, u_int=ui, C=C,
                h=h, mu=mu)


def u_fluid_deepbed(y, sol):
    G, mu, ui, C = sol['G'], sol['mu'], sol['u_int'], sol['C']
    return ui + C * y - 0.5 * G * y * y / mu


def u_poiseuille_noslip(y, G, mu, h):
    """Reference no-slip Poiseuille (both walls solid) for contrast."""
    return 0.5 * G * y * (h - y) / mu


if __name__ == '__main__':
    s = solve_composite(4.1667e-4)
    print('uD=%.6e u_int=%.6e C=%.6e' % (s['uD'], s['u_int'], s['C']))
    print('kappa*Hp=%.2f  kappa*h=%.2f' % (s['kappa'] * s['Hp'],
                                           s['kappa'] * s['h']))
