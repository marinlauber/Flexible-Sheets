import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse, optimize, linalg

def r(u, dx=1, nu=0.000001):
    
    # convection & diffusion
    un = np.empty_like(u)
    un[1:-1] = -1*u[1:-1]*(u[2:] - u[:-2])/(2*dx) + nu*((u[:-2] -2*u[1:-1] + u[2:])/dx**2)
    
    # periodic BC's
    un[0] = -1*u[0]*(u[1] - u[-1])/(2*dx) + nu*((u[-1] -2*u[0] + u[1])/dx**2)
    un[-1] = -1*u[-1]*(u[0] - u[-2])/(2*dx) + nu*((u[0] -2*u[-1] + u[-2])/dx**2)
    
    return un


def div(v, dx=1):
    
    divv = np.empty_like(v)
    # out - in / size
    divv[:-1] = (v[1:] - v[:-1]) / dx
    divv[ -1] = (v[ 0] - v[ -1]) / dx
    
    return divv


def grad(f, dx=1):
    
    grad = np.empty(len(f))
    # dfdx|x = f(x+1/2) - f(x-1/2) / h
    grad[1:] = (f[1:] - f[:-1]) / dx
    grad[ 0] = (f[ 0] - f[ -1]) / dx
    
    return grad


def build_pressure_mat(beta, sigma, dx):
    # lower diagonals are just the beta coefficients
    L = sparse.diags([beta[1:],0,beta[1:]],[-1,0,1],shape=(len(sigma),len(sigma))).toarray()
    L[0,-1] = beta[0]; L[-1,0] = beta[-1]
    
    # diag = -sum(non-diags)
    for i in range(len(sigma)):
        L[i, i] = -(np.sum(L[i, :i])+np.sum(L[i, i:]))
    
    # sclae by grid size
    L *= dx**(-2)

    return L


def solve_pressure(beta, sigma, dx=1, Np=5000):

    # biuld pressure mat
    L = build_pressure_mat(beta, sigma, dx)
    
    # solve using jacobi
    p = Jacobi(L, sigma, Np)
    
    return p - p.mean() # zero mean


def Jacobi(A, b, Nint):
    
    # storage arrays
    x = np.zeros_like(b)
    x_n = np.zeros_like(b)
    
    # inverse manually
    inv_ii = inv(A)
    
    # solve
    for _ in range(Nint):
        
        for i in range(A.shape[0]):
            s = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x_n[i] = (b[i] - s) * inv_ii[i]
        x = x_n
    
    return x


def inv(A):
    Ai = np.diag(A).copy()
    for i in range(len(Ai)):
        if(abs(A[i, i])>1e-9):
            Ai[i] = 1./A[i, i]
    return Ai


def ddn(f, dx):
    
    dn = np.empty(len(f))
    # dfdn|x = f(x+1) - f(x-1) / 2h
    dn[1:-1] = (f[2:] - f[:-2]) / dx
    dn[ 0] = (f[1] - f[-1]) / dx
    dn[-1] = (f[0] - f[-2]) / dx
    
    return .5*dn #double dx
