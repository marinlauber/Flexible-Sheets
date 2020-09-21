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


def build_pressure_mat(beta, sigma, dx, show=False):
    # lower diagonals are just the beta coefficients
    L = sparse.diags([beta[1:],0,beta[1:]],[-1,0,1],shape=(len(sigma),len(sigma))).toarray()
    L[0,-1] = beta[0]; L[-1,0] = beta[-1]
    
    # diag = -sum(non-diags)
    for i in range(len(sigma)):
        L[i, i] = -(np.sum(L[i, :i])+np.sum(L[i, i:]))
    
    # sclae by grid size
    L *= dx**(-2)

    if show:
        p = plt.imshow(abs(L),cmap='binary')
        plt.colorbar(p)
        plt.show()

    return L


def solve_pressure(beta, sigma, dx=1, tol=1e-10, Nsteps=0, verbose=False):

    # build pressure mat
    L = build_pressure_mat(beta, sigma, dx)
    
    # solve using jacobi
    p = Jacobi(L, sigma, tol, Nsteps, verbose)

    return p - p.mean() # zero mean


def Jacobi(A, b, tol, Nsteps, verbose):
    
    # storage arrays
    x = np.zeros_like(b)
    x_n = np.zeros_like(b)
    
    # inverse manually
    inv_ii = inv(A)

    # residuals for convergence
    res0 = residuals(A, x, b)
    res = res0; k=0; step=True
    r = res0*tol
    # solve
    while((res>r) & step):
        
        for i in range(A.shape[0]):
            s = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x_n[i] = (b[i] - s) * inv_ii[i]
        x = x_n

        res = residuals(A, x, b); k+=1
        if((k>=Nsteps) & (Nsteps!=0)): step=False

    if verbose:
        print('Jacobi solver:')
        print('\tres0: %.3e\n' % res0, '\tres: %.3e\n' % res, '\titer: %d' % k)
    return x


def residuals(A, x, b):
    return np.sqrt((len(b))**(-2)*np.einsum('i->', (np.abs(np.matmul(A, x) - b))**2))


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


def kernel(d, eps=2):
    return np.where(abs(d)<=eps, 0.25*(1+np.cos(np.pi*d/eps)), 0)


def draw_piston(plt, X, t=1, V=1):
    plt.fill([X-t/2,X-t/2,X+t/2,X+t/2],[-2,2,2,-2],fill=False,hatch='///')
    plt.quiver([X,X],[1.25,-1.25],[V,V],[0.,0.])


def draw_results(x, xs, X, u0, u, p, sigma, V=1, fname='None'):
    plt.plot(x, u0, '-.k', lw=1, label=r"$u^n$")
    plt.plot(xs, sigma, '-.', lw=1, label=r"div($u^n$)")
    plt.plot(xs, p, '--ok',  lw=1, label=r"$p^n$")
    plt.plot(x, u, '-m', lw=1, label=r"$u^{n+1}$")
    draw_piston(plt, X, 0.015, 1)
    plt.xlim(-1,1); plt.ylim(-2,2)
    plt.xlabel(r'$x/L$')
    plt.legend()
    if fname!='None':
        plt.savefig(fname, dpi=900)