# Solvers

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse.linalg import bicgstab
import random




class Solver:

    def __init__(self, model, N, tmax, timestep):
        self.model = model
        self.N = N
        self.tmax = tmax
        self.timestep = timestep
        self.h = np.abs(self.model.x0 - self.model.xf) / (self.N - 1)

        # The Laplacian and the sparse unitary matrix
        N = self.N
        e = np.ones(N ** 2)
        e2 = ([1] * (N - 1) + [0]) * N
        e3 = ([0] + [1] * (N - 1)) * N
        e4 = ([1] + [0] * (N - 1)) * N
        e5 = ([0] * (N - 1) + [1]) * N
        self.M = spdiags([-4 * e, e2, e3, e, e, e4, e5, e, e],
                         [0, -1, 1, -N, N, -(N - 1), N - 1, -(N * N - N), (N * N - N)], N ** 2, N ** 2)
        self.Id = spdiags([e], [0], N ** 2, N ** 2)

        def integrate(self):
            pass

        # Integrate automatically when the class is called
        self.integrate()

    # Initial values
    def generateGrid(self):
        np.random.seed(1234)
        self.X, self.Y = np.meshgrid(np.linspace(self.model.x0, self.model.xf, self.N),
                                     np.linspace(self.model.y0, self.model.yf, self.N))
        N2 = self.N // 2
        r = int(self.N / 20.0)
        self.U = self.model.u0 * np.ones((self.N, self.N))
        self.V = self.model.v0 * np.ones((self.N, self.N))

        self.U[N2 - r:N2 + r, N2 - r:N2 + r] = self.model.up
        self.V[N2 - r:N2 + r, N2 - r:N2 + r] = self.model.vp

        self.U = self.U + np.multiply(self.U, np.random.choice([-self.model.pert, self.model.pert], (self.N, self.N)))
        self.V = self.V + np.multiply(self.V, np.random.choice([-self.model.pert, self.model.pert], (self.N, self.N)))

    # Initial values Giraffe
    def generateGrid_G(self):
        np.random.seed(1234)
        self.X, self.Y = np.meshgrid(np.linspace(self.model.x0, self.model.xf, self.N),
                                     np.linspace(self.model.y0, self.model.yf, self.N))
        N2 = self.N // 2
        r = int(self.N / 20.0)
        self.U = self.model.u0 * np.ones((self.N, self.N))
        self.V = self.model.v0 * np.ones((self.N, self.N))
        self.y_grid = self.model.y_0 * np.ones((self.N, self.N))

        for k in range(0, N2):
            i = np.random.randint(0, self.N - 1)
            j = np.random.randint(0, self.N - 1)
            self.U[i, j] = self.model.up
            self.y_grid[i, j] = 1.0

        # Plotting

    def configPlot(self):
        u = self.U
        v = self.V
        f = plt.figure(figsize=(25, 10), dpi=400, facecolor='w', edgecolor='k');
        sp = f.add_subplot(1, 2, 1);
        plt.pcolor(u.reshape((self.N, self.N)), cmap=plt.cm.copper)
        plt.axis('tight')
        # RdBu
        sp = f.add_subplot(1, 2, 2);
        plt.pcolor(v.reshape((self.N, self.N)), cmap=plt.cm.copper)
        plt.axis('tight')
        plt.show()


class Euler(Solver):
    def __init__(self, model, N, tmax, timestep):
        super(Euler, self).__init__(model, N, tmax, timestep)

    def integrate(self):
        # Short references
        N = self.N
        d1 = self.model.const1
        d2 = self.model.const2
        dt = self.timestep
        dh = self.h
        # Initial values
        self.generateGrid()
        u = self.U.reshape((N * N))
        v = self.V.reshape((N * N))
        # Nr of timesteps
        Nt = int(self.tmax / dt)

        for i in range(Nt):
            u += d1 * (self.M).dot(u) * dt / (dh ** 2) + self.model.f(u, v) * dt
            v += d2 * (self.M).dot(v) * dt / (dh ** 2) + self.model.g(u, v) * dt

        self.U = u
        self.V = v
        return


class Euler1(Solver):
    def __init__(self, model, N, tmax, timestep):
        super(Euler1, self).__init__(model, N, tmax, timestep)

    def integrate(self):
        # Short references
        N = self.N
        d1 = self.model.const1
        d2 = self.model.const2
        dt = self.timestep
        dh = self.h
        # Initial values
        self.generateGrid_G()
        u = self.U.reshape((N * N))
        v = self.V.reshape((N * N))
        y = self.y_grid.reshape((N * N))
        # Nr of timesteps
        Nt = int(self.tmax / dt)

        for i in range(Nt):
            u += d1 * (self.M).dot(u) * dt / (dh ** 2) + self.model.f(u, v) * dt

            v += d2 * (self.M).dot(v) * dt / (dh ** 2) + self.model.g(u, v, y) * dt
            y += self.model.yt(y, u) * dt
            # print(y)
            # print(y0)
        # print(u)
        self.U = u
        self.V = v
        return


class Crank_Nicolson(Solver):
    def __init__(self, model, N, tmax, timestep):
        super(Crank_Nicolson, self).__init__(model, N, tmax, timestep)

    def integrate(self):
        # Short references
        N = self.N
        d1 = self.model.const1
        d2 = self.model.const2
        dt = self.timestep
        dh = self.h
        # Initial values
        self.generateGrid()
        u = self.U.reshape((N * N))
        v = self.V.reshape((N * N))
        # Nr of timesteps
        Nt = round(self.tmax / dt)

        for i in range(int(Nt)):
            RHS1 = u + d1 * (self.M).dot(u) * dt / (2 * dh ** 2) + self.model.f(u, v) * dt / 2
            RHS2 = v + d2 * (self.M).dot(v) * dt / (2 * dh ** 2) + self.model.g(u, v) * dt / 2

            tmp1 = RHS1 + self.model.f(RHS1, RHS2) * dt / 2
            tmp2 = RHS2 + self.model.g(RHS1, RHS2) * dt / 2

            u_f = bicgstab((self.Id - d1 * (self.M) * dt / (2 * dh ** 2)), tmp1, x0=u, tol=1e-5, maxiter=None, M=None,
                           callback=None, atol=None)
            v_f = bicgstab((self.Id - d2 * (self.M) * dt / (2 * dh ** 2)), tmp2, x0=v, tol=1e-5, maxiter=None, M=None,
                           callback=None, atol=None)

            u = u_f[0]
            v = v_f[0]

        self.U = u
        self.V = v
        return


class Crank_NicolsonG(Solver):
    def __init__(self, model, N, tmax, timestep):
        super(Crank_NicolsonG, self).__init__(model, N, tmax, timestep)

    def integrate(self):
        # Short references
        N = self.N
        d1 = self.model.const1
        d2 = self.model.const2
        dt = self.timestep
        dh = self.h
        # Initial values
        self.generateGrid_G()
        u = self.U.reshape((N * N))
        v = self.V.reshape((N * N))
        y = self.y_grid.reshape((N * N))
        # Nr of timesteps
        Nt = int(self.tmax / dt)

        for i in range(Nt):
            RHS1 = u + d1 * (self.M).dot(u) * dt / (2 * dh ** 2) + self.model.f(u, v) * dt / 2
            RHS2 = v + d2 * (self.M).dot(v) * dt / (2 * dh ** 2) + self.model.g(u, v, y) * dt / 2

            tmp1 = RHS1 + self.model.f(RHS1, RHS2) * dt / 2
            tmp2 = RHS2 + self.model.g(RHS1, RHS2, y) * dt / 2

            u_f = bicgstab((self.Id - d1 * (self.M) * dt / (2 * dh ** 2)), tmp1, x0=u, tol=1e-5, maxiter=None, M=None,
                           callback=None, atol=None)
            v_f = bicgstab((self.Id - d2 * (self.M) * dt / (2 * dh ** 2)), tmp2, x0=v, tol=1e-5, maxiter=None, M=None,
                           callback=None, atol=None)

            y += self.model.yt(y, u) * dt

            u = u_f[0]
            v = v_f[0]

        self.U = u
        self.V = v
        return
