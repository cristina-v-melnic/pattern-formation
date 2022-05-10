# Models
class GrayScott:
    def __init__(self, alpha, beta):
        # self.alpha = 0.022
        # self.beta = 0.061

        # Model constants
        self.alpha = alpha
        self.beta = beta

        # Diffusion coefficients
        self.const1 = 0.00001335
        self.const2 = 0.00000572

        # Initial values
        self.u0 = 1.0
        self.v0 = 0
        self.up = 0.50
        self.vp = 0.25

        # Domain size
        self.x0 = 0
        self.xf = 2.5
        self.y0 = 0
        self.yf = 2.5

        # Perturbation
        self.pert = 0.01

    def f(self, u, v):
        return -u * v * v + self.alpha * (1 - u)

    def g(self, u, v):
        return u * v * v - (self.alpha + self.beta) * v


class GrayScott_Pearson:
    def __init__(self, alpha, beta):
        # self.alpha = 0.022
        # self.beta = 0.061

        # Model constants
        self.alpha = alpha
        self.beta = beta

        # Diffusion coefficients
        self.const1 = 0.00002
        self.const2 = 0.00001

        # Initial values
        self.u0 = 1.0
        self.v0 = 0
        self.up = 0.50
        self.vp = 0.25

        # Domain size
        self.x0 = 0
        self.xf = 2.5
        self.y0 = 0
        self.yf = 2.5

        # Perturbation
        self.pert = 0.01

    def f(self, u, v):
        return -u * v * v + self.alpha * (1 - u)

    def g(self, u, v):
        return u * v * v - (self.alpha + self.beta) * v


class Giraffe:
    def __init__(self):
        # Model constants
        self.Ka = 0.1
        self.Ks = 20.0
        self.Ky = 22.0
        self.Rho_a = 0.025
        self.Rho_s = 0.0025
        self.Rho_y = 0.03
        self.sigma_s = 0.00225
        self.sigma_y = 0.00015
        self.mu_s = 0.00075
        self.mu_y = 0.003

        # Diffusion coefficients
        self.const1 = 0.015
        self.const2 = 0.03

        # Initial values
        self.u0 = 0
        self.v0 = 3
        self.up = 5
        self.vp = 3
        self.y_0 = 0

        # Domain size
        self.x0 = 0
        self.xf = 120
        self.y0 = 0
        self.yf = 120
        # Perturbation
        self.pert = 0.25

    def f(self, a, s):
        return self.Rho_a * ((a ** 2 * s) / (1 + self.Ka * a ** 2) - a)

    def g(self, a, s, y):
        return self.sigma_s / (1 + self.Ks * y) - (self.Rho_s * a ** 2 * s) / (1 + self.Ka * a ** 2) - self.mu_s * s

    def yt(self, y, a):
        return self.Rho_y * (y ** 2 / (1 + self.Ky * y ** 2)) - self.mu_y * y + self.sigma_y * a


class Leopard:
    def __init__(self):
        # Model constants
        self.Ka = 0.5
        self.Ks = 0.3
        self.Ky = 22.0
        self.Rho_a = 0.05
        self.Rho_s = 0.0035
        self.Rho_y = 0.03
        self.sigma_s = 0.0075
        self.sigma_y = 0.00007
        self.mu_s = 0.003
        self.mu_y = 0.003

        # Diffusion coefficients
        self.const1 = 0.01
        self.const2 = 0.1

        # Initial values
        self.u0 = 0
        self.v0 = 2.5
        self.up = 2
        self.vp = 2.5
        self.y_0 = 0

        # Domain size
        self.x0 = 0
        self.xf = 150
        self.y0 = 0
        self.yf = 150

        # Perturbation
        self.pert = 0.25

    def f(self, a, s):
        return self.Rho_a * ((a ** 2 * s) / (1 + self.Ka * a ** 2) - a)

    def g(self, a, s, y):
        return self.sigma_s / (1 + self.Ks * y) - (self.Rho_s * a ** 2 * s) / (1 + self.Ka * a ** 2) - self.mu_s * s

    def yt(self, y, a):
        return self.Rho_y * (y ** 2 / (1 + self.Ky * y ** 2)) - self.mu_y * y + self.sigma_y * a

