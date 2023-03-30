class MathFuncs:
    def __init__(self):
        # константы неизменяемые
        self.k_0 = 0.008  # 0.018
        self.m = 0.786
        self.R = 0.35
        self.T_w = 2000
        self.T_0 = 10000
        self.c = 3e10
        self._p = 4
        self.STEP_RK = 1e-4
        self.z_max = 1
        self.z_min = 0
        self.EPS = 1e-6
        self.e = 2.718281828459045

    def T(self, z):
        return (self.T_w - self.T_0) * (z ** self._p) + self.T_0

    def k(self, z):
        return self.k_0 * ((self.T(z) / 300) ** 2)

    def u_p(self, z):
        return 3.084e-4 / (self.e ** (4.799e+4 / self.T(z)) - 1)

    def U_z(self, z, f):
        return -(3 * self.R * f * self.k(z)) / self.c

    def F_z(self, z, f, u):
        if abs(z - 0) < 1e-4:
            return ((self.R * self.c) / 2) * self.k(z) * (self.u_p(z) - u)
        else:
            return self.R * self.c * self.k(z) * (self.u_p(z) - u) - (f / z)

    def divF(self, z, u):
        return self.c * self.k(z) * (self.u_p(z) - u)
