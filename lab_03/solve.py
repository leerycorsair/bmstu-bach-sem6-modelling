from utils import MathFuncs


class Solve(MathFuncs):
    def k_n(self, z):
        return self.c / (3 * self.R * self.k(z))

    def half_kappa(self, z):
        return (self.k_n(z + self.STEP_RK / 2) + self.k_n(z - self.STEP_RK / 2)) / 2

    def f_n(self, z):
        return self.c * self.k(z) * self.u_p(z)

    def p_n(self, z):
        return self.c * self.k(z)

    def V_n(self, z):
        return ((z + self.STEP_RK / 2) ** 2 - (z - self.STEP_RK / 2) ** 2) / 2

    def V_n_plus(self, z):
        return ((z + self.STEP_RK / 2) ** 2 - z ** 2) / 2

    def V_n_minus(self, z):
        return (z ** 2 - (z - self.STEP_RK / 2) ** 2) / 2

    # Простая аппроксимация
    def approx_plus_half(self, func, n):
        return (func(n) + func(n + self.STEP_RK)) / 2

    def approx_minus_half(self, func, n):
        return (func(n - self.STEP_RK) + func(n)) / 2

    def A(self, z):
        return (z - self.STEP_RK / 2) * (self.half_kappa(z - self.STEP_RK / 2))

    def C(self, z):
        return (z + self.STEP_RK / 2) * self.half_kappa(z + self.STEP_RK / 2)

    def B(self, z):
        return self.A(z) + self.C(z) + self.p_n(
            z) * z * self.STEP_RK ** 2 * self.R

    def D(self, z):
        return self.f_n(z) * z * self.STEP_RK ** 2 * self.R

    def left_boundary_condition(self, z0, F0, h):
        K0 = -self.half_kappa(z0 + h / 2) * (z0 + h / 2) + self.c * self.R * h * h / 8 * self.k(z0 + h / 2) * (
            z0 + h / 2)
        M0 = self.half_kappa(z0 + h / 2) * (z0 + h / 2) + self.c * self.R * h * h / 8 * self.k(z0 + h / 2) * (
            z0 + h / 2)
        P0 = self.c * self.R * h * h / 4 * \
            self.k(z0 + h / 2) * self.u_p(z0 + h / 2) * (z0 + h / 2)
        return K0, M0, P0

    def right_boundary_condition(self, z, h):
        KN = self.half_kappa(z - h / 2) * (
            z - h / 2) + self.m * self.c * z * h / 2 + self.c * self.R * h * h / 8 * self.k(z - h / 2) * (
            z - h / 2) + self.R * self.c * h * h * self.k(z) / 4
        MN = -self.half_kappa(z - h / 2) * (z - h / 2) + self.c * \
            self.R * h * h / 8 * self.k(z - h / 2) * (z - h / 2)
        PN = self.c * self.R * h * h / 4 * (
            self.k(z - h / 2) * self.u_p(z - h / 2) * (z - h / 2) + self.k(z) * self.u_p(z))
        return KN, MN, PN

    def right_hod(self):
        # Прямой ход
        h = self.STEP_RK
        K0, M0, P0 = self.left_boundary_condition(0, 0, self.STEP_RK)
        KN, MN, PN = self.right_boundary_condition(1, self.STEP_RK)
        eps = [0, -K0 / M0]
        eta = [0, P0 / M0]

        x = h
        n = 1
        while x < self.z_max:
            eps.append(self.C(x) / (self.B(x) - self.A(x) * eps[n]))
            eta.append((self.A(x) * eta[n] + self.D(x)) /
                       (self.B(x) - self.A(x) * eps[n]))
            n += 1
            x += h

        # Обратный ход
        u = [0] * (n)
        u[n - 1] = (PN - MN * eta[n]) / (KN + MN * eps[n])

        for i in range(n - 2, -1, -1):
            u[i] = eps[i + 1] * u[i + 1] + eta[i + 1]

        return u

    def get_center(self, y, z, h):
        res = [(-3 * y[0] + 4 * y[1] - y[2]) / 2 / h]
        for i in range(1, len(y) - 1):
            r = (y[i + 1] - y[i - 1]) / 2 / h
            res.append(r)
        res.append((3 * y[-1] - 4 * y[-2] + y[-3]) / 2 / h)
        return res

    def F_res_deriv(self, u, z):
        f = [0]
        u_res = self.get_center(u, z, self.STEP_RK)
        for i in range(1, len(z)):
            r = -self.c / 3 / self.R / self.k(z[i]) * u_res[i]
            f.append(r)
        return f

    def F_res_integ(self, z, un, un1, f):
        if abs(z - 1) < 1e-4:
            return self.m * self.c * un / 2
        return self.half_kappa(z - self.STEP_RK / 2) * (un - un1) / self.STEP_RK
