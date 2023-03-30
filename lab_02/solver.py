
from math import exp
from config import DEFAULT_K_0, DEFAULT_M, DEFAULT_P, DEFAULT_PRECISION, DEFAULT_R, DEFAULT_T_0, DEFAULT_T_W


DEFAULT_C = 3 * 10e8


def T(z: float) -> float:
    return (DEFAULT_T_W - DEFAULT_T_0) * z ** DEFAULT_P + DEFAULT_T_0


def U_P(z: float) -> float:
    return 3.084e-4 / (exp(4.799e4/T(z)) - 1)


def K(z: float) -> float:
    return DEFAULT_K_0 * (T(z)/300)**2


def U(z: float, F: float) -> float:
    return -F * 3 * DEFAULT_R * K(z) / DEFAULT_C


def F(z: float, F: float, u: float) -> float:
    if (abs(z) < DEFAULT_PRECISION):
        return DEFAULT_C * DEFAULT_R / 2 * K(z) * (U_P(z) - u)
    else:
        return DEFAULT_R * DEFAULT_C * K(z) * (U_P(z) - u) - F / z


def Xi(F: float, u: float) -> float:
    return F - DEFAULT_M * DEFAULT_C * u / 2


def runge_kutta_4th(z_start: float, z_end: float, iters: int, f_start: float, u_start: float) -> list[list[float]]:
    res = [[z_start], [f_start], [u_start]]
    z_cur = z_start
    f_cur = f_start
    u_cur = u_start
    h = 1 / (iters - 1)
    while (abs(z_cur - z_end) > 1e-7):
        # k - for f, q - for u
        k1 = h * F(z_cur, f_cur, u_cur)
        q1 = h * U(z_cur, f_cur)

        k2 = h * F(z_cur + h/2, f_cur + k1/2, u_cur + q1/2)
        q2 = h * U(z_cur + h/2, f_cur + k1/2)

        k3 = h * F(z_cur + h/2, f_cur + k2/2, u_cur + q2/2)
        q3 = h * U(z_cur + h/2, f_cur + k2/2)

        k4 = h * F(z_cur + h, f_cur + k3, u_cur + q3)
        q4 = h * U(z_cur + h, f_cur + k3)

        f_cur += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        u_cur += (q1 + 2 * q1 + 2 * q3 + q4) / 6
        z_cur += h

        res[0].append(z_cur)
        res[1].append(f_cur)
        res[2].append(u_cur)
    return res


def calc_xi(z_start: float, z_stop: float, iters: int, f_start: float, up: float):
    xi_l = 0.0
    xi_r = 1.0
    while (abs(xi_l - xi_r) > DEFAULT_PRECISION):
        xi_c = (xi_l + xi_r) / 2
        r_c = runge_kutta_4th(z_start, z_stop, iters, f_start, xi_c * up)
        r_l = runge_kutta_4th(z_start, z_stop, iters, f_start, xi_l * up)

        psi_c = Xi(r_c[1][-1], r_c[2][-1])
        psi_l = Xi(r_l[1][-1], r_l[2][-1])

        if (psi_c * psi_l < 0):
            xi_r = xi_c
        else:
            xi_l = xi_c
    return (xi_l + xi_r) / 2
