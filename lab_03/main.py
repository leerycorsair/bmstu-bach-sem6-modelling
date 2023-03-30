from matplotlib import pyplot as plt
from numpy import arange
from prettytable import PrettyTable
import solve


class Show(solve.Solve):
    def draw(self):
        a = solve.Solve()
        # График
        u_res = a.right_hod()
        z_res = [i for i in arange(0, 1 + self.STEP_RK, a.STEP_RK)]

        f_res = [0] * len(z_res)
        up_res = [0] * len(z_res)
        divF = [0] * len(z_res)

        f3_res = self.F_res_deriv(u_res, z_res)

        for i in range(0, len(z_res) - 1):
            up_res[i] = self.u_p(z_res[i])
            divF[i] = self.divF(z_res[i], u_res[i])

        for i in range(1, len(z_res)):
            f_res[i] = self.F_res_integ(
                z_res[i], u_res[i - 1], u_res[i], f_res[i - 1])

        tb = PrettyTable()
        tb.add_column("Z", z_res)
        tb.add_column("F", f_res)
        tb.add_column("F deriv", f3_res)
        tb.add_column("U", u_res)
        tb.add_column("divF", divF)

        with open('result.txt', 'w') as f:
            f.write(str(tb))

        plt.figure(100)
        plt.plot(z_res, u_res,  label='u')
        plt.plot(z_res, up_res, label='u_p')
        plt.legend()

        plt.figure(200)
        plt.plot(z_res, f_res)
        plt.title('F(z)')

        plt.figure(300)
        plt.plot(z_res, divF)
        plt.title("divF")

        plt.figure(400)
        plt.plot(z_res, f3_res)
        plt.title("F(z) deriv")
        plt.show()


def main():
    res = Show()
    res.draw()


if __name__ == "__main__":
    main()
