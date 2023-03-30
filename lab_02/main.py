
import prettytable
from config import DEFAULT_F_0, DEFAULT_ITERS, DEFAULT_Z_0, DEFAULT_Z_MAX
from solver import U_P, calc_xi, runge_kutta_4th
import seaborn as sns
import matplotlib.pyplot as plt


def res_log(points, f_points, u_points, u_p_points) -> None:

    table = prettytable.PrettyTable(
        ["No", "z", "F", "u", "u_p"])

    for i in range(len(points)):
        table.add_row([i, points[i], 10*f_points[i], u_points[i], u_p_points[i]])

    f = open("log.txt", "w")
    print(table, file=f)
    f.close()


def res_graph(points, f_points, u_points, u_p_points) -> None:
    fig, axs = plt.subplots(2)
    fig.suptitle("Result")
    axs[0].plot(points, [elem * 10 for elem in f_points], label="F(z)")
    axs[0].legend()
    axs[1].plot(points, u_points, label="u(z)")
    axs[1].plot(points, u_p_points, label="u_p(z)")
    axs[1].legend()
    plt.show()


def main():
    xi = calc_xi(DEFAULT_Z_0, DEFAULT_Z_MAX, DEFAULT_ITERS,
                 DEFAULT_F_0, U_P(DEFAULT_Z_0))
    print("Xi value = ", xi)

    res = runge_kutta_4th(DEFAULT_Z_0, DEFAULT_Z_MAX,
                          DEFAULT_ITERS, DEFAULT_F_0, xi * U_P(DEFAULT_Z_0))
    res.append([U_P(z) for z in res[0]])

    res_log(*res)
    res_graph(*res)


if __name__ == "__main__":
    main()
