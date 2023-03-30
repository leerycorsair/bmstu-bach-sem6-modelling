
import matplotlib.pyplot as plt
import math
import prettytable
from cauchy import euler, runge_kutta, picard


def diffeq(x: float, y: float):
    return x*x + y*y


def res_log(points, euler_sol, runge_kutta_sol, picard_sol_1,
            picard_sol_2, picard_sol_3, picard_sol_4) -> None:

    table = prettytable.PrettyTable(
        ["No", "X", "Euler", "RungeKutta", "Picard 1", "Picard 2", "Picard 3", "Picard 4"])

    for i in range(len(points)):
        table.add_row([i, points[i], euler_sol[i], runge_kutta_sol[i],
                      picard_sol_1[i], picard_sol_2[i], picard_sol_3[i], picard_sol_4[i]])
    table.float_format = ".5"

    f = open("log.txt", "w")
    print(table, file=f)
    f.close()


def res_graph(points, euler_sol, runge_kutta_sol, picard_sol_1,
              picard_sol_2, picard_sol_3, picard_sol_4) -> None:


    def _extend():
        lists = [points, euler_sol, runge_kutta_sol, picard_sol_1,
                 picard_sol_2, picard_sol_3, picard_sol_4]

        for l in lists:
            l.extend([-1 * elem for elem in l])
            l.sort()

    _extend()

    plt.title("Graph")
    plt.xlabel("X - values")
    plt.ylabel("Y - values")
    plt.plot(points, euler_sol, label="Euler")
    plt.plot(points, runge_kutta_sol, label="RungeKutta")
    plt.plot(points, picard_sol_1, label="Picard 1")
    plt.plot(points, picard_sol_2, label="Picard 2")
    plt.plot(points, picard_sol_3, label="Picard 3")
    plt.plot(points, picard_sol_4, label="Picard 4")
    plt.legend()
    plt.show()


def main():
    x_0 = 0
    y_0 = 0

    x_max = float(input('\n\nEnter X_MAX value:'))
    step = float(input('Enter STEP value:'))

    iterations = math.ceil(abs(x_max - x_0) / step) + 1

    points = [(x_0 + i * step) for i in range(iterations)]
    euler_sol = euler.euler(diffeq, x_0, y_0, step, iterations)
    runge_kutta_sol = runge_kutta.runge_kutta(
        diffeq, x_0, y_0, step, iterations)
    picard_sol_1 = picard.picard_1(x_0, step, iterations)
    picard_sol_2 = picard.picard_2(x_0, step, iterations)
    picard_sol_3 = picard.picard_3(x_0, step, iterations)
    picard_sol_4 = picard.picard_4(x_0, step, iterations)

    res_log(points, euler_sol, runge_kutta_sol, picard_sol_1,
            picard_sol_2, picard_sol_3, picard_sol_4)
    res_graph(points, euler_sol, runge_kutta_sol, picard_sol_1,
              picard_sol_2, picard_sol_3, picard_sol_4)


if __name__ == "__main__":
    main()
