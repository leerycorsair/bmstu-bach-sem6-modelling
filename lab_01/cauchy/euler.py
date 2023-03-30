
def euler(diffeq, x_0: float, y_0: float, step: float, iters: int) -> list[float]:
    res = list()

    for i in range(iters):
        res.append(y_0)
        y_0 += step * diffeq(x_0, y_0)
        x_0 += step

    return res
