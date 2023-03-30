
def _pic_1(x: float) -> float:
    return (x**3) / 3


def _pic_2(x: float) -> float:
    return (x ** 3) / 3 + (x ** 7) / 63


def _pic_3(x: float) -> float:
    return (x ** 3) / 3 + (x ** 7) / 63 + 2 * (x ** 11) / 2079 + (x ** 15) / 59535


def _pic_4(x: float) -> float:
    return (x ** 3) / 3 + (x ** 7) / 63 + 2 * (x ** 11) / 2079 + 13* (x ** 15) / 218295 + 82 * (x ** 19) / 37328445 + \
        662 * (x ** 23) / 10438212015 + 4 * (x ** 27) / 3341878155 + (x ** 31) / 109876902975


def _picard(x_0: float, step: float, iters: int, pic_f) -> list[float]:
    res = list()

    for i in range(iters):
        y_0 = pic_f(x_0)
        res.append(y_0)
        x_0 += step

    return res


def picard_1(x_0: float, step: float, iters: int) -> list[float]:
    return _picard(x_0, step, iters, _pic_1)


def picard_2(x_0: float, step: float, iters: int) -> list[float]:
    return _picard(x_0, step, iters, _pic_2)


def picard_3(x_0: float, step: float, iters: int) -> list[float]:
    return _picard(x_0, step, iters, _pic_3)


def picard_4(x_0: float, step: float, iters: int) -> list[float]:
    return _picard(x_0, step, iters, _pic_4)
