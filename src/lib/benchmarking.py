import time


class Time:
    """
    Calculates the time inside a function. It is necessary to embed it
    inside the code as in the example.

    Example
    -------
    with Timer('fill S2'):
        for i in range(0, N):
            for j in range(0, N):
                for k in range(0, N):
                    print(i, j, k)
    """

    def __init__(self, name=''):
        self.name = name

    def __enter__(self):
        """ Runs when the <with> clause is used """
        self.tstart = time.perf_counter()

    def __exit__(self, type, value, traceback):
        """ Runs when the <with> clause is finished """
        print('%s %.6f s' % (self.name, time.perf_counter() - self.tstart))
