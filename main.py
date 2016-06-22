__author__ = 'penny'

from draw_fila import DrawFilament
import numpy as np
import time

# for filament 1

for i in np.linspace(0.01, 0.2, 20):  # threshold
    start = time.time()
    fila = DrawFilament(min_l=25.5, max_l=28.5,
                        x_box=[26.8, 27.1, 27.1, 26.8, 26.8],
                        y_box=[-0.23, -0.23, -0.4, -0.4, -0.23],
                        line_b=[0.06, -0.19, -0.44], threshold=i)
    fila.save_tree()
    end = time.time()
    print(i, ' time is', end-start)

