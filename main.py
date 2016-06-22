__author__ = 'penny'

from draw_fila import DrawFilament
import numpy as np
import time

# filament pool filament_number: [min_l, max_l, x_box, y_box, line_b]
fil_pool = {1: [25.5, 28.5, [26.8, 27.1, 27.1, 26.8, 26.8], [-0.23, -0.23, -0.4, -0.4, -0.23], [0.06, -0.19, -0.44]],
            2: [22.5, 25.5, [24.9, 25.5, 25.5, 24.9, 24.9], [-0.18, -0.18, -0.6, -0.6, -0.18], [0.06, -0.21, -0.48]],
            3: [22, 26, [24.6, 25.15, 25.15, 24.6, 24.6], [-0.05, -0.05, -0.26, -0.26, -0.05], [0.05, -0.24, -0.53]],
            4: [19.5, 22.5, [21.0, 21.5, 21.5, 21.0, 21.0], [-0.1, -0.1, -0.25, -0.25, -0.1], [0.05, -0.24, -0.53]],
            5: [16.5, 19.5, [18.2, 19.5, 19.5, 18.2, 18.2], [0.14, 0.14, -0.22, -0.22, 0.14], [0.05, -0.26, -0.57]],
            6: [10.5, 13.5, [10.9, 11.4, 11.4, 10.9, 10.9], [0, 0, -0.2, -0.2, 0], [0.04, -0.31, -0.66]],
            9: [334.5, 337.5, [335.0, 335.62, 335.62, 335.0, 335.0], [-0.2, -0.2, -0.4, -0.4, -0.2], [0.03, -0.33, -0.69]],
            10: [331.5, 334.5, [331.5, 333, 333, 331.5, 331.5], [0.1, 0.1, -0.25, -0.25, 0.1], [0.03, -0.32, -0.67]]}

error_log = open('.error_log.txt', 'w')
for fil_i in range(2, 11):
    if fil_i not in fil_pool.keys():
        continue
    for i in np.linspace(0.01, 0.2, 20):  # threshold
        try:
            start = time.time()
            fila = DrawFilament(fila_n=fil_i,
                                min_l=fil_pool[fil_i][0], max_l=fil_pool[fil_i][1],
                                x_box=fil_pool[fil_i][2], y_box=fil_pool[fil_i][3],
                                line_b=fil_pool[fil_i][4], threshold=i)
            fila.save_tree()
            end = time.time()
            print(i, ' time is', end-start)
        except ZeroDivisionError:
            error_log.write('filament%d with threshold %.2f got an error!' % (fil_i, i))
            continue
error_log.close()

