__author__ = 'penny'

from draw_all_data import DrawAllData
import numpy as np
import time

# TODO: figure the time bottleneck, be either draw_bg_figure or analysis tree

error_log = open('.error_log.txt', 'w')
for lon_range in np.linspace(12, 60, 17):
# for fil_i in range(1, 10):  # test with filament1
    # get bg figure draw
    # for i in np.linspace(0.01, 0.2, 20):  # threshold
    for i in np.linspace(0.04, 0.09, 6):  # threshold
        try:
            start = time.time()
            all_data = DrawAllData(lon_range=lon_range, threshold=i)
            all_data.save_tree()
            end = time.time()
            print(i, 'total time is', end-start)
        except ZeroDivisionError:
            error_log.write('GLM_0%d with threshold %.2f got an error!' % (lon_range, i))
            continue
error_log.close()

