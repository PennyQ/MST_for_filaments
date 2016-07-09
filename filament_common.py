__author__ = 'penny'

import networkx as nx
import matplotlib.pyplot as plt
# import csv
import numpy as np
import matplotlib as mpl
import aplpy
from ds9norm import DS9Normalize
# import matplotlib.patches as mpatches
# import matplotlib.lines as mlines
import math
# import decimal
from astropy.table import Table
from astropy.io import fits
from xlwt import *
import time
import os
from math import *
from get_tree import GraphMST

# TODO: divide draw figures from MST generation/analysis


class DrawFilament():
    """
    This class is to initiate filament bg_figure and set graph nodes,
    which is common for each filament.
    """


    SCUTUM_FINAL = '/Users/penny/Works/MST_filaments/Extinction_filaments_data' \
                   '/Scutum_final.txt'
    SCUTUM_FINAL_4TH_QUAD = '/Users/penny/Works/MST_filaments/Extinction_filaments_data' \
                   '/Scutum_final_4th_quad.txt'

    def __init__(self, fila_n, min_l, max_l, x_box, y_box, line_b):

        self.min_l = min_l
        self.max_l = max_l
        self.x_box = x_box
        self.y_box = y_box

        assert len(line_b) == 3
        self.line_b = line_b
        self.fil_n = fila_n

        # init graph
        # self.graph = self.init_graph()

        # get bg figure
        _fil_candidate_file = '/Users/penny/Works/MST_filaments/' \
                              'Extinction_filaments_data/Candid%d_cropped_final.fits' % self.fil_n
        self.bg_fig = self.get_bg_figure(_fil_candidate_file)

        # init output directory
        if not os.path.exists('./fil%d_output' % self.fil_n):
            os.makedirs('./fil%d_output' % self.fil_n)

    def get_bg_figure(self, fil_candidate_file):
        myfig = plt.figure(figsize=(20, 20))

        # TODO: set the file directory as function parameter
        fig = aplpy.FITSFigure(fil_candidate_file, figure=myfig, zorder=0)
        norm = DS9Normalize(stretch='sqrt', clip_hi=99, clip_lo=1,
                            contrast=1.56, bias=0.65)
        norm.update_clip(fig._data)

        fig.show_grayscale()
        fig.image.set_norm(norm)
        fig.show_grayscale(vmin=0, vmax=145)
        fig.image.set_norm(norm)

        fig.tick_labels.set_font(family='sans-serif', size='large')
        fig.tick_labels.set_xformat("dd")
        fig.tick_labels.set_yformat("d.dd")
        fig.axis_labels.set_font(family='sans-serif',size='large')
        fig.axis_labels.set_xtext("Galactic Longitude [deg]")
        fig.axis_labels.set_ytext("Galactic Latitude [deg]")
        # plt.title("Filament%d Kruskal MST, threshold %.2fdeg" % (self.fil_n, self.threshold))
        fig.axis_labels.set_font(size=20)
        fig.tick_labels.set_font(size=20)

        # the galactic longitude range of the output figure
        if self.fil_n in [1, 2, 3]:
            scutumdata = Table.read(self.SCUTUM_FINAL, format='ascii',
                                    delimiter='\t', guess=False)
        else:
            scutumdata = Table.read(self.SCUTUM_FINAL_4TH_QUAD, format='ascii',
                                    delimiter='\t', guess=False)
        ws = (scutumdata['long_final'] > self.min_l) & (scutumdata['long_final'] < self.max_l)
        arml = scutumdata['long_final'][ws]
        num_rows = arml.size

        cmap = plt.cm.get_cmap("gist_rainbow")
        cmap.set_under(color="k")
        cmap.set_over(color="k")
        norm = mpl.colors.Normalize(clip=False,vmin=65,vmax=85)
        m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

        # Filament coordinates
        upperb = [self.line_b[0]]*num_rows
        armb = [self.line_b[1]]*num_rows
        lowerb = [self.line_b[2]]*num_rows
        # armv = scutumdata['Vlsr'][ws]

        box = [np.vstack((self.x_box, self.y_box))]

        # Q~for draw three parallel lines, why it's called rectangles??
        # it's for setting different cmap
        fig.show_rectangles(arml, armb, 0.1, 0.003, color="magenta", cmap=cmap,  #0.03, 0.002 for fila1
                            norm=norm, edgecolors="none", linewidth=2)
        fig.show_rectangles(arml, upperb, 0.3, 0.001, color="red", cmap=cmap,
                            norm=norm, edgecolors="none", linewidth=2)
        fig.show_rectangles(arml, lowerb, 0.3, 0.001, color="red", cmap=cmap,
                            norm=norm, edgecolors="none", linewidth=2)

        fig.show_lines(box, color="crimson", linewidth=3, zorder=10)

        return fig

    def output_graph(self, threshold):  # threshold = np.linspace(0.04, 0.09, 1)
        """
        For each threshold there is a new graph
        """
        self.bg_fig.remove_layer('lines')
        for i in threshold:
            graph_mst = GraphMST(i, self.max_l, self.min_l, self.fil_n, self.bg_fig)
            graph_mst.save_tree()




