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

# TODO: divide draw figures from MST generation/analysis


class DrawFilament():
    # path to file

    PF_IRDC = '/Users/penny/Works/MST_filaments/Peretto_Fuller_data/peretto_fuller_irdc.fit'
    SCUTUM_FINAL = '/Users/penny/Works/MST_filaments/Extinction_filaments_data' \
                   '/Scutum_final.txt'
    SCUTUM_FINAL_4TH_QUAD = '/Users/penny/Works/MST_filaments/Extinction_filaments_data' \
                   '/Scutum_final_4th_quad.txt'

    def __init__(self, fila_n, min_l, max_l, x_box, y_box, line_b, threshold):
        self.lat = []
        self.lng = []
        self.min_l = min_l
        self.max_l = max_l
        self.graph = nx.Graph()
        self.threshold = threshold

        self.x_box = x_box
        self.y_box = y_box

        assert len(line_b) == 3
        self.line_b = line_b

        # init the graph
        self.create_graph()

        self.fil_n = fila_n

        self.fil_candidate_file = '/Users/penny/Works/MST_filaments/' \
                                  'Extinction_filaments_data/Candid%d_cropped_final.fits' % self.fil_n
        if not os.path.exists('./fil%d_output' % self.fil_n):
            os.makedirs('./fil%d_output' % self.fil_n)

    def create_graph(self):
        # create arrays for storing coordinates of molecular clouds
        hdulist = fits.open(self.PF_IRDC)
        tbdata = hdulist[1].data

        mask = np.where((self.max_l >= tbdata['_Glon']) & (tbdata['_Glon'] >= self.min_l)
                        & (1 >= tbdata['_Glat']) & (tbdata['_Glat'] >= -1))
        covered_data = tbdata[mask]
        self.lat = covered_data['_Glat']
        self.lng = covered_data['_Glon']
        for i in range(len(covered_data)):
            self.graph.add_node(i, posxy=(self.lng[i], self.lat[i]))

        positions = nx.get_node_attributes(self.graph, 'posxy')

        # build edges of graph, within a certain distance treshold
        for n, d in self.graph.nodes_iter(data=True):
            xn, yn = positions[n]
            for m, g in self.graph.nodes_iter(data=True):
                xm, ym = positions[m]
                dist = (math.sqrt(math.pow((xm - xn), 2) + math.pow((ym - yn), 2)))
                if dist <= self.threshold:

                    # tau_av is larger then the opacity is larger, so the weight is lower
                    # IRDC yes then weight is lower -> +0 no ->+1
                    if covered_data['IRDC'][m] is 'y' or covered_data['IRDC'][n] is 'y' or \
                                    covered_data['tau_p'][m] > 1.5 or covered_data['tau_p'][n] > 1.5:
                        has_key_node = True
                    else:
                        has_key_node = False

                    # set attribute for each edge
                    self.graph.add_edge(n, m, weight=dist, has_key_node=has_key_node)

    def get_bg_figure(self):
        myfig = plt.figure(figsize=(20, 20))

        # TODO: set the file directory as function parameter
        fig = aplpy.FITSFigure(self.fil_candidate_file, figure=myfig, zorder=0)
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
        plt.title("Filament%d Kruskal MST, threshold %.2fdeg" % (self.fil_n, self.threshold))
        fig.axis_labels.set_font(size=20)
        fig.tick_labels.set_font(size=20)

        # the galactic longitude range of the output figure
        start = time.time()
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

        # Filament 1 coordinates
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
        fig.show_markers(self.lng, self.lat, s=30, alpha=1, zorder=2,
                         marker='^', c='yellow')
        fig.show_lines(box, color="crimson", linewidth=3, zorder=10)

        return fig

    def get_tree_list(self):
        # using Kruskal's algorithm
        # If the graph is not connected a spanning forest is constructed.
        # A spanning forest is a union of the spanning trees for each connected component of the graph.
        T = nx.minimum_spanning_tree(self.graph)

        # using Prim's algorithm
        # T = nx.prim_mst(self.graph)

        # generate lists of individual trees
        tree_list = nx.connected_component_subgraphs(T)

        # #largest tree in the region
        # max_tree = max(tree_list, key = len)
        return tree_list

    def get_distance(self, delta1, delta2):
        return sqrt(delta1**2 + delta2**2)

    def save_tree(self):
        i = 1
        workbook = Workbook()
        ws = workbook.add_sheet('Filament'+str(self.fil_n))

        # add title for sheet
        title = ['size', 'center_point_xpos', 'center_point_ypos', 'tree_weight', 'edge_ave_weight']
        for col, value in enumerate(title):
            ws.write(0, col, value)

        tree_list = self.get_tree_list()

        # draw MST on top of bg figure
        fig = self.get_bg_figure()

        # analysis each tree
        for each_tree in tree_list:
            # skip one node
            if each_tree.size() == 0:  # number of edges
                continue

            # node centralities: measuring length/ width ratio of tree
            node_pos = nx.get_node_attributes(each_tree, 'posxy')
            # tree_node = np.array(pos_tree.keys())
            tree_posxy = np.array(node_pos.values())
            x_pos = tree_posxy[:, 0]
            y_pos = tree_posxy[:, 1]

            # for storing center point
            delta_x = max(x_pos) - min(x_pos)
            delta_y = max(y_pos) - min(y_pos)
            each_tree_weight = 0

            # draw mst trees here
            start = time.time()
            for each in sorted(each_tree.edges(data='has_key_node')):
                # print('what is each edge', each, each[2])
                # ('what is each edge', (323, 337, 2.2539938916387352), 2.2539938916387352)
                each_tree_weight += int(each[2])

                x1, y1 = node_pos[each[0]]
                x2, y2 = node_pos[each[1]]
                mst_long = [x1, x2]
                mst_lat = [y1, y2]
                edge = [np.vstack((mst_long, mst_lat))]

                # TODO: can I use a new fig as the bg_fig, and draw upon it?
                if each[2]:
                    fig.show_lines(edge, color='deeppink', alpha=0.6, linewidth=4)
                else:
                    fig.show_lines(edge, color='aquamarine', alpha=0.6, linewidth=4)

            ave_weight = float(each_tree_weight)/float(each_tree.size())
            tree = [each_tree.size(), delta_x/2., delta_y/2., each_tree_weight, ave_weight]

            for col, col_value in enumerate(tree):
                ws.write(i, col, col_value)
            i += 1
            end = time.time()
            print('draw mst tree time', end -start)

        start = time.time()
        fig_name = './fil%d_output/Filament%d_%.2f_MST.png' % (self.fil_n, self.fil_n, self.threshold)
        plt.savefig(fig_name)

        end = time.time()
        print('save fig time', end - start)  # TODO: here is 10 sec, almost half of time

        wb_name = './fil%d_output/tree1_%.2f.xls' % (self.fil_n, self.threshold)
        workbook.save(wb_name)
