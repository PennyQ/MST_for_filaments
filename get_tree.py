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


class GraphMST():
    PF_IRDC = '/Users/penny/Works/MST_filaments/Peretto_Fuller_data/peretto_fuller_irdc.fit'

    def __init__(self, threshold, max_l, min_l, fil_n, bg_fig):
        self.threshold = threshold
        self.lat = []
        self.lng = []
        self.max_l = max_l
        self.min_l = min_l
        self.fil_n = fil_n

        self.bg_fig = bg_fig

        # init graph
        self.graph = self.init_graph()

    def init_graph(self):
        graph = nx.Graph()

        # create arrays for storing coordinates of molecular clouds
        hdulist = fits.open(self.PF_IRDC)
        tbdata = hdulist[1].data

        mask = np.where((self.max_l >= tbdata['_Glon']) & (tbdata['_Glon'] >= self.min_l)
                        & (1 >= tbdata['_Glat']) & (tbdata['_Glat'] >= -1))
        covered_data = tbdata[mask]
        self.lat = covered_data['_Glat']
        self.lng = covered_data['_Glon']
        for i in range(len(covered_data)):
            graph.add_node(i, posxy=(self.lng[i], self.lat[i]))
        positions = nx.get_node_attributes(graph, 'posxy')

        # build edges of graph, within a certain distance treshold
        for n, d in graph.nodes_iter(data=True):
            xn, yn = positions[n]
            for m, g in graph.nodes_iter(data=True):
                xm, ym = positions[m]
                dist = (math.sqrt(math.pow((xm - xn), 2) + math.pow((ym - yn), 2)))
                if dist <= self.threshold:   # dist is the threshold here
                    # tau_av is larger then the opacity is larger, so the weight is lower
                    # IRDC yes then weight is lower -> +0 no ->+1
                    is_irdc = 2.
                    if covered_data['IRDC'][n] is 'y':
                        is_irdc -= 1.
                    if covered_data['IRDC'][m] is 'y':
                        is_irdc -= 1.
                    # TODO: choose a more meaningful way of defining weight here
                    weight = dist*20. - (covered_data['tau_p'][n] + covered_data['tau_p'][m])/2. + is_irdc
                    if weight < 0:
                        weight = 0.
                    # print('new weight and old', weight, dist, is_irdc)
                    graph.add_edge(n, m, weight=weight)
        return graph

    def get_tree_list(self):
        # using Kruskal's algorithm
        # If the graph is not connected a spanning forest is constructed.
        # A spanning forest is a union of the spanning trees for each connected component of the graph.
        T = nx.minimum_spanning_tree(self.graph)

        # using Prim's algorithm
        # T = nx.prim_mst(graph)

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
        title = ['size', 'center_point_xpos', 'center_point_ypos', 'length_ratio', 'likelihood']
        for col, value in enumerate(title):
            ws.write(0, col, value)

        start = time.time()
        tree_list = self.get_tree_list()
        end = time.time()
        print('get_tree_list', end-start)  # 0.014

        # analysis each tree
        for each_tree in tree_list:
            # skip one node
            if each_tree.size() == 0:  # number of edges
                continue

            # if it's parallel
            # node centralities: measuring length/ width ratio of tree
            pos_tree = nx.get_node_attributes(each_tree, 'posxy')
            x_val = []
            y_val = []
            # TODO: check networkx to find a better get access to all nodes
            for n, d in each_tree.nodes_iter(data=True):
                xn, yn = pos_tree[n]
                x_val.append(xn)
                y_val.append(yn)
            x_max = max(x_val)
            x_min = min(x_val)
            y_max = max(y_val)
            y_min = min(y_val)

            delta_x = x_max - x_min
            delta_y = y_max - y_min

            # draw mst trees here
            for each in sorted(each_tree.edges(data='weight')):
                # print('what is each edge', each, each[2])
                # ('what is each edge', (323, 337, 2.2539938916387352), 2.2539938916387352)

                x1, y1 = pos_tree[each[0]]
                x2, y2 = pos_tree[each[1]]
                mst_long = [x1, x2]
                mst_lat = [y1, y2]
                edge = [np.vstack((mst_long, mst_lat))]

                # show markers of nodes
                self.bg_fig.show_markers(self.lng, self.lat, s=30, alpha=1, zorder=2, marker='^', c='yellow')

                # TODO: can I use a new self.bg_fig as the bg_fig, and draw upon it?
                if each[2] > 2:
                    # TODO: show a box around this edge and its neighbour edge
                    self.bg_fig.show_lines(edge, color='blue', alpha=0.6, linewidth=4)
                elif each[2] > 1:
                    self.bg_fig.show_lines(edge, color='aquamarine', alpha=0.6, linewidth=4)
                else:
                    # these are what we want
                    self.bg_fig.show_lines(edge, color='yellow', alpha=0.6, linewidth=4)
                print('self.bg_fig after filament layers', self.bg_fig.list_layers())

            # draw box, avoid wide/high box
            if x_max-x_min < 0.04:
                fil_x_box = [x_min-0.02, x_max+0.02, x_max+0.02, x_min-0.02, x_min-0.02]
                fil_y_box = [y_max, y_max, y_min, y_min, y_max]
            elif y_max-y_min < 0.04:
                fil_x_box = [x_min, x_max, x_max, x_min, x_min]
                fil_y_box = [y_max+0.02, y_max+0.02, y_min-0.02, y_min-0.02, y_max+0.02]
            else:
                fil_x_box = [x_min, x_max, x_max, x_min, x_min]
                fil_y_box = [y_max, y_max, y_min, y_min, y_max]

            box = [np.vstack((fil_x_box, fil_y_box))]

            # check the filament or not ====================

            # break complex structures
            # TODO: modify the threshold later and do a low_threshold MST algorithm for it
            if each_tree.size() > 20:
                continue

            # skip tiny pieces
            # TODO:  modify the threshold later and do a high_threshold MST algorithm for it
            if self.get_distance(delta_x, delta_y) < 0.05:
                continue

            # length_ratio = atan2(delta_y, delta_x)  # return between -pi and pi
            length_ratio = atan(delta_y/delta_x)

            # if length_ratioo ~ 0 -> parallel   -> PI/9 as a threshold
            if fabs(length_ratio) < pi/9 or fabs(length_ratio) > (8./9.)*pi:
                # parallel and skinny, show as yellow
                likelihood = 1.
                # self.bg_fig.show_lines(box, color="yellow", linewidth=2, alpha=0.8, zorder=10)

            #  elif (11./18.)*pi > fabs(length_ratio) > (9./18.)*pi:  # vertical and skinny
            elif fabs(length_ratio) > (7./18.)*pi:  # vertical and skinny
                # self.bg_fig.show_lines(box, color="blue", linewidth=2, alpha=0.8, zorder=10)
                likelihood = 0.5

            else:
                angle = []
                for line in sorted(each_tree.edges(data=True)):  # sorted start from 1st dim
                    x1, y1 = pos_tree[line[0]]
                    x2, y2 = pos_tree[line[1]]
                    line_long = x2 - x1
                    line_lat = y2 - y1
                    angle.append(atan(line_lat/line_long))

                mean_angle = np.mean(angle)
                angle = np.array(angle)
                # skinny inclination range pi/9. around mean_angle
                num_fil = np.where((angle > (mean_angle-pi/9.)) & (angle < (mean_angle+pi/9.)))[0].size

                if float(num_fil)/float(len(angle)) > 0.7:  # not parallel/vertical but skinny
                    likelihood = 0.5
                    # self.bg_fig.show_lines(box, color="blue", linewidth=2, alpha=0.8, zorder=10)
                else:
                    likelihood = 0.  # not parallel or skinny

            # (size, center_point_xpos, center_point_ypos, length_ratio, likelihood)
            tree = [each_tree.size(), delta_x/2., delta_y/2., length_ratio, likelihood]
            for col, col_value in enumerate(tree):
                ws.write(i, col, col_value)
            i += 1

        fig_name = './fil%d_output/Filament%d_%.2f_MST.png' % (self.fil_n, self.fil_n, self.threshold)
        plt.savefig(fig_name)

        wb_name = './fil%d_output/tree1_%.2f.xls' % (self.fil_n, self.threshold)
        workbook.save(wb_name)

    def analysis_filament_likelihood(self):
        pass
