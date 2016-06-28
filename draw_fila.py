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
from xlwt import *
import time
import os
from math import *


class DrawFilament():
    # path to file
    CLOUD_LB = '/Users/penny/Works/MST_filaments/Peretto_Fuller_data/cloudfitsLB.txt'
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
        start = time.time()
        self.create_graph()
        end = time.time()
        print('create_graph time is', end-start)

        self.fil_n = fila_n

        self.fil_candidate_file = '/Users/penny/Works/MST_filaments/' \
                                  'Extinction_filaments_data/Candid%d_cropped_final.fits' % self.fil_n
        if not os.path.exists('./fil%d_output' % self.fil_n):
            os.makedirs('./fil%d_output' % self.fil_n)

    def create_graph(self):
        # create arrays for storing coordinates of molecular clouds
        cloudfitsLB = open(self.CLOUD_LB)
        for n, line in enumerate(cloudfitsLB):
            columns = line.split(", ")

            if ((float(columns[0]) >= self.min_l) and (float(columns[0]) <= self.max_l) and
                        (float(columns[1])) <= 1 and (float(columns[1]) >= (-1))):
                self.lat.append(float(columns[0]))
                self.lng.append(float(columns[1]))
                self.graph.add_node(str(n), posxy=(self.lat[-1], self.lng[-1]))

        positions = nx.get_node_attributes(self.graph, 'posxy')

        # build edges of graph, within a certain distance treshold
        for n, d in self.graph.nodes_iter(data=True):
            xn, yn = positions[n]
            for m, g in self.graph.nodes_iter(data=True):
                xm, ym = positions[m]
                dist = (math.sqrt(math.pow((xm - xn),2) + math.pow((ym - yn), 2)))
                if dist <= self.threshold:   # dist is the threshold here
                    self.graph.add_edge(n, m, weight=dist)

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

        # TODO: Q~for draw three parallel lines, why it's called rectangles??
        fig.show_rectangles(arml, armb, 0.1, 0.003, color="magenta", cmap=cmap,  #0.03, 0.002 for fila1
                            norm=norm, edgecolors="none", linewidth=2)
        fig.show_rectangles(arml, upperb, 0.3, 0.001, color="red", cmap=cmap,
                            norm=norm, edgecolors="none", linewidth=2)
        fig.show_rectangles(arml, lowerb, 0.3, 0.001, color="red", cmap=cmap,
                            norm=norm, edgecolors="none", linewidth=2)
        fig.show_markers(self.lat, self.lng, s=30, alpha=1, zorder=2,
                         marker='^', c='yellow')
        fig.show_lines(box, color="crimson", linewidth=3, zorder=10)

        return fig

    # def draw_tree_dandidates(self, each_tree):

    def get_tree_list(self):
        # using Kruskal's algorithm
        T = nx.minimum_spanning_tree(self.graph)

        # using Prim's algorithm
        # T = nx.prim_mst(self.graph)

        # start = time.time()
        # fig = self.draw_figure()
        # end = time.time()
        # print('draw_figure time is', end-start)  # 11.197342157363892

        # Tree Diagnosis
        # number of trees in the graph
        # ntrees = nx.number_connected_components(T)

        #generate lists of individual trees
        tree_list = nx.connected_component_subgraphs(T)

        # #largest tree in the region
        # max_tree = max(tree_list, key = len)
        # print "largest tree in the region: " + str(max_tree)
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
        print('get_tree_list+draw figure time is ', end-start)  # 24.365610122680664

        # print edges of the MST and overlay them on top of GLIMPSE file
        # pos_mst = nx.get_node_attributes(tree_list, 'posxy')
        # draw MST upon the bone figure
        fig = self.get_bg_figure()

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


            #  draw here
            for each in sorted(each_tree.edges(data=True)):
                x1, y1 = pos_tree[each[0]]
                x2, y2 = pos_tree[each[1]]
                mst_long = [x1, x2]
                mst_lat = [y1, y2]
                edge = [np.vstack((mst_long, mst_lat))]
                fig.show_lines(edge, color='aquamarine', alpha=0.6, linewidth=4)

            # skip tiny pieces
            # TODO: the 0.1 threshold here should be replaced with others later
            if self.get_distance(delta_x, delta_y) < 0.1:
                continue

            length_ratio = atan2(delta_y, delta_x)  #return between -pi and pi

            # if length_ratioo ~ 0 -> parallel   -> PI/9 as a threshold first
            #                  ~ PI/2 -> vertical
            if fabs(length_ratio) < pi/9 or fabs(length_ratio) > (8./9.)*pi:
                # save to file  and draw
                likelihood = 1

                # box example [26.8, 27.1, 27.1, 26.8, 26.8], [-0.23, -0.23, -0.4, -0.4, -0.23]
                # draw box
                fil_x_box = [x_min, x_max, x_max, x_min, x_min]
                fil_y_box = [y_max, y_max, y_min, y_min, y_max]
                box = [np.vstack((fil_x_box, fil_y_box))]
                fig.show_lines(box, color="yellow", linewidth=2, alpha=0.6, zorder=10)

            # elif (11./18.)*pi > fabs(length_ratio) > (9./18.)*pi: # vertical
            #     likelihood = 0.5

            else:
                likelihood = 0.

            # skinny 

            print(nx.info(each_tree))
            # (size, center_point_xpos, center_point_ypos, length_ratio, likelihood)
            tree = [each_tree.size(), delta_x/2., delta_y/2., length_ratio, likelihood]
            for col, col_value in enumerate(tree):
                ws.write(i, col, col_value)
            i += 1



            # number of nodes
            # nnodes = each.number_of_nodes()

            # size similar to above
            # tree_size = each.size()

            # tree total length
            # tree_length = each.size(weight='weight')

            # tree diameter :calculated in number of vertices which must be visited
            # in order to travel from one node to another
            # diam = nx.diameter(each)



            # average degree
            # avg_deg = float(nnodes) / tree_size

            # http://networkx.readthedocs.io/en/stable/reference/algorithms.clustering.html?highlight=clustering ??
            # what is clustering coefficient: a measure of the degree to which nodes tend to cluster
            # cc = nx.clustering(each)
            # avg_cc = sum(cc.values()) / len(cc)

            # tree density??
            # for undirected graphs, tree density = 2m/ (n*(n-1)), m = # of edges, n = # of nodes
            # rho = nx.density(each)



            # average inclination angle
            '''angle = []
            for line in sorted(each.edges(data=True)):  # sorted start from 1st dim
                x1, y1 = pos_tree[line[0]]
                x2, y2 = pos_tree[line[1]]
                line_long = x2 - x1
                line_lat = y2 - y1
                angle.append(math.degrees(math.atan(line_lat/line_long)))
            avg_angle = np.mean(angle)'''
        fig_name = './fil%d_output/Filament%d_%.2f_MST.png' % (self.fil_n, self.fil_n, self.threshold)
        plt.savefig(fig_name)

        wb_name = './fil%d_output/tree1_%.2f.xls' % (self.fil_n, self.threshold)
        workbook.save(wb_name)
