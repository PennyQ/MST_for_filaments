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
        start = time.time()
        self.create_graph()
        end = time.time()
        print('create_graph time is', end-start) # this is ok

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

        '''g_lon = tbdata['_Glon']
        g_lat = tbdata['_Glat']
        tau_av = tbdata['tau_av']

        # TODO: we could try to reduce the loop here later :)
        for i in range(len(g_lon)):
            if self.max_l >= g_lon[i] >= self.min_l and 1 >= g_lat[i] >= -1:
                self.lat.append(g_lat[i])
                self.lng.append(g_lon[i])
                self.graph.add_node(i, posxy=(self.lng[-1], self.lat[-1]))'''

        positions = nx.get_node_attributes(self.graph, 'posxy')

        # build edges of graph, within a certain distance treshold
        for n, d in self.graph.nodes_iter(data=True):
            xn, yn = positions[n]
            for m, g in self.graph.nodes_iter(data=True):
                xm, ym = positions[m]
                dist = (math.sqrt(math.pow((xm - xn), 2) + math.pow((ym - yn), 2)))
                if dist <= self.threshold:   # dist is the threshold here
                    # tau_av is larger then the opacity is larger, so the weight is lower
                    # IRDC yes then weight is lower -> +0 no ->+1
                    is_irdc = 1.
                    if covered_data['IRDC'][n] is 'y':
                        is_irdc -= .5
                    if covered_data['IRDC'][m] is 'y':
                        is_irdc -= .5
                    # TODO: choose a more meaningful way of defining weight here
                    weight = dist*20. - (covered_data['tau_p'][n] + covered_data['tau_p'][m])/2. + is_irdc
                    if weight < 0:
                        weight = 0.
                    print('new weight and old', weight, dist, is_irdc)
                    self.graph.add_edge(n, m, weight=weight)

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

        start = time.time()
        tree_list = self.get_tree_list()
        end = time.time()
        print('get_tree_list', end-start)  # 0.014

        # draw MST on top of bg figure
        start = time.time()
        fig = self.get_bg_figure()
        end = time.time()
        print('get_bg_figure', end-start)

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
            each_tree_weight = 0
            # draw mst trees here
            for each in sorted(each_tree.edges(data='weight')):
                # print('what is each edge', each, each[2])
                # ('what is each edge', (323, 337, 2.2539938916387352), 2.2539938916387352)
                each_tree_weight += each[2]

                x1, y1 = pos_tree[each[0]]
                x2, y2 = pos_tree[each[1]]
                mst_long = [x1, x2]
                mst_lat = [y1, y2]
                edge = [np.vstack((mst_long, mst_lat))]

                # TODO: can I use a new fig as the bg_fig, and draw upon it?
                if each[2] > 2:
                    # TODO: show a box around this edge and its neighbour edge
                    fig.show_lines(edge, color='blue', alpha=0.6, linewidth=4)
                    print('did I draw here?')
                elif each[2] > 1:
                    fig.show_lines(edge, color='aquamarine', alpha=0.6, linewidth=4)
                else:
                    # these are what we want
                    fig.show_lines(edge, color='yellow', alpha=0.6, linewidth=4)

            '''# draw box, avoid wide/high box
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
                fig.show_lines(box, color="yellow", linewidth=2, alpha=0.8, zorder=10)

            #  elif (11./18.)*pi > fabs(length_ratio) > (9./18.)*pi:  # vertical and skinny
            elif fabs(length_ratio) > (7./18.)*pi:  # vertical and skinny
                fig.show_lines(box, color="blue", linewidth=2, alpha=0.8, zorder=10)
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
                    fig.show_lines(box, color="blue", linewidth=2, alpha=0.8, zorder=10)
                else:
                    likelihood = 0.  # not parallel or skinny'''

            # (size, center_point_xpos, center_point_ypos, length_ratio, likelihood)
            # tree = [each_tree.size(), delta_x/2., delta_y/2., length_ratio, likelihood]
            ave_weight = float(each_tree_weight)/float(each_tree.size())
            tree = [each_tree.size(), delta_x/2., delta_y/2., each_tree_weight, ave_weight]

            for col, col_value in enumerate(tree):
                ws.write(i, col, col_value)
            i += 1

        fig_name = './fil%d_output/Filament%d_%.2f_MST.png' % (self.fil_n, self.fil_n, self.threshold)
        plt.savefig(fig_name)

        wb_name = './fil%d_output/tree1_%.2f.xls' % (self.fil_n, self.threshold)
        workbook.save(wb_name)
