import networkx as nx
import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib as mpl
import aplpy
from ds9norm import DS9Normalize
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import math
import decimal
from astropy.table import Table
import xlwt
from xlwt import *

# path to file containing galactic coordinates of Peretto & Fuller
flb = open('/Users/penny/Works/MST_filaments/Peretto_Fuller_data/cloudfitsLB.txt')

#create arrays for storing coordinates of molecular clouds
lat = []
lng = []
N = 0

#create networkx Graph, to be populated with points from Peretto & Fuller catalog
G = nx.Graph()

for line in flb:
	columns = line.split(", ")

	if ((float(columns[0]) >= 25.5) and (float(columns[0]) <= 28.5) and (float(columns[1])) <= 1 and (float(columns[1]) >= (-1))):
		lat.append(float(columns[0]))
		lng.append(float(columns[1]))
		G.add_node(str(N), posxy = (lat[-1], lng[-1]))
		N += 1
positions = nx.get_node_attributes(G,'posxy')

# build edges of graph, within a certain distance treshold
for n,d in G.nodes_iter(data=True):
	xn, yn = positions[n]
	for m,g in G.nodes_iter(data=True):
		xm, ym = positions[m]
		dist = (math.sqrt(math.pow((xm - xn),2) + math.pow((ym - yn),2)))
		if (dist <= 0.05):   #dist is the threshold here
			G.add_edge(n, m, weight = dist)


#===============draw figures start========================
myfig = plt.figure(figsize=(20,20))
fig=aplpy.FITSFigure("/Users/penny/Works/MST_filaments/Extinction_filaments_data/Candid1_cropped_final.fits",figure=myfig, zorder = 0)
norm = DS9Normalize(stretch='sqrt', clip_hi=99,clip_lo=1, contrast=1.56, bias=0.65)
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
plt.title("Filament1 Kruskal MST, threshold 0.1deg")
fig.axis_labels.set_font(size=20)
fig.tick_labels.set_font(size=20)

# the galactic longitude range of the output figure
scutumdata=Table.read("/Users/penny/Works/MST_filaments/Extinction_filaments_data/Scutum_final.txt", format='ascii', delimiter='\t', guess=False)
ws = (scutumdata['long_final'] > 25.5) & (scutumdata['long_final'] < 28.5)
arml=scutumdata['long_final'][ws]
num_rows=arml.size

cmap=plt.cm.get_cmap("gist_rainbow")
cmap.set_under(color="k")
cmap.set_over(color="k")
norm=mpl.colors.Normalize(clip=False,vmin=65,vmax=85)
m=plt.cm.ScalarMappable(cmap=cmap,norm=norm)

#Filament 1 coordinates
upperb=[0.06]*num_rows
armb=[-0.19]*num_rows
lowerb=[-0.44]*num_rows
armv = scutumdata['Vlsr'][ws]

xbox = [26.8, 27.1, 27.1, 26.8, 26.8]
ybox = [-0.23, -0.23, -0.4, -0.4, -0.23]
box = [np.vstack((xbox, ybox))]

fig.show_rectangles(arml, armb, 0.03, 0.002, color = "magenta",cmap = cmap, norm=norm, edgecolors="none", linewidth=2)
fig.show_rectangles(arml, upperb, 0.3, 0.001, color="red",cmap = cmap, norm=norm, edgecolors="none", linewidth = 2)
fig.show_rectangles(arml, lowerb, 0.3, 0.001, color="red",cmap = cmap, norm=norm, edgecolors="none", linewidth = 2)
fig.show_markers(lat, lng, s = 30, alpha = 1, zorder = 2, marker = '^', c = 'yellow')
fig.show_lines(box, color = "crimson", linewidth = 3)

# using Kruskal's algorithm
T = nx.minimum_spanning_tree(G)

# using Prim's algorithm
# T = nx.prim_mst(G)

#print edges of the MST and overlay them on top of GLIMPSE file
pos_mst = nx.get_node_attributes(T,'posxy')
for each in sorted(T.edges(data=True)):
	x1, y1 = pos_mst[each[0]]
	x2, y2 = pos_mst[each[1]]
	mst_long = [x1, x2]
	mst_lat = [y1, y2]
	edge = [np.vstack((mst_long, mst_lat))]
	fig.show_lines(edge, color = 'aquamarine', linewidth = 2.5)

# Tree Diagnosis
# number of trees in the graph
ntrees = nx.number_connected_components(T)
print "number of trees = " + str(ntrees)

#generate lists of individual trees
tree_list = nx.connected_component_subgraphs(T)

# #largest tree in the region
# max_tree = max(tree_list, key = len)
# print "largest tree in the region: " + str(max_tree)

i = 1
w = Workbook()
ws = w.add_sheet('Filament1')
for each in tree_list:
	#number of nodes
	nnodes = each.number_of_nodes()

	#size
	tree_size = each.size()

	#tree total length
	tree_length = each.size(weight = 'weight')

	#tree diameter
	#calculated in number of vertices which must be visited in order to travel from one node to another
	diam = nx.diameter(each)

	#average degree
	avg_deg = float(nnodes) / tree_size

	#clustering coefficient
	cc = nx.clustering(each)
	avg_cc = sum(cc.values()) / len(cc)

	#tree density
	# for undirected graphs, tree density = 2m/ (n*(n-1)), m = # of edges, n = # of nodes
	rho = nx.density(each)

	#node centralities
	#measuring length/ width ratio of tree
	pos_tree = nx.get_node_attributes(each,'posxy')
	x_val = []
	y_val = []
	for n,d in each.nodes_iter(data=True):
		xn, yn = pos_tree[n]
		x_val.append(xn)
		y_val.append(yn)
	x_max = max(x_val)
	x_min = min(x_val)
	y_max = max(y_val)
	y_min = min(y_val)
	if (x_min <= 26.94 <= x_max) and (y_min <= -0.3 <= y_max):
		bone_likelihood = 1
	else:
		bone_likelihood = 0
	delta_x = x_max - x_min
	delta_y = y_max - y_min
	length_ratio = delta_x / delta_y

	#average inclination angle
	angle = []
	for line in sorted(each.edges(data=True)):
		x1, y1 = pos_tree[line[0]]
		x2, y2 = pos_tree[line[1]]
		line_long = x2 - x1
		line_lat = y2 - y1
		angle.append(math.degrees(math.atan(line_lat/ line_long)))
	avg_angle = np.mean(angle)

	tree = [nnodes, tree_size, tree_length, diam, avg_deg, avg_cc, rho, delta_x, delta_y, length_ratio, avg_angle, bone_likelihood]
	for col, col_value in enumerate(tree):
		ws.write(i, col, col_value)
	i += 1

w.save('tree1.xls')
plt.savefig("Filament1 Kruskal 0.1 MST.png")
plt.show()