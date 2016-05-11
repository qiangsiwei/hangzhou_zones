# -*- coding: utf-8 -*-

from math import exp
import networkx as nx
import numpy as np
import pylab as pl,random
from mpl_toolkits.mplot3d import Axes3D

def draw_figure(fignum,data,pos):
    fig = pl.figure(fignum)
    pl.clf()
    ax = fig.add_subplot(111)
    x,y = zip(*data)
    x1,y1 = zip(*pos)
    ax.scatter(x,y,c='b',alpha=0.1,s=5)
    ax.scatter(x1,y1,c="r",alpha=1,s=20)
    ax.set_title('Self Organizing Map')
    pl.savefig(str(fignum)+'.png')
    pl.close(fignum)
    fignum = fignum + 1

def read_coordinate(filename):
    import fileinput
    data = []
    for line in fileinput.input(filename):
        part = line.strip().split(" ")
        if int(part[2])!=0:
            for i in xrange(int(part[2])/int(1e2)+1):
                data.append([int(part[0]),int(part[1])])
    return data

class CSom():
    def __init__(self,data,num_of_nodes = 100):
        self.a = None
        self.data = data
        self.codebook_vectors = []
        self.graph = nx.Graph()
        self.count_of_neurons = 0
        self.pos = None
        self.tmax = None
        self.eps_e = None
        self.imagecount = 0
        temp_list = np.random.random_integers(0, len(data),num_of_nodes)
        temp_list = sorted(temp_list)
        for i in temp_list:
            self.graph.add_node(self.count_of_neurons, pos=(self.data[i]))
            self.count_of_neurons+=1

    def get_distance_matrix(self):
        vertices = self.graph.nodes()
        distance_matrix = []
        size_of_mat = len(vertices)
        for i in xrange(0,size_of_mat):
            distance_matrix.append([])
            for j in xrange(0,size_of_mat):
                path_len = len(nx.shortest_path(self.graph,i,j))-1
                distance_matrix[i].append(path_len)
        return distance_matrix

    def euclidean_distance(self, a,b):
        return ((a[0]-b[0])**2+(a[1]-b[1])**2)**0.5
        
    def find_closest_codebookvector(self,x):
        min_dist = float("inf")
        self.pos = nx.get_node_attributes(self.graph,'pos')
        for node,coordinates in self.pos.iteritems():
            dist = self.euclidean_distance(x,coordinates)
            if  dist < min_dist:
                min_dist = dist
                min_node = node
                coord_of_bmu = coordinates
        return min_node, coord_of_bmu

    def get_new_pos(self,sample_x,vector,Dij):
        mulfactor1 = [(sample_x[0]-vector[0]),(sample_x[1]-vector[1])]
        mulfactor2 = exp(-float(Dij))
        net_move = [mulfactor1[0]*mulfactor2*self.eps_e, mulfactor1[1]*mulfactor2*self.eps_e]
        new_position = [vector[0]+net_move[0], vector[1]+net_move[1]]
        return new_position

    def update_vectors(self,sample_x):
        winnernode, win_node_coord = self.find_closest_codebookvector(sample_x)
        self.pos = nx.get_node_attributes(self.graph,'pos')
        for node in self.graph.nodes():
            pos_of_node = self.pos[node]
            newpos = self.get_new_pos(sample_x,pos_of_node,self.euclidean_distance(win_node_coord,pos_of_node))
            self.graph.add_node(node, pos=newpos)

    def train(self,max_iterations=10000):
        self.tmax = max_iterations
        fignum = 0
        stepsize = 10
        for i in xrange(0,max_iterations+1):
            print "Iterating..{0:d}".format(i)
            div_factor = (float(i)/float(self.tmax))
            self.eps_e = (1-div_factor)
            random.shuffle(self.data)
            for x in self.data:
                self.update_vectors(x)
            check = i%stepsize
            if check == 0:
                fignum += 1
                self.draw_updated_fig(fignum)

    def draw_updated_fig(self,fignum):
        dictn= nx.get_node_attributes(self.graph,'pos')
        pos =[]
        for node,coordinates in dictn.iteritems():
            pos.append(coordinates)
        pos.append(dictn[0])
        draw_figure(fignum,data,pos)

if __name__ == '__main__':
    import fileinput

    data = read_coordinate("all_user_count.txt")
    obj = CSom(data, num_of_nodes=400)
    if obj is not None:
        obj.train(20)
    with open("output_si.txt",'w') as f:
        positions = nx.get_node_attributes(obj.graph,'pos')
        for i,j in positions.iteritems():
            f.write(str(j)+"\n")

    sommap = {}
    for line in fileinput.input("output_si.txt"):
        part = line.strip().replace("[","").replace("]","").split(",")
        px, py = float(part[0])*5,1000-float(part[1])*5
        sommap[tuple([px,py])] = 1
    fileinput.close()
    data = [list(item) for item in sommap.keys()]
    print len(data)

    import pytesselate
    import itertools
    import pydraw

    polygons = pytesselate.voronoi(data)
    crs = pydraw.CoordinateSystem([0,0,1000,1000])
    img = pydraw.Image(1000,1000)
    for i,(center,poly) in enumerate(polygons):
        if i%8 == 0:
            img.drawpolygon(poly,fillcolor=(0,0,0),outlinecolor=(0,0,0))
        if i%8 == 1:
            img.drawpolygon(poly,fillcolor=(255,0,0),outlinecolor=(0,0,0))
        if i%8 == 2:
            img.drawpolygon(poly,fillcolor=(0,255,0),outlinecolor=(0,0,0))
        if i%8 == 3:
            img.drawpolygon(poly,fillcolor=(0,0,255),outlinecolor=(0,0,0))
        if i%8 == 4:
            img.drawpolygon(poly,fillcolor=(0,255,255),outlinecolor=(0,0,0))
        if i%8 == 5:
            img.drawpolygon(poly,fillcolor=(255,0,255),outlinecolor=(0,0,0))
        if i%8 == 6:
            img.drawpolygon(poly,fillcolor=(255,255,0),outlinecolor=(0,0,0))
        if i%8 == 7:
            img.drawpolygon(poly,fillcolor=(255,255,255),outlinecolor=(0,0,0))
        if center:
            img.drawsquare(*center[:2], fillsize=2, fillcolor=(255,0,0))
    img.save("soms.png")
    img.view()

