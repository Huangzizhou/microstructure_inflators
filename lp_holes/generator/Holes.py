import numpy as np
from numpy.linalg import norm
from math import pi, sqrt

class Holes(object):
    class Hole:
        def __init__(self, p, radius, num_samples=72):
            if p == "inf":
                p = np.inf;
            self.measure = p;
            self.radius = radius;
            self.num_samples = num_samples;
            self.boundary = self.__generate();
            self.edges = np.array([[i, (i+1)%num_samples] for i in
                range(num_samples)], dtype=int);

        def __generate(self):
            angles = np.arange(-pi, pi, 2*pi/self.num_samples)[:self.num_samples];
            x = np.cos(angles);
            y = np.sin(angles);
            lp_norm = norm([x,y], self.measure, axis=0).T;
            x = np.divide(x, lp_norm);
            y = np.divide(y, lp_norm);
            return np.array([x,y]).T * self.radius;

        def __compute_hole_area(self):
            area = 0.0;
            for i in range(self.num_samples):
                vi = self.boundary[i,:];
                vj = self.boundary[(i+1)%self.num_samples, :];
                area += norm(np.cross(vi, vj)) / 2.0;
            return area;

        def __normalize_hole_area(self):
            target_area = pi*radius*radius; # Use circular hole area as target
            hole_area = self.__compute_hole_area();
            self.boundary *= sqrt(target_area/hole_area);


    def __init__(self, width, height, radius, p):
        self.width = width;
        self.height = height;
        self.hole_profile = self.Hole(p, radius);
        self.centers = [];

    def drill_holes(self, centers):
        self.centers = centers;

    @property
    def vertices(self):
        num_sample_per_hole = self.hole_profile.num_samples;
        num_holes = len(self.centers);
        vtx = np.zeros((num_holes * num_sample_per_hole, 2));
        for i,center in enumerate(self.centers):
            hole_boundary = np.copy(self.hole_profile.boundary);
            hole_boundary = np.add(hole_boundary, center);
            vtx[i*num_sample_per_hole:(i+1)*num_sample_per_hole, :] =\
                    hole_boundary;
        return vtx;

    @property
    def edges(self):
        num_sample_per_hole = self.hole_profile.num_samples;
        num_holes = len(self.centers);
        edges = np.zeros((num_holes * num_sample_per_hole, 2), dtype=int);
        for i,center in enumerate(self.centers):
            hole_edges = np.copy(self.hole_profile.edges);
            hole_edges += i*num_sample_per_hole;
            edges[i*num_sample_per_hole:(i+1)*num_sample_per_hole, :] =\
                    hole_edges;
        return edges;

