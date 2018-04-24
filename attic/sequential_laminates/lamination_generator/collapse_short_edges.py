#!/usr/bin/env python

import argparse
from heapq import heapify, heappop
import logging
import numpy as np
from numpy.linalg import norm

from mesh_io import load_mesh, form_mesh, save_mesh_raw
from PyMeshUtils import ShortEdgeRemoval

class EdgeCollapser:
    def __init__(self, vertices, faces):
        assert(faces.shape[1] == 3);
        self.vertices = vertices;
        self.faces = faces;
        self.logger = logging.getLogger(" {:>15} ".format(
            self.__class__.__name__[0:15]));

    def collapse(self, abs_threshold, rel_threshold):
        """ Note this method remove all edges with length less than threshold.
        This could result in a non-manifold mesh.
        """
        min_edge_length = abs_threshold;
        if rel_threshold is not None:
            ave_edge_len = self.__get_ave_edge_length()
            min_edge_length = rel_threshold  * ave_edge_len;
        self.logger.info("Minimum edge threshold: {:.3}".format(min_edge_length));

        num_collapsed = self.__collapse_C(min_edge_length);
        self.logger.info("{} edges collapsed".format(num_collapsed));

    def __compute_edge_lengths(self):
        mesh = form_mesh(self.vertices, self.faces);
        mesh.add_attribute("edge_length");
        edge_lengths = np.copy(mesh.get_face_attribute("edge_length"));

        v1 = self.faces[:,0].reshape((-1, 1));
        v2 = self.faces[:,1].reshape((-1, 1));
        v3 = self.faces[:,2].reshape((-1, 1));

        edge1 = np.hstack([np.minimum(v1, v2), np.maximum(v1, v2)]);
        edge2 = np.hstack([np.minimum(v2, v3), np.maximum(v2, v3)]);
        edge3 = np.hstack([np.minimum(v3, v1), np.maximum(v3, v1)]);

        self.edge_lengths = dict(
                [(tuple(edge), l) for edge,l in zip(edge1, edge_lengths[:,0])] +
                [(tuple(edge), l) for edge,l in zip(edge2, edge_lengths[:,1])] +
                [(tuple(edge), l) for edge,l in zip(edge3, edge_lengths[:,2])]);

    def __get_ave_edge_length(self):
        if not hasattr(self, "edge_lengths"):
            self.__compute_edge_lengths();
        return np.mean(self.edge_lengths.values());

    def __collapse_C(self, min_edge_length):
        collapser = ShortEdgeRemoval(self.vertices, self.faces);
        num_collapsed = collapser.run(min_edge_length);
        self.vertices = collapser.get_vertices();
        self.faces = collapser.get_faces();
        return num_collapsed;



def collapse_short_edges(vertices, faces, abs_threshold=0.0,
        rel_threshold=None):
    collapser = EdgeCollapser(vertices, faces);
    collapser.collapse(abs_threshold, rel_threshold);
    return collapser.vertices, collapser.faces;

def parse_args():
    parser = argparse.ArgumentParser(description="Collapse short edges");
    parser.add_argument("in_mesh", help="input mesh file");
    parser.add_argument("out_mesh", help="output mesh file");
    parser.add_argument("--rel-threshold",\
            help="relative minimal edge length comparing to ave",
            default=None, type=float);
    parser.add_argument("--abs-threshold",\
            help="absolute minimal edge length",
            default=0.1, type=float);
    parser.add_argument("--verbose",\
            help="print more info out", action="store_true");
    parser.add_argument("--timing",\
            help="printing timing info", action="store_true");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    if args.verbose:
        logging.basicConfig(level=logging.INFO);
    mesh = load_mesh(args.in_mesh);

    vertices = mesh.vertices;
    faces = mesh.faces;
    vertices, faces = collapse_short_edges(vertices, faces,
            abs_threshold = args.abs_threshold,
            rel_threshold = args.rel_threshold);
    save_mesh_raw(args.out_mesh, vertices, faces);

if __name__ == "__main__":
    main();
