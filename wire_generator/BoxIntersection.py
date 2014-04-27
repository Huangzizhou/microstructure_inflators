import numpy as np
from numpy.linalg import norm

from UniquePointExtractor import UniquePointExtractor

class BoxIntersection:
    def __init__(self, bbox_min, bbox_max):
        self.bbox_min = bbox_min;
        self.bbox_max = bbox_max;

    def intersect(self, vertices, faces):
        self.vertices = vertices;
        self.faces = faces;
        self.init_edges();

        self.intersections = [];
        for f in self.faces:
            self.__bbox_face_intersection(f);
        for e in self.edges:
            self.__bbox_edge_intersection(
                    self.vertices[e[0]], self.vertices[e[1]]);

        return self.__extract_unique_points(self.intersections, 1e-7);

    def init_edges(self):
        edges = set();
        for f in self.faces:
            edges.add(tuple(sorted([f[0], f[1]])));
            edges.add(tuple(sorted([f[1], f[2]])));
            edges.add(tuple(sorted([f[2], f[0]])));
        self.edges = np.array(list(edges));

    def __bbox_edge_intersection(self, p1, p2):
        eps = 1e-6;
        bbox = np.array([self.bbox_min, self.bbox_max]);
        direction = p2 - p1;
        effective_coord = direction != 0.0;
        step = np.divide((bbox-p1)[:,effective_coord],
                direction[np.newaxis, effective_coord]).ravel();
        candidates = np.outer(step, direction) + p1[np.newaxis, :];
        inside = np.logical_and(
                np.all(candidates < self.bbox_max+eps, axis=1), 
                np.all(candidates > self.bbox_min-eps, axis=1) );
        inside = np.logical_and(inside, step >= -eps);
        inside = np.logical_and(inside, step <= 1.0+eps);
        if np.any(inside):
            candidates = candidates[inside];
            for p in candidates:
                self.intersections.append(p);

    def __bbox_face_intersection(self, face):
        bbox_min = self.bbox_min;
        bbox_max = self.bbox_max;
        corners = np.array([
            [bbox_min[0], bbox_min[1], bbox_min[2]],
            [bbox_max[0], bbox_min[1], bbox_min[2]],
            [bbox_max[0], bbox_max[1], bbox_min[2]],
            [bbox_min[0], bbox_max[1], bbox_min[2]],
            [bbox_min[0], bbox_min[1], bbox_max[2]],
            [bbox_max[0], bbox_min[1], bbox_max[2]],
            [bbox_max[0], bbox_max[1], bbox_max[2]],
            [bbox_min[0], bbox_max[1], bbox_max[2]] ]);
        intersections = [];

        self.__edge_face_intersection(face, corners[0], corners[1]);
        self.__edge_face_intersection(face, corners[1], corners[2]);
        self.__edge_face_intersection(face, corners[2], corners[3]);
        self.__edge_face_intersection(face, corners[3], corners[0]);

        self.__edge_face_intersection(face, corners[4], corners[5]);
        self.__edge_face_intersection(face, corners[5], corners[6]);
        self.__edge_face_intersection(face, corners[6], corners[7]);
        self.__edge_face_intersection(face, corners[7], corners[4]);

        self.__edge_face_intersection(face, corners[0], corners[4]);
        self.__edge_face_intersection(face, corners[1], corners[5]);
        self.__edge_face_intersection(face, corners[2], corners[6]);
        self.__edge_face_intersection(face, corners[3], corners[7]);

    def __edge_face_intersection(self, face, p1, p2):
        eps = 1e-6;
        t0 = self.vertices[face[0]];
        t1 = self.vertices[face[1]];
        t2 = self.vertices[face[2]];
        normal = np.cross(t1 - t0, t2 - t0);
        normal = normal / norm(normal);

        def is_inside(p):
            e0 = t0 - p;
            e1 = t1 - p;
            e2 = t2 - p;
            a0 = np.dot(np.cross(e0, e1), normal);
            a1 = np.dot(np.cross(e1, e2), normal);
            a2 = np.dot(np.cross(e2, e0), normal);
            area = np.array([a0, a1, a2]);
            return np.all(area > -eps);


        d1 = np.dot(p1 - t0, normal);
        proj_len = np.dot(p1 - p2, normal);
        if abs(proj_len) < eps:
            p_normal_dist = np.dot(p1, normal);
            t_normal_dist = np.dot(t1, normal);
            if abs(p_normal_dist - t_normal_dist) < eps:
                if is_inside(p1):
                    self.intersections.append(p1);
                if is_inside(p2):
                    self.intersections.append(p2);
        else:
            frac = d1 / proj_len;
            if frac >= -eps and frac < 1.0 + eps:
                intersection = p1 * (1.0 - frac) + p2 * frac;
                if is_inside(intersection):
                    self.intersections.append(intersection);

    def __extract_unique_points(self, pts, eps):
        return UniquePointExtractor.extract(pts, eps);

