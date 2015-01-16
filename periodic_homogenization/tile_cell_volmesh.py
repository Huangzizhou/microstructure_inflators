#!/usr/bin/env python

import sys
import argparse
import os
import os.path
import numpy as np

PYMESH_PATH = os.environ.get("PYMESH_PATH");
if PYMESH_PATH is None:
    raise ImportError("Please set PYMESH_PATH to the correct lib path.")
sys.path.append(os.path.join(PYMESH_PATH, "lib"));
sys.path.append(os.path.join(PYMESH_PATH, "swig"));
LINEAR_ELASTICITY_PATH = os.environ.get("LINEAR_ELASTICITY_PATH");
sys.path.append(LINEAR_ELASTICITY_PATH);
sys.path.append(os.path.join(LINEAR_ELASTICITY_PATH, "PyUtils"));

from mesh_io import *
from Mesh import Mesh
import PyMesh

def parse_args():
    parser = argparse.ArgumentParser(description="Base Cell Tiler");
    parser.add_argument("in_mesh", help="input mesh file (.msh)");
    parser.add_argument("out_mesh", help="tiled surface mesh output file (.obj, .msh)");
    parser.add_argument("num_tile_x", help="number of tiles in x direction");
    parser.add_argument("num_tile_y", help="number of tiles in y direction");
    parser.add_argument("num_tile_z", help="number of tiles in z direction");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    mesh = load_mesh(args.in_mesh);

    tets = mesh.voxels    # N x 4 matrix of integer indices
    vtxs = mesh.vertices  # N x 3 matrix of floats

    # Get pairs of vertices on opposite periodic cell boundary faces
    import PeriodicBoundaryCondition as pbc
    bmin,bmax = mesh.bbox;
    vtxs = mesh.vertices;
    # List of each dimension's indentified vertex pairs across the periodic
    # boundaries
    identifiedFacePairs = []
    for d in range(mesh.dim):
        minFaceVertices, = np.where(np.abs(vtxs[:,d] - bmin[d]) < 1e-5);
        maxFaceVertices, = np.where(np.abs(vtxs[:,d] - bmax[d]) < 1e-5);
        # print len(minFaceVertices), " in min face ", d, ": "
        # print mesh.vertices[minFaceVertices]
        # print "max face ", d, ": "
        # print len(maxFaceVertices), " in max face ", d, ": "
        # print mesh.vertices[maxFaceVertices]
        assert(minFaceVertices.shape == maxFaceVertices.shape);
        pairs = pbc.match_boundaries( mesh, minFaceVertices, maxFaceVertices);
        identifiedFacePairs.append(pairs);

    numCellVertices = vtxs.shape[0];
    numCellTets = tets.shape[0];

    fullTets = np.zeros((0, 4), dtype=int);
    fullVerts = np.zeros((0, 3));
    nx, ny, nz = map(int, [args.num_tile_x, args.num_tile_y, args.num_tile_z]);
    vtxOffset = lambda i, j, k: (i * nz * ny + j * nz + k) * numCellVertices

    # Glue in (i, j, k) scanline order, attaching a copy of the base cell to the
    # existing structure by merging duplicated vertices. These duplicated
    # vertices for neighbors along dimension d are specified by the entries in
    # identifiedFacePairs[d]. In this scanline order, only neighbors
    # (i - 1, j, k), (i, j - 1, k), (i, j, k - 1) >= (0, 0, 0) exist.
    #      *-------*
    #     /   /   /|
    #    /---+---/ |
    #   /   /   /|/|     ^ j
    #  *---+---* + *     |
    #  |   |   |/|/      |     i
    #  |---+---| /       *----->
    #  |   |   |/       /
    #  *-------*       v k
    # Gluing is done by overwriting tets' "tiled vertex indices" (i.e. index
    # before merging duplicates) with the corresponding "glued vertex indices"
    # (index after merging duplicates). For vertices that end up merging, the
    # final glued index is stored in a dictionary mapping the vertices'
    # original tiled vertex indices to the index they are merged to.
    #
    # Note: This scanline approach only works on geometry that is actually
    # triply periodic. For instance, it failed on James' zigzag shape that is
    # missing geometry in one of its corners.

    # Update a global tiled vertex index -> glued vertex index map.
    # Tiled vertex indieces are of the form:
    #   base_vertex_index + vtxOffset(i, j, k)
    # We verify with assets that edge/corner cases (where multiple merges
    # happen) are resolved properly.
    gluedVertexIndex = {};
    def updateGluedVertexIndices(i, j, k):
        myOffset = vtxOffset(i, j, k);
        newGlue = {}
        # Vertices on the "min face" in each dimension are linked to their
        # paired vertex on the corresponding neighbor's "max face"
        if (i > 0):
            for vmin, vmax in identifiedFacePairs[0]:
                newGlue[vmin + myOffset] = \
                    gluedVertexIndex.get(vmax + vtxOffset(i - 1, j, k),
                            vmax + vtxOffset(i - 1, j, k));
        if (j > 0):
            for vmin, vmax in identifiedFacePairs[1]:
                rep = gluedVertexIndex.get(vmax + vtxOffset(i, j - 1, k),
                            vmax + vtxOffset(i, j - 1, k));
                if (vmin + myOffset in newGlue):
                    if (newGlue[vmin + myOffset] != rep):
                        print "ERROR: v", vmin + myOffset, "wants to link to ",\
                              rep, " and ", newGlue[vmin + myOffset]
                        assert(newGlue[vmin + myOffset] == rep);
                else:
                    newGlue[vmin + myOffset] = rep;
        if (k > 0):
            for vmin, vmax in identifiedFacePairs[2]:
                rep = gluedVertexIndex.get(vmax + vtxOffset(i, j, k - 1),
                            vmax + vtxOffset(i, j, k - 1));
                if (vmin + myOffset in newGlue):
                    if (newGlue[vmin + myOffset] != rep):
                        print "ERROR: v", vmin + myOffset, "wants to link to ",\
                              rep, " and ", newGlue[vmin + myOffset]
                        assert(newGlue[vmin + myOffset] == rep);
                else:
                    newGlue[vmin + myOffset] = rep;
        gluedVertexIndex.update(newGlue);

    # Glue everything together
    cellSize = bmax - bmin;
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                updateGluedVertexIndices(i, j, k);
                myOffset = vtxOffset(i, j, k);
                xyzOffset = cellSize * np.array([i, j, k]);
                fullVerts = np.vstack([fullVerts, vtxs + xyzOffset]);
                gluedTets = np.copy(tets);
                vIdxs = gluedTets.ravel();
                for vi in range(len(vIdxs)):
                    tiledVIdx = vIdxs[vi] + myOffset;
                    vIdxs[vi] = gluedVertexIndex.get(tiledVIdx, tiledVIdx);
                fullTets = np.vstack([fullTets, gluedTets]);
    # print "verts: ", fullVerts
    # print "tets: ", fullTets

    # print max(fullTets.ravel())

    gluedMesh = form_mesh(fullVerts, np.array([]), fullTets);
    save_mesh(args.out_mesh, gluedMesh);

if __name__ == "__main__":
    main();

