<!-- MarkdownTOC autolink="true" bracket="round" depth=3 -->
<!-- /MarkdownTOC -->


Dependencies included directly or as submodules:

- MeshFEM
- [Ceres](https://github.com/ceres-solver/ceres-solver) + [glog](https://github.com/google/glog), optional (TODO)
- [NLopt](https://github.com/stevengj/nlopt), optional (TODO)

Dependencies *not* included:

- CGAL, required

TODOS:

- [ ] Replace bfgs from Dlib by [CppNumericalSolvers](https://github.com/PatWie/CppNumericalSolvers)
- [ ] Restore the gitbook doc and clean it up + use a single doc for both MeshFEM and this repo

### Folders

| Folder | Description |
|--------|-------------|
| `Documentation/`        | Draft of documentation. |
| `isosurface_inflator/`  | Generate tet-mesh from a graph labeled with vertex thicknesses. |
| `pattern_optimization/` | Optimize the parameters of a pattern to attain target elastic properties. |
| `patterns/`             | List possible pattern topologies in 2D/3D. |
| `worst_case_stress/`    | Compute the worst-case stress for a given pattern. |
| `slice_supporter/`      | Generate raft + support structures below a microstructure (operates on slice images). |
| `topology_enumeration/` | Combinatorially generate cubic topologies on the base tetrahedron |
