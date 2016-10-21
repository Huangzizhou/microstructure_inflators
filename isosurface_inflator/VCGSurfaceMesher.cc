#include "VCGSurfaceMesherImpl.hh"

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>>>;
// Enable for slower builds...
// template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>>>;
