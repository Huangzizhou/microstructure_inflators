################################################################################
# Source files
################################################################################

set(SOURCES
		BaseCellOperations.hh
		IterateFactory.hh
		IterateManager.hh
		LookupTable.hh
		ObjectiveTerm.hh
		ObjectiveTermNormalizations.hh
		PatternOptimizationJob.cc
		PatternOptimizationJob.hh
		PatternOptimizationIterate.hh
		PatternOptimizationConfig.hh
		SDConversions.hh
		SDConversionsOrthoCell.hh
		ShapeVelocityInterpolatorOrthoCell.hh
		objective_terms/IsotropicFit.hh
		objective_terms/IsotropicFitRel.hh
		objective_terms/ProximityRegularization.hh
		objective_terms/TargetVolume.hh
		objective_terms/TensorFit.hh
		constraints/Printability.hh
		constraints/TensorFit.hh
)

add_library(pattern_optimization ${SOURCES})
micro_copy_headers(pattern_optimization ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(pattern_optimization
	PUBLIC
		inflators
		optimizers
	PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
	add_sanitizers(pattern_optimization)
endif()
