# Initial creator (and more general one)
add_executable(HexaPillarsCreator HexaPillarsCreator.cpp)
set_target_properties(HexaPillarsCreator PROPERTIES CXX_STANDARD 11)

# Eigen3 library
find_package(Eigen3 REQUIRED)
target_link_libraries(HexaPillarsCreator PUBLIC Eigen3::Eigen ${Boost_LIBRARIES})

# Json
target_link_libraries(HexaPillarsCreator PUBLIC cli11 json)

# Boost library
find_package(Boost 1.55 REQUIRED COMPONENTS filesystem system program_options thread)
target_include_directories(HexaPillarsCreator SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(HexaPillarsCreator PUBLIC ${Boost_LIBRARIES})

# creator generating simpler wire structure
add_executable(HexaPillarsCreatorSimpler HexaPillarsCreatorSimpler.cpp)
set_target_properties(HexaPillarsCreatorSimpler PROPERTIES CXX_STANDARD 11)

# Eigen3 library
find_package(Eigen3 REQUIRED)
target_link_libraries(HexaPillarsCreatorSimpler PUBLIC Eigen3::Eigen ${Boost_LIBRARIES})

# Json
target_link_libraries(HexaPillarsCreatorSimpler PUBLIC cli11 json)

# Boost library
find_package(Boost 1.55 REQUIRED COMPONENTS filesystem system program_options thread)
target_include_directories(HexaPillarsCreatorSimpler SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(HexaPillarsCreatorSimpler PUBLIC ${Boost_LIBRARIES})