if (AMReX_SPACEDIM EQUAL 3)
    add_executable(SphereTurbIC testSphereTurbIC.cpp ../../turbulence/TurbDataReader.cpp ${QuokkaObjSources})
endif()
