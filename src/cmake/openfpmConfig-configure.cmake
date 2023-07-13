get_property(OPENFPM_INCLUDES TARGET openfpm::binary_config_ PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(OPENFPM_DEFINITION TARGET openfpm::binary_config_ PROPERTY INTERFACE_COMPILE_DEFINITIONS)
get_property(OPENFPM_LIBS TARGET openfpm::binary_config_ PROPERTY INTERFACE_LINK_LIBRARIES)
get_property(OPENFPM_COMPILE_OPTIONS TARGET openfpm::binary_config_ PROPERTY INTERFACE_COMPILE_OPTIONS)

if (openfpm_USE_SHARED_LIBS)
    target_link_libraries(openfpm::binary_config INTERFACE openfpm::binary_config_shared)
else()
    target_link_libraries(openfpm::binary_config INTERFACE openfpm::binary_config_static)
endif()

target_link_libraries(openfpm::binary_config INTERFACE openfpm::binary_config_)