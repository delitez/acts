add_library(
    ActsExamplesDetectorContextual
    SHARED
    src/AlignedDetector.cpp
    src/InternalAlignmentDecorator.cpp
    src/ExternalAlignmentDecorator.cpp
)

target_include_directories(
    ActsExamplesDetectorContextual
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetectorContextual
    PUBLIC
        ActsCore
        ActsExamplesFramework
        ActsExamplesDetectorsCommon
        ActsExamplesDetectorGeneric
)

acts_compile_headers(ActsExamplesDetectorContextual GLOB "include/**/*.hpp")

install(
    TARGETS ActsExamplesDetectorContextual
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
