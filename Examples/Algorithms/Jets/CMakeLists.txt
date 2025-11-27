acts_add_library(ExamplesJets SHARED src/TruthJetAlgorithm.cpp)
target_include_directories(
    ActsExamplesJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesJets
    PUBLIC Acts::Fatras Acts::PluginFastJet Acts::ExamplesFramework Acts::ExamplesIoHepMC3 HepMC3::HepMC3 Acts::ExamplesIoEDM4hep EDM4HEP::edm4hep
)

acts_compile_headers(ExamplesJets GLOB "include/**/*.hpp")
