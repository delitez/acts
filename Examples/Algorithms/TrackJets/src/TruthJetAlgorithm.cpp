// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackJets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <ostream>
#include <stdexcept>

ActsExamples::TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }

  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
  m_outputJets.initialize(m_cfg.outputJets);
  
}

ActsExamples::ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& truthParticles = m_inputTruthParticles(ctx);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputParticles;
  int particleIndex = -1;
  for (const auto& particle : truthParticles) {

    particleIndex++;

    fastjet::PseudoJet pseudoJet( particle.momentum().x(),
                                  particle.momentum().y(),
                                  particle.momentum().z(),
                                  particle.energy() );
    pseudoJet.set_user_index(particleIndex);
    inputParticles.push_back(pseudoJet);

    }

    // Run the jet clustering
    fastjet::ClusterSequence clusterSeq(inputParticles, DefaultJetDefinition);

    // Get the jets above pT 20 GeV
    std::vector<fastjet::PseudoJet> jets = clusterSeq.inclusive_jets(20 * Acts::UnitConstants::GeV);

    // Store the jets in the output data handle
    m_outputJets(ctx, std::move(jets));
//   }
  return ProcessCode::SUCCESS;
} 