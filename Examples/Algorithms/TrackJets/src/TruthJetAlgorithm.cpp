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

#include <ostream>
#include <stdexcept>

ActsExamples::TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }

  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
}

ActsExamples::ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& truthParticles = m_inputTruthParticles(ctx);

  ACTS_INFO("event " << ctx.eventNumber << " collection '"
                     << m_cfg.inputTruthParticles << "' contains "
                     << truthParticles.size() << " truthParticles");
  for (const auto& particle : truthParticles) {
    ACTS_INFO("  particle " << particle);
    ACTS_INFO("    process_type: " << particle.process());
    ACTS_INFO("    position:     " << particle.position().transpose() / 1_mm
                                   << " mm");
    ACTS_INFO("    direction:    " << particle.direction().transpose());
    ACTS_INFO("    time:         " << particle.time() / 1_ns << " ns");
    ACTS_INFO("    |p|:          " << particle.absoluteMomentum() / 1_GeV
                                   << " GeV");
  }
  return ProcessCode::SUCCESS;
} 