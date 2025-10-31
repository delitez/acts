// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <mutex>
#include <ostream>
#include <ranges>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/Print.h>
#include <boost/container/flat_map.hpp>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsExamples {

TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }
  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
  if (m_cfg.doJetLabeling && !m_cfg.inputHepMC3Event.has_value()) {
    throw std::invalid_argument("Input HepMC3 event is not configured");
  }
  m_inputHepMC3Event.initialize(m_cfg.inputHepMC3Event.value());
  m_outputJets.initialize(m_cfg.outputJets);
}

namespace {
ActsPlugins::FastJet::JetLabel jetLabelFromHadronType(
    Acts::HadronType hadronType) {
  using enum Acts::HadronType;
  switch (hadronType) {
    case BBbarMeson:
    case BottomMeson:
    case BottomBaryon:
      return ActsPlugins::FastJet::JetLabel::BJet;
    case CCbarMeson:
    case CharmedMeson:
    case CharmedBaryon:
      return ActsPlugins::FastJet::JetLabel::CJet;
    case StrangeMeson:
    case StrangeBaryon:
    case LightMeson:
    case LightBaryon:
      return ActsPlugins::FastJet::JetLabel::LightJet;
    default:
      return ActsPlugins::FastJet::JetLabel::Unknown;
  }
}

}  // namespace

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const auto& truthParticles = m_inputTruthParticles(ctx);

  // Initialize the output container
  std::vector<ActsPlugins::FastJet::TruthJet<TrackContainer>>
      outputJetContainer{};

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition =
      fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  int particleIndex = 0;
  for (const auto& particle : truthParticles) {
    fastjet::PseudoJet pseudoJet(particle.momentum().x(),
                                 particle.momentum().y(),
                                 particle.momentum().z(), particle.energy());

    pseudoJet.set_user_index(particleIndex);
    inputPseudoJets.push_back(pseudoJet);
    particleIndex++;
  }
  ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());
  // Run the jet clustering
  fastjet::ClusterSequence clusterSeq(inputPseudoJets, defaultJetDefinition);
  // Get the jets above a certain pt threshold
  std::vector<fastjet::PseudoJet> jets =
      sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));

  ACTS_DEBUG("Number of clustered jets: " << jets.size());

  // Convert the fastjet::PseudoJet objects back to
  // ActsPlugins::FastJet::TruthJet objects
  for (const auto& jet : jets) {
    Acts::Vector4 jetFourMom(jet.px(), jet.py(), jet.pz(), jet.e());
    ActsPlugins::FastJet::TruthJet<TrackContainer> truthJet(jetFourMom);
    outputJetContainer.push_back(truthJet);
  }

  m_outputJets(ctx, std::move(outputJetContainer));

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::finalize() {
  ACTS_INFO("Finalizing truth jet clustering");
  return ProcessCode::SUCCESS;
}

};  // namespace ActsExamples
