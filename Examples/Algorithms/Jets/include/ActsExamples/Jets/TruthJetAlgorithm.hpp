// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/FastJet/Jets.hpp"

#include <string>

namespace fastjet {
class PseudoJet;
}

namespace HepMC3 {
  class GenEvent;
}

namespace ActsExamples {
struct AlgorithmContext;

using TruthJetContainer =
    std::vector<ActsPlugins::FastJet::TruthJet<TrackContainer>>;

class TruthJetAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputTruthParticles;
    /// Output jets collection.
    std::string outputJets;
    /// Minimum jet pT.
    double jetPtMin = 20 * Acts::UnitConstants::GeV;
    /// Jet eta range
    std::pair<std::optional<double>, std::optional<double>> jetEtaRange = {
        std::nullopt, std::nullopt};
    /// Jet clustering radius
    double jetClusteringRadius = 0.4;
    /// Only cluster HS particles
    bool clusterHSParticlesOnly = true;
    /// input HepMC3 event
    std::optional<std::string> inputHepMC3Event;
    /// Do jet labeling
    bool doJetLabeling = true;
    /// Delta R for labeling
    double jetLabelingDeltaR = 0.4;
    /// Minimum hadron pT for labeling
    double jetLabelingHadronPtMin = 5 * Acts::UnitConstants::GeV;
    /// Only label HS hadrons
    bool jetLabelingHSHadronsOnly = true;
    /// Output for debugging
    bool debugCsvOutput = false;
  };

  TruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode initialize() override;
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<SimParticleContainer> m_inputTruthParticles{
      this, "inputTruthParticles"};
  ReadDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_inputHepMC3Event{
      this, "inputHepMC3Event"};
  WriteDataHandle<TruthJetContainer> m_outputJets{this, "outputJets"};

  /// Statistics for jets
  mutable std::atomic<std::size_t> m_numJets = 0;
  mutable std::atomic<std::size_t> m_numJetsAfterOverlapRemoval = 0;
  mutable std::atomic<std::size_t> m_numLightJets = 0;
  mutable std::atomic<std::size_t> m_numCJets = 0;
  mutable std::atomic<std::size_t> m_numBJets = 0;
};

}  // namespace ActsExamples
