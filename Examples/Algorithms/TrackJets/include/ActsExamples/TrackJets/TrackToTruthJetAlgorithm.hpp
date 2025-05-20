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
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

#include <fastjet/PseudoJet.hh>

namespace ActsExamples {
struct AlgorithmContext;

class TrackToTruthJetAlgorithm : public IAlgorithm {
 public:
  struct Config {
    /// Input tracks collection.
    std::string inputTracks;
    /// Input jets collection.
    std::string inputJets;
    /// Output jets collection matched to tracks.
    std::string outputTrackJets;
  };

  TrackToTruthJetAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() const;

  const Config& config() const { return m_cfg; }

 private:
    Config m_cfg;
    ReadDataHandle<SimParticleContainer> m_inputTracks{
        this, "inputTracks"};
    ReadDataHandle<std::vector<fastjet::PseudoJet>> m_inputjets{
        this, "inputJets"};
    // FIX: A collection of tracks and jets - need to think about this
    WriteDataHandle<std::vector<fastjet::PseudoJet>> m_outputJets{
        this, "outputTrackJets"};
};

}  // namespace ActsExamples


