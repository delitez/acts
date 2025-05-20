// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackJets/TrackToTruthJetAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <ostream>
#include <stdexcept>


namespace ActsExamples {

    TrackToTruthJetAlgorithm::TrackToTruthJetAlgorithm(const Config& cfg,
                                                       Acts::Logging::Level lvl)
        : IAlgorithm("TrackToTruthJetAlgorithm", lvl), m_cfg(cfg) {
        if (m_cfg.inputTracks.empty()) {
            throw std::invalid_argument("Input tracks collection is not configured");
        }
        if (m_cfg.inputJets.empty()) {
            throw std::invalid_argument("Input jets collection is not configured");
        }

        m_inputTracks.initialize(m_cfg.inputTracks);
        m_inputjets.initialize(m_cfg.inputJets);
        m_outputJets.initialize(m_cfg.outputTrackJets);
    }

} // namespace ActsExamples