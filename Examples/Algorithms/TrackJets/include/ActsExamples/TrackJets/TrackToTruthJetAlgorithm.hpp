// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

    class TrackToTruthJetAlgorithm final : public IAlgorithm {
        public:
        struct Config {
            /// Input (truth) track collection
            std::string inputTracks;
            /// Input truth jets collection
            std::string inputTruthJets;
            /// Output track-to-truth jets collection
            std::string outputTrackToTruthJets;
        };

        /// Constructor of the truth jet algorithm
        ///
        /// @param config is the config struct to configure the algorithm
        /// @param level is the logging level
        TruthJetAlgorithm(Config config, Acts::Logging::Level level);

        /// Framework execute method of the truth jet algorithm
        ///
        /// @param ctx is the algorithm context
        /// @return a process code to steer the algporithm flow
        ActsExamples::ProcessCode execute(AlgorithmContext& ctx) const override;

        const Config& config() const { return m_config; }

        private:
        Config m_config;

        ReadDataHandle<> m_inputTracks;
        ReadDataHandle<> m_inputTruthJets;
        WriteDataHandle<> m_outputTrackToTruthJets;

    };
} // namespace ActsExamples