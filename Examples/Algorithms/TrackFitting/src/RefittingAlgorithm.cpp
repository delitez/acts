// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/RefittingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"

#include <algorithm>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

namespace ActsExamples {

RefittingAlgorithm::RefittingAlgorithm(
    Config config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("RefittingAlgorithm", std::move(logger)),
      m_cfg(std::move(config)) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode RefittingAlgorithm::execute(const AlgorithmContext& ctx) const {
  const auto& inputTracks = m_inputTracks(ctx);

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  if (m_cfg.addBeamSpotMeasurement) {
    const Acts::Vector3 beamSpotCenter{0., 0., 0.};

    auto beamSpotVectorTrackStateContainer =
        std::make_shared<Acts::VectorMultiTrajectory>();
    auto beamSpotTrackState =
        beamSpotVectorTrackStateContainer->makeTrackState();
    // set it to the reference surface of the track and not to beamSpot
    ACTS_VERBOSE("beamSpot track state reference surface: "
                 << beamSpotTrackState.hasReferenceSurface());

    const Acts::Vector2 beamSpotMeasValue{0., 0.};
    Acts::SquareMatrix2 inflatedCov =
        Acts::SquareMatrix2::Zero();  //* 12.5 * Acts::UnitConstants::um;
    inflatedCov(0, 0) =
        12.5 *
        Acts::UnitConstants::um;  // 17.68 * Acts::UnitConstants::um; //156.25 *
                                  // Acts::UnitConstants::um;   // square
                                  // this 12.5 (used to be 17.68)
    inflatedCov(1, 1) =
        55.5 *
        Acts::UnitConstants::mm;  // 78.49 * Acts::UnitConstants::mm; //3080.25
                                  // * Acts::UnitConstants::mm;   // square
                                  // this 55.5 (used to be 3080.25)

    if (inputTracks.size() == 0) {
      ACTS_INFO("Input tracks collection is empty");
      return ProcessCode::SKIP;
    }
    auto trackRefSurfacePtr =
        inputTracks.at(0).referenceSurface().getSharedPtr();
    beamSpotTrackState.setReferenceSurface(trackRefSurfacePtr);
    beamSpotTrackState.allocateCalibrated(beamSpotMeasValue, inflatedCov);

    ACTS_VERBOSE(
        "Set uncalibrated source link for beamSpot track state with surface "
        << beamSpotTrackState.referenceSurface().geometryId());

    Acts::SourceLink testSL{42};
    beamSpotTrackState.setUncalibratedSourceLink(std::move(testSL));
    ACTS_VERBOSE("Get uncalibrated source link for beamSpot track state ");
    Acts::SourceLink uncalibSL = beamSpotTrackState.getUncalibratedSourceLink();

    auto beamSpotConstVectorTrackStateContainer =
        std::make_shared<Acts::ConstVectorMultiTrajectory>(
            std::move(*beamSpotVectorTrackStateContainer));

    auto beamSpotConstTrackState =
        beamSpotConstVectorTrackStateContainer->getTrackState(
            beamSpotTrackState.index());
    Acts::SourceLink uncalibSLconst =
        beamSpotConstTrackState.getUncalibratedSourceLink();
    ACTS_VERBOSE(
        "Got uncalibrated source link for beamSpot const track state ");

    RefittingCalibrator::RefittingSourceLink beamSpotSL{
        beamSpotConstTrackState};
  }

  // Perform the fit for each input track
  std::vector<Acts::SourceLink> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;
  RefittingCalibrator calibrator;

  auto itrack = 0ul;
  for (const auto& track : inputTracks) {
    // Check if you are not in picking mode
    ++itrack;
    if (m_cfg.pickTrack > -1 &&
        static_cast<std::size_t>(m_cfg.pickTrack) != itrack - 1) {
      continue;
    }

    if (!track.hasReferenceSurface()) {
      ACTS_VERBOSE("Skip track " << itrack << ": missing ref surface");
      continue;
    }

    TrackFitterFunction::GeneralFitterOptions options{
        ctx.geoContext,
        ctx.magFieldContext,
        ctx.calibContext,
        &track.referenceSurface(),
        Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext),
        true};

    Acts::BoundTrackParameters initialParams(
        track.referenceSurface().getSharedPtr(), track.parameters(),
        track.covariance(), track.particleHypothesis());

    if (initialParams.covariance()) {
      for (auto i = 0ul; i < m_cfg.initialVarInflation.size(); ++i) {
        (*initialParams.covariance())(i, i) *= m_cfg.initialVarInflation.at(i);
      }
    }

    trackSourceLinks.clear();
    surfSequence.clear();

    for (auto state : track.trackStatesReversed()) {
      surfSequence.push_back(&state.referenceSurface());

      if (!state.hasCalibrated()) {
        continue;
      }

      auto sl = RefittingCalibrator::RefittingSourceLink{state};
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }

    if (surfSequence.empty()) {
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    if (m_cfg.addBeamSpotMeasurement) {
      beamSpotSL =
          RefittingCalibrator::RefittingSourceLink{beamSpotConstTrackState};
      trackSourceLinks.push_back(Acts::SourceLink{beamSpotSL});
      surfSequence.push_back(&beamSpotConstTrackState.referenceSurface());
    }

    std::ranges::reverse(surfSequence);

    ACTS_VERBOSE("Initial parameters: " << track.parameters().transpose());

    ACTS_DEBUG("Invoke direct fitter for track " << itrack);
    auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
                               calibrator, surfSequence, tracks);

    if (result.ok()) {
      // Get the fit output object
      const auto& refittedTrack = result.value();
      if (refittedTrack.hasReferenceSurface()) {
        ACTS_VERBOSE("Refitted parameters for track " << itrack);
        ACTS_VERBOSE("  " << track.parameters().transpose());
        ACTS_VERBOSE("Measurements: " << refittedTrack.nMeasurements());
        ACTS_VERBOSE("Outliers: " << refittedTrack.nOutliers());
      } else {
        ACTS_DEBUG("No refitted parameters for track " << itrack);
      }
    } else {
      ACTS_WARNING("Fit failed for event "
                   << ctx.eventNumber << " track " << itrack << " with error: "
                   << result.error() << ", " << result.error().message());
    }
    ++itrack;
  }

  ACTS_DEBUG("Fitted tracks: " << trackContainer->size());

  if (logger().doPrint(Acts::Logging::DEBUG)) {
    std::stringstream ss;
    trackStateContainer->statistics().toStream(ss);
    ACTS_DEBUG(ss.str());
  }

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracks(ctx, std::move(constTracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
