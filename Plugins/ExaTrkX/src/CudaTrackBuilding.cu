// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/CudaTrackBuilding.hpp"
#include "Acts/Plugins/ExaTrkX/detail/ConnectedComponents.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <c10/cuda/CUDAGuard.h>
#include <c10/cuda/CUDAStream.h>
#include <torch/torch.h>

namespace Acts {

std::vector<std::vector<int>> CudaTrackBuilding::operator()(
    std::any /*nodes*/, std::any edges, std::any weights,
    std::vector<int>& spacepointIDs, const ExecutionContext& execContext) {
  ACTS_VERBOSE("Start CUDA track building");
  c10::cuda::CUDAStreamGuard guard(execContext.stream.value());

  const auto edgeTensor = std::any_cast<torch::Tensor>(edges).to(torch::kCUDA);
  assert(edgeTensor.size(0) == 2);

  const auto numSpacepoints = spacepointIDs.size();
  const auto numEdges = static_cast<std::size_t>(edgeTensor.size(1));

  if (numEdges == 0) {
    ACTS_WARNING("No edges remained after edge classification");
    return {};
  }

  auto stream = execContext.stream->stream();

  auto cudaSrcPtr = edgeTensor.data_ptr<std::int64_t>();
  auto cudaTgtPtr = edgeTensor.data_ptr<std::int64_t>() + numEdges;

  auto ms = [](auto t0, auto t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
        .count();
  };

  int* cudaLabels;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaLabels, numSpacepoints * sizeof(int), stream));

  auto t0 = std::chrono::high_resolution_clock::now();
  std::size_t numberLabels = detail::connectedComponentsCuda(
      numEdges, cudaSrcPtr, cudaTgtPtr, numSpacepoints, cudaLabels, stream,
      m_cfg.useOneBlockImplementation);
  auto t1 = std::chrono::high_resolution_clock::now();
  ACTS_DEBUG("Connected components took " << ms(t0, t1) << " ms");

  // TODO not sure why there is an issue that is not detected in the unit tests
  numberLabels += 1;

  std::vector<int> trackLabels(numSpacepoints);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(trackLabels.data(), cudaLabels,
                                  numSpacepoints * sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cudaLabels, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("Found " << numberLabels << " track candidates");

  std::vector<std::vector<int>> trackCandidates(numberLabels);

  for (const auto [label, id] : Acts::zip(trackLabels, spacepointIDs)) {
    trackCandidates[label].push_back(id);
  }

  return trackCandidates;
}

}  // namespace Acts
