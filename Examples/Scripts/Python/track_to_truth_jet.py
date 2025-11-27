# #!/usr/bin/env python3

# from pathlib import Path
# from typing import Optional

# import acts
# import acts.examples

# from truth_tracking_kalman import runTruthTrackingKalman

# u = acts.UnitConstants

# from acts.examples.simulation import (
#         addParticleGun,
#         ParticleConfig,
#         EtaConfig,
#         PhiConfig,
#         MomentumConfig,
#         TruthJetConfig,
#         # TrackToTruthJetConfig,
#         addFatras,
#         addPythia8,
#         addTruthJetAlg,
#         # addTrackToTruthJetAlg,
#         addDigitization,
#         ParticleSelectorConfig,
#         addDigiParticleSelection,
# )
# from acts.examples.reconstruction import (
#         addSeeding,
#         SeedingAlgorithm,
#         addKalmanTracks,
# )

# s = acts.examples.Sequencer(
#         events=1, numThreads=1, logLevel=acts.logging.INFO
# )
# outputDir = "/Users/delitez/atlas/acts-spack/ci-dependencies/truth_jet_test_output"

# from acts.examples.odd import getOpenDataDetector

# detector = getOpenDataDetector()
# trackingGeometry = detector.trackingGeometry()
# digiConfigFile = "/Users/delitez/atlas/acts-spack/ci-dependencies/acts/Examples/Configs/odd-digi-smearing-config.json"

# field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


# rnd = acts.examples.RandomNumbers(seed=42)
# outputDir = Path(outputDir)

# addPythia8(
#         s,
#         nhard=1,
#         npileup=1,
#         hardProcess=["Top:qqbar2ttbar=on"],
#         vtxGen=acts.examples.GaussianVertexGenerator(
#             mean=acts.Vector4(0, 0, 0, 0),
#             stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
#         ),
#         rnd=rnd,
#         outputDirRoot=outputDir,
#         outputDirCsv=outputDir,
#         writeHepMC3=outputDir
#         )

# # addFatras(
# #             s,
# #             trackingGeometry,
# #             field,
# #             rnd=rnd,
# #             enableInteractions=True,
# # )


# # addDigitization(
# #         s,
# #         trackingGeometry,
# #         field,
# #         digiConfigFile=digiConfigFile,
# #         rnd=rnd,
# #     )

# # addDigiParticleSelection(
# #         s,
# #         ParticleSelectorConfig(
# #             pt=(0.9 * u.GeV, None),
# #             measurements=(7, None),
# #             removeNeutral=True,
# #             removeSecondaries=True,
# #         ),
# #     )

# # addSeeding(
# #         s,
# #         trackingGeometry,
# #         field,
# #         rnd=rnd,
# #         inputParticles="particles_generated",
# #         seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
# #         particleHypothesis=acts.ParticleHypothesis.muon,
# #     )

# # reverseFilteringMomThreshold=0 * u.GeV

# # addKalmanTracks(
# #         s,
# #         trackingGeometry,
# #         field,
# #         reverseFilteringMomThreshold,
# #     )

# # s.addAlgorithm(
# #         acts.examples.TrackSelectorAlgorithm(
# #             level=acts.logging.INFO,
# #             inputTracks="tracks",
# #             outputTracks="selected-tracks",
# #             selectorConfig=acts.TrackSelector.Config(
# #                 minMeasurements=7,
# #             ),
# #         )
# #     )
# # s.addWhiteboardAlias("tracks", "selected-tracks")

# # s.addWriter(
# #         acts.examples.RootTrackStatesWriter(
# #             level=acts.logging.INFO,
# #             inputTracks="tracks",
# #             inputParticles="particles_selected",
# #             inputTrackParticleMatching="track_particle_matching",
# #             inputSimHits="simhits",
# #             inputMeasurementSimHitsMap="measurement_simhits_map",
# #             filePath=str(outputDir / "trackstates_kf.root"),
# #         )
# #     )

# # s.addWriter(
# #         acts.examples.RootTrackSummaryWriter(
# #             level=acts.logging.INFO,
# #             inputTracks="tracks",
# #             inputParticles="particles_selected",
# #             inputTrackParticleMatching="track_particle_matching",
# #             filePath=str(outputDir / "tracksummary_kf.root"),
# #         )
# #     )

# # s.addWriter(
# #         acts.examples.TrackFitterPerformanceWriter(
# #             level=acts.logging.INFO,
# #             inputTracks="tracks",
# #             inputParticles="particles_selected",
# #             inputTrackParticleMatching="track_particle_matching",
# #             filePath=str(outputDir / "performance_kf.root"),
# #         )
# #     )

# addTruthJetAlg(
#     s,
#     TruthJetConfig(
#         inputTruthParticles="particles_generated",
#         inputHepMC3Event="pythia8-event",
#         outputJets="output_jets",
#     ),
#         loglevel=acts.logging.DEBUG
# )

# # addTrackToTruthJetAlg(
# #     s,
# #     TrackToTruthJetConfig(
# #         inputTracks="tracks",
# #         inputJets="truth_jets",
# #         outputTrackJets="track_jets",
# #         maxDeltaR=0.4
# #     ),
# #     loglevel=acts.logging.DEBUG
# # )

# s.run()

#!/usr/bin/env python3

import os
from pathlib import Path
from typing import Optional

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    TruthJetConfig,
    addFatras,
    addPythia8,
    addTruthJetAlg,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
    addSimParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
)

from acts.examples.edm4hep import (
    EDM4hepParticleOutputConverter,
    PodioWriter,
    PodioReader,
)

s = acts.examples.Sequencer(events=1, numThreads=1, logLevel=acts.logging.INFO)
outputDir = "/Users/delitez/atlas/acts-spack/ci-dependencies/jetAlg_output"

from acts.examples.odd import getOpenDataDetector

detector = getOpenDataDetector()
trackingGeometry = detector.trackingGeometry()
digiConfigFile = "/Users/delitez/atlas/acts-spack/ci-dependencies/acts/Examples/Configs/odd-digi-smearing-config.json"

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


rnd = acts.examples.RandomNumbers(seed=42)
outputDir = Path(outputDir)
out = outputDir / "particles_edm4hep.root"

addPythia8(
    s,
    nhard=1,
    npileup=1,
    hardProcess=["Top:qqbar2ttbar=on"],
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
    ),
    rnd=rnd,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir,
    writeHepMC3=outputDir,
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    enableInteractions=True,
)


addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=digiConfigFile,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(0.9 * u.GeV, None),
        measurements=(7, None),
        removeNeutral=True,
        removeSecondaries=True,
    ),
)

edm4hepParticleConverter = EDM4hepParticleOutputConverter(
    acts.logging.INFO,
    inputParticles="particles_generated_selected",
    outputParticles="MCParticles",
)

s.addAlgorithm(edm4hepParticleConverter)

s.addWriter(
    PodioWriter(
        level=acts.logging.INFO,
        outputPath=str(out),
        category="events",
        collections=edm4hepParticleConverter.collections,
        separateFilesPerThread=True,
    )
)

s.addReader(
    PodioReader(
        level=acts.logging.DEBUG,
        inputPath=str(out),
        outputFrame="events",
        category="events",
    )
)

edm4hepReader = acts.examples.edm4hep.EDM4hepSimInputConverter(
    inputFrame="events",
    inputSimHits=[
        "PixelBarrelReadout",
        "PixelEndcapReadout",
        "ShortStripBarrelReadout",
        "ShortStripEndcapReadout",
        "LongStripBarrelReadout",
        "LongStripEndcapReadout",
    ],
    outputParticlesGenerator="particles_generated",
    outputParticlesSimulation="particles_simulated",
    outputSimHits="simhits",
    outputSimVertices="vertices_truth",
    dd4hepDetector=detector,
    trackingGeometry=trackingGeometry,
    sortSimHitsInTime=False,
    particleRMax=1080 * u.mm,
    particleZ=(-3030 * u.mm, 3030 * u.mm),
    particlePtMin=150 * u.MeV,
    level=acts.logging.DEBUG,
)
s.addAlgorithm(edm4hepReader)

s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesSimulation)

addSimParticleSelection(
    s,
    ParticleSelectorConfig(
        rho=(0.0, 24 * u.mm),
        absZ=(0.0, 1.0 * u.m),
        eta=(-3.0, 3.0),
        removeNeutral=True,
    ),
)

addTruthJetAlg(
    s,
    TruthJetConfig(
        inputTruthParticles="particles_selected",
        inputHepMC3Event="pythia8-event",
        inputEDM4HepParticles="MCParticles",
        outputJets="output_jets",
        jetPtMin=20 * u.GeV,
    ),
    loglevel=acts.logging.INFO,
)

s.run()

