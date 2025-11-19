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

s = acts.examples.Sequencer(events=10, numThreads=1, logLevel=acts.logging.INFO)
outputDir = "/Users/delitez/atlas/acts-spack/ci-dependencies/jetAlg_output"

from acts.examples.odd import getOpenDataDetector

detector = getOpenDataDetector()
trackingGeometry = detector.trackingGeometry()
digiConfigFile = "/Users/delitez/atlas/acts-spack/ci-dependencies/acts/Examples/Configs/odd-digi-smearing-config.json"

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


rnd = acts.examples.RandomNumbers(seed=42)
outputDir = Path(outputDir)
out = outputDir / "particles_edm4hep.root"

addParticleGun(
            s,
            MomentumConfig(
                0 * u.GeV,
                10 * u.GeV,
                transverse=True,
            ),
            EtaConfig(-2.5, 2.5),
            PhiConfig(0.0, 360.0 * u.degree),
            ParticleConfig(
                1, acts.PdgParticle.eB0, randomizeCharge=True
            ),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns
                ),
            ),
            multiplicity=1,
            rnd=rnd,
        )

# addPythia8(
#     s,
#     nhard=1,
#     npileup=1,
#     hardProcess=["Top:qqbar2ttbar=on"],
#     vtxGen=acts.examples.GaussianVertexGenerator(
#         mean=acts.Vector4(0, 0, 0, 0),
#         stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
#     ),
#     rnd=rnd,
#     outputDirRoot=outputDir,
#     outputDirCsv=outputDir,
#     # writeHepMC3=outputDir,
# )

# addFatras(
#     s,
#     trackingGeometry,
#     field,
#     rnd=rnd,
#     enableInteractions=True,
# )


# addDigitization(
#     s,
#     trackingGeometry,
#     field,
#     digiConfigFile=digiConfigFile,
#     rnd=rnd,
# )

# addDigiParticleSelection(
#     s,
#     ParticleSelectorConfig(
#         pt=(0.9 * u.GeV, None),
#         measurements=(7, None),
#         removeNeutral=True,
#         removeSecondaries=True,
#     ),
# )

edm4hepParticleConverter = EDM4hepParticleOutputConverter(
    acts.logging.INFO,
    inputParticles="particles",
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

# s.addReader(
#     PodioReader(
#         level=acts.logging.DEBUG,
#         inputPath=str(out),
#         outputFrame="events",
#         category="events",
#     )
# )

# edm4hepReader = acts.examples.edm4hep.EDM4hepSimInputConverter(
#     inputFrame="events",
#     inputSimHits=[
#         "PixelBarrelReadout",
#         "PixelEndcapReadout",
#         "ShortStripBarrelReadout",
#         "ShortStripEndcapReadout",
#         "LongStripBarrelReadout",
#         "LongStripEndcapReadout",
#     ],
#     outputParticlesGenerator="particles_generated",
#     outputParticlesSimulation="particles_simulated",
#     outputSimHits="simhits",
#     outputSimVertices="vertices_truth",
#     dd4hepDetector=detector,
#     trackingGeometry=trackingGeometry,
#     sortSimHitsInTime=False,
#     particleRMax=1080 * u.mm,
#     particleZ=(-3030 * u.mm, 3030 * u.mm),
#     particlePtMin=150 * u.MeV,
#     level=acts.logging.DEBUG,
# )
# s.addAlgorithm(edm4hepReader)

# s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesSimulation)

# addSimParticleSelection(
#     s,
#     ParticleSelectorConfig(
#         rho=(0.0, 24 * u.mm),
#         absZ=(0.0, 1.0 * u.m),
#         eta=(-3.0, 3.0),
#         removeNeutral=True,
#     ),
# )


s.run()

# addTruthJetAlg(
#     s,
#     TruthJetConfig(
#         inputTruthParticles="particles_selected",
#         inputEDM4HepParticles="MCParticles",
#         outputJets="output_jets",
#         jetPtMin=20 * u.GeV,
#     ),
#     loglevel=acts.logging.INFO,
# )
