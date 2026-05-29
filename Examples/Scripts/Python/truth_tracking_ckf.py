#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import acts
import acts.examples
from acts.examples.edm4hep import (
    EDM4hepParticleOutputConverter,
    PodioWriter,
    PodioReader,
    EDM4hepSimInputConverter,
)

u = acts.UnitConstants


def runTruthTrackingCKF(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputHitsPath: Optional[Path] = None,
    decorators=[],
    generatedParticleType: acts.PdgParticle = acts.PdgParticle.eMuon,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=100.0,
    numParticles=1,
    linkForward: bool = False,
    useJosephFormulation: bool = False,
    s: acts.examples.Sequencer = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addFatras,
        addPythia8,
        addDigitization,
        addSimParticleSelection,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )

    from acts.examples.root import (
        RootParticleReader,
        RootSimHitReader,
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootTrackFitterPerformanceWriter,
    )

    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TrackSmearingSigmas,
        addKalmanTracks,
        addCKFTracks,
        TrackSelectorConfig,
        CkfConfig,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    logger = acts.getDefaultLogger("Truth tracking example", acts.logging.INFO)

    if inputHitsPath is not [None]:

        from acts.examples.odd import getOpenDataDetector

        detector = getOpenDataDetector()
        trackingGeometry = detector.trackingGeometry()
        digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"
        inputFile = Path(inputHitsPath)
        s.addReader(
            PodioReader(
                level=acts.logging.INFO,
                inputPath=str(inputFile),
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
            level=acts.logging.INFO,
        )
        s.addAlgorithm(edm4hepReader)

        s.addWhiteboardAlias(
            "particles", edm4hepReader.config.outputParticlesSimulation
        )

        addSimParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                removeNeutral=True,
            ),
        )
    else:
        if inputParticlePath is None:
            addPythia8(
                s,
                hardProcess=["Top:qqbar2ttbar=on"],
                npileup=200,
                vtxGen=acts.examples.GaussianVertexGenerator(
                    mean=acts.Vector4(0, 0, 0, 0),
                    stddev=acts.Vector4(
                        0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                    ),
                ),
                rnd=rnd,
            )
            # addParticleGun(
            #     s,
            #     ParticleConfig(
            #         num=numParticles, pdg=generatedParticleType, randomizeCharge=True
            #     ),
            #     EtaConfig(-3.0, 3.0, uniform=True),
            #     MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            #     PhiConfig(0.0, 360.0 * u.degree),
            #     vtxGen=acts.examples.GaussianVertexGenerator(
            #         mean=acts.Vector4(0, 0, 0, 0),
            #         stddev=acts.Vector4(0, 0, 0, 0),
            #     ),
            #     multiplicity=1,
            #     rnd=rnd,
            # )
        else:
            logger.info("Reading particles from {}", inputParticlePath.resolve())
            assert inputParticlePath.exists()
            s.addReader(
                RootParticleReader(
                    level=acts.logging.INFO,
                    filePath=str(inputParticlePath.resolve()),
                    outputParticles="particles_generated",
                )
            )
            s.addWhiteboardAlias("particles", "particles_generated")

        if inputHitsPath is None:
            addFatras(
                s,
                trackingGeometry,
                field,
                rnd=rnd,
                enableInteractions=True,
            )
        else:
            logger.info("Reading hits from {}", inputHitsPath.resolve())
            assert inputHitsPath.exists()
            s.addReader(
                RootSimHitReader(
                    level=acts.logging.INFO,
                    filePath=str(inputHitsPath.resolve()),
                    outputSimHits="simhits",
                )
            )
            s.addWhiteboardAlias("particles_simulated_selected", "particles_generated")

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

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        particleHypothesis=acts.ParticleHypothesis.muon,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        geoSelectionConfigFile=srcdir / "Examples/Configs/odd-seeding-config.json",
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(
            pt=(1.0 * u.GeV, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=7,
            maxHoles=2,
            maxOutliers=2,
        ),
        CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=2,
            # seedDeduplication=True,
            # stayOnSeed=True,
            pixelVolumes=[16, 17, 18],
            stripVolumes=[23, 24, 25],
            maxPixelHoles=1,
            maxStripHoles=2,
            constrainToVolumes=[
                2,  # beam pipe
                32,
                4,  # beam pip gap
                16,
                17,
                18,  # pixel
                20,  # PST
                23,
                24,
                25,  # short strip
                26,
                8,  # long strip gap
                28,
                29,
                30,  # long strip
            ],
        ),
        outputDirRoot=outputDir,
        writeCovMat=True,
    )

    s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
    s.addWhiteboardAlias("tracks", "selected-tracks")

    s.addWriter(
        RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_ckf.root"),
        )
    )

    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_ckf.root"),
        )
    )

    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_ckf.root"),
        )
    )

    return s


if "__main__" == __name__:
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the truth-tracking CKF example on one or more EDM4hep files"
    )
    parser.add_argument(
        "--edm4hep",
        nargs="+",
        type=Path,
        default=[Path("/Users/delitez/atlas/ttbar_edm4hep/0/edm4hep.root")],
        help="One or more EDM4hep input files",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path.cwd(),
        help="Output directory",
    )

    cli_args = parser.parse_args()

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    # ODD
    from acts.examples.odd import getOpenDataDetector

    detector = getOpenDataDetector()
    trackingGeometry = detector.trackingGeometry()
    digiConfigFile = srcdir / "Examples/Configs/odd-digi-smearing-config.json"

    # field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    field = detector.field
    # field = acts.SolenoidBField(
    #     radius=1200 * u.mm,
    #     length=6000 * u.mm,
    #     bMagCenter=3 * u.T,
    #     nCoils=1194,
    # )

    if cli_args.edm4hep is not None:
        for edm4hepInput in cli_args.edm4hep:
            outputDir = (
                cli_args.output
                if len(cli_args.edm4hep) == 1
                else cli_args.output / edm4hepInput.stem
            )

            runTruthTrackingCKF(
                trackingGeometry=trackingGeometry,
                field=field,
                digiConfigFile=digiConfigFile,
                detector=detector,
                outputDir=outputDir,
                inputHitsPath=edm4hepInput,
            ).run()
    else:
        runTruthTrackingCKF(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=digiConfigFile,
            detector=detector,
            outputDir=outputDir,
        ).run()
