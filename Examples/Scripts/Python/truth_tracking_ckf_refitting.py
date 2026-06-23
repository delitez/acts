#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional

import acts
import acts.examples
from acts.examples.root import (
    RootTrackStatesWriter,
    RootTrackSummaryWriter,
    RootTrackFitterPerformanceWriter,
)

from truth_tracking_ckf import runTruthTrackingCKF

u = acts.UnitConstants


def runRefittingCkf(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputHitsPath: Optional[list[Path]] = None,
    multipleScattering: bool = True,
    energyLoss: bool = True,
    reverseFilteringMomThreshold=float("inf"),
    reverseFilteringCovarianceScaling=100.0,
    useJosephFormulation: bool = False,
    s: acts.examples.Sequencer = None,
):
    outputDir.mkdir(parents=True, exist_ok=True)
    s = runTruthTrackingCKF(
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDir=outputDir,
        inputHitsPath=inputHitsPath,
        reverseFilteringMomThreshold=0 * u.GeV,  # use direct smoothing
        reverseFilteringCovarianceScaling=reverseFilteringCovarianceScaling,
        useJosephFormulation=useJosephFormulation,
        s=s,
    )

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "reverseFilteringCovarianceScaling": reverseFilteringCovarianceScaling,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": acts.logging.INFO,
        "chi2Cut": float("inf"),
        "useJosephFormulation": useJosephFormulation,
    }

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            level=acts.logging.INFO,
            inputTracks="ckf_tracks",
            outputTracks="ckf_refit_tracks",
            # inputTracks="kf_tracks",
            # outputTracks="kf_refit_tracks",
            initialVarInflation=6 * [100.0],
            fit=acts.examples.makeKalmanFitterFunction(
                trackingGeometry, field, **kalmanOptions
            ),
            beamSpotConstraint=acts.SquareMatrix2(
                [[0.0125 * u.mm, 0], [0, 55.5 * u.mm]]
            ),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="ckf_refit_tracks",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    s.addWriter(
        RootTrackStatesWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_ckf_refit.root"),
        )
    )

    s.addWriter(
        RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "tracksummary_ckf_refit.root"),
        )
    )

    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="ckf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "performance_ckf_refit.root"),
        )
    )

    return s


def collectEdm4hepInputs(inputPaths: list[Path]) -> list[tuple[Path, Path]]:
    resolvedInputs: list[tuple[Path, Path]] = []

    for inputPath in inputPaths:
        if inputPath.is_dir():
            for edm4hepFile in sorted(inputPath.rglob("edm4hep.root")):
                if edm4hepFile.is_file():
                    relativeParent = edm4hepFile.parent.relative_to(inputPath)
                    outputSuffix = (
                        relativeParent
                        if relativeParent != Path(".")
                        else edm4hepFile.parent.name
                    )
                    resolvedInputs.append((edm4hepFile, outputSuffix))
        else:
            outputSuffix = (
                inputPath.parent.name if inputPath.parent.name else inputPath.stem
            )
            resolvedInputs.append((inputPath, Path(outputSuffix)))

    return resolvedInputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run truth tracking CKF refitting on one or more EDM4hep files"
    )
    parser.add_argument(
        "--edm4hep",
        nargs="+",
        type=Path,
        default=[None],
        help="One or more EDM4hep input files or directories containing edm4hep.root files",
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

    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    # Get detector and field
    geoDir = getOpenDataDetectorDirectory()

    # Load material map

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"

    oddSeedingSel = geoDir / "config/odd-seeding-config.json"
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    # Get detector
    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = detector.field

    # field = acts.SolenoidBField(
    #     radius=1200 * u.mm,
    #     length=6000 * u.mm,
    #     bMagCenter=3 * u.T,
    #     nCoils=1194,
    # )

    # solenoid = acts.SolenoidBField(
    #     radius=1200 * u.mm,
    #     length=6000 * u.mm,
    #     bMagCenter=3 * u.T,
    #     nCoils=1194,
    # )

    # field = acts.solenoidFieldMap(
    #     rlim=(0, 1200 * u.mm),
    #     zlim=(-5000 * u.mm, 5000 * u.mm),
    #     nbins=(50, 50),
    #     field=solenoid,
    # )

    if cli_args.edm4hep:
        runRefittingCkf(
            trackingGeometry=trackingGeometry,
            field=field,
            digiConfigFile=oddDigiConfig,
            outputDir=cli_args.output,
            inputHitsPath=cli_args.edm4hep,
        ).run()
