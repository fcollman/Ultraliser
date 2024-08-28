/***************************************************************************************************
 * Copyright (c) 2016 - 2024
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include <Ultraliser.h>
#include <AppCommon.h>
#include <AppArguments.h>
#include "SkeletonizeNeuron.h"

#include <geometry/Circle.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< AppArguments > args = std::make_unique< AppArguments >(
        argc, argv, COPYRIGHT
        "This application reconstructs a high quality neuronal morphology skeleton (directed "
        "acyclic graph) from an input mesh model of the neuron. If this mesh contains spines, the "
        "application is capable of segmenting those spines and reconstruct high quality skeletons "
        "of the individual spines and providing some information on their types and locations. "
        "The application requires a neuronal mesh to create a valid morphology skeleton. "
        "The scale of the input mesh must be microns.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addNeuronalMorphologyExportArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addProcessingArguments();
    args->addGenericSkeletonizationArguments();
    args->addSkeletonizationAccelerationArgument();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMeshPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();
    options->verifyNeuronalMorphologyExportArguments();
    options->verifyProcessingArguments();

    // Initialize context, once everything is in place and all the options are verified
    options->initializeContext();

    // Return the executable options
    return options;
}


void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Load the input mesh of the neuron
    auto inputNeuronMesh = loadInputMesh(options);

    // Run the soma repair operations in case of missing slices
    runFixSomaSlicingArtifactsOperations(options, inputNeuronMesh);

    // Construct the neuron volume
    auto neuronVolume = createNeuronVolume(inputNeuronMesh, options, SILENT);

    // Project the volume created from the input mesh
    projectInputNeuronVolume(options, neuronVolume);

    // Create a NeuronSkeletonizer process
    auto skeletonizer = createNeuronSkeletonizer(options, neuronVolume);

    // Run the skeletonization operations
    runNeuronSkeletonizationOperations(options, skeletonizer);

    // Run the spine segmentation operations
    runSpineSegmentationOperations(options, skeletonizer, inputNeuronMesh);

    // Create high quality mesh
    runHighQualityMeshGeneration(options, skeletonizer, inputNeuronMesh);

    // Write the statistics
    runStatisticsGenerator(options, skeletonizer);
}

}

int main(int argc , const char** argv)
{
    TIMER_SET;

    Ultraliser::run(argc, argv);

    LOG_STATUS_IMPORTANT("Ultralization Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    ULTRALISER_DONE;
}
