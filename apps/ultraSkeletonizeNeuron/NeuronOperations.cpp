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

#include "SkeletonizeNeuron.h"
#include "AppCommon.h"

// Defines
#define NEURON_SMOOTHING_ITERATIONS 25

namespace Ultraliser
{

Volume* createNeuronVolume(Mesh* neuronMesh, const AppOptions* options, const bool verbose)
{
    // Create the volume from the mesh
    auto neuronVolume = createVolumeGrid(neuronMesh, options, verbose);

    // Adaptive and conservative Voxelization
    neuronVolume->surfaceVoxelization(neuronMesh, verbose, false, 1.0);
    neuronVolume->solidVoxelization(options->voxelizationAxis, SILENT);

    // Remove the border voxels that span less than half the voxel
    // TODO: VERIFY neuronVolume->surfaceVoxelization(neuronMesh, false, false, 0.5);
    auto bordeVoxels = neuronVolume->searchForBorderVoxels(verbose);
    for (size_t i = 0; i < bordeVoxels.size(); ++i)
    {
        for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
        {
            auto voxel = bordeVoxels[i][j];
            neuronVolume->clear(voxel.x(), voxel.y(), voxel.z());
        }
        bordeVoxels[i].clear();
    }
    bordeVoxels.clear();
    neuronVolume->surfaceVoxelization(neuronMesh, false, false, 0.5);

    // Return the volume
    return neuronVolume;
}

Mesh* reconstructNeuronMeshFromVolume(Volume* neuronVolume, AppOptions* options, const bool verbose)
{
    // Construct the mesh using the DMC technique
    auto reconstructedNeuronMesh = DualMarchingCubes::generateMeshFromVolume(neuronVolume);

    // Smooth the resulting surface mesh
    reconstructedNeuronMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, verbose);

    // Return a pointer to the resulting neuron
    return reconstructedNeuronMesh;
}

Mesh* remeshNeuron(Mesh* inputNeuronMesh, AppOptions* options, const bool verbose)
{
    // Compute the bounding box of the input neuron mesh
    Vector3f pMinInput, pMaxInput;
    inputNeuronMesh->computeBoundingBox(pMinInput, pMaxInput);
    const auto& meshBoundingBox = pMaxInput - pMinInput;

    // Compute the resolution of the volume
    const auto largestDimension = meshBoundingBox.getLargestDimension();
    size_t resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);

    // Construct the volume from the input mesh
    auto volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                             VolumeGrid::getType(options->volumeType), verbose);

    // Apply surface and solid voxelization to the input neuron mesh
    volume->surfaceVoxelization(inputNeuronMesh, false, false, 1.0);
    volume->solidVoxelization(options->voxelizationAxis);

    // Remove the border voxels that span less than half the voxel
    auto bordeVoxels = volume->searchForBorderVoxels(verbose);
    for (size_t i = 0; i < bordeVoxels.size(); ++i)
    {
        for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
        {
            auto voxel = bordeVoxels[i][j]; volume->clear(voxel.x(), voxel.y(), voxel.z());
        }
        bordeVoxels[i].clear();
    }
    bordeVoxels.clear();
    volume->surfaceVoxelization(inputNeuronMesh, false, false, 0.5);

    // Construct the mesh using the DMC technique
    auto reconstructedNeuronMesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Smooth the resulting surface mesh
    reconstructedNeuronMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, verbose);

    // Return a pointer to the resulting neuron
    return reconstructedNeuronMesh;
}

void projectInputNeuronVolume(const AppOptions* options, const Volume* neuronVolume)
{
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        neuronVolume->project(options->projectionPrefix + "-input",
                             options->projectXY, options->projectXZ, options->projectZY);
    }
}

NeuronSkeletonizer* createNeuronSkeletonizer(const AppOptions* options, Volume* neuronVolume)
{
    // Create a skeletonization object
    return new NeuronSkeletonizer(neuronVolume,
                                  options->removeSpinesFromSkeleton,
                                  options->useAccelerationStructures,
                                  options->debugSkeletonization,
                                  options->debuggingPrefix);
}

void projectSkeletonVolume(const AppOptions* options, const NeuronSkeletonizer* skeletonizer)
{
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        const std::string prefix = options->projectionPrefix + SKELETON_SUFFIX;
        skeletonizer->getVolume()->project(
                    prefix, options->projectXY, options->projectXZ, options->projectZY);
    }
}

void exportNeuronMorphology(const AppOptions* options,
                            NeuronSkeletonizer* skeletonizer)
{
    // Export the SWC file of the neuron
    if (options->exportSWCNeuron)
    {
        skeletonizer->exportSWCFile(options->morphologyPrefix, options->resampleSkeleton, VERBOSE);
    }
}

void runNeuronSkeletonizationOperations(const AppOptions* options,
                                        NeuronSkeletonizer* skeletonizer)
{
    // Initialize the skeletonizer
    skeletonizer->initialize(VERBOSE);

    // Skeletonize the volume to obtain the center-lines
    skeletonizer->skeletonizeVolumeToCenterLines(VERBOSE);

    // Project the center-lines of the skeleton before constructing the graph
    projectSkeletonVolume(options, skeletonizer);

    // Construct the neuron graph from the volume
    skeletonizer->constructGraph(VERBOSE);

    // Segment the different components of the graph
    skeletonizer->segmentComponents(VERBOSE);

    // Export the morphology of the neuron
    exportNeuronMorphology(options, skeletonizer);
}

}
