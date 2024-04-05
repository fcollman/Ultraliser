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

Mesh* remeshSpine(Mesh* inputSpineMesh, const float voxelsPerMicron = 50,
                  const bool verbose = SILENT)
{
    // Compute the bounding box of the input neuron mesh
    Vector3f pMinInput, pMaxInput;
    inputSpineMesh->computeBoundingBox(pMinInput, pMaxInput);
    const auto& meshBoundingBox = pMaxInput - pMinInput;

    if (meshBoundingBox.x() < 1e-16 || meshBoundingBox.y() < 1e-16 || meshBoundingBox.z() < 1e-16)
    {
        return nullptr;
    }

    // Compute the resolution of the volume
    const auto largestDimension = meshBoundingBox.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, 0.1,
                                VOLUME_TYPE::BIT, SILENT);

    // Apply surface and solid voxelization to the input neuron mesh
    volume->surfaceVoxelization(inputSpineMesh, SILENT, false, 1.0);
    volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ, SILENT);
    volume->surfaceVoxelization(inputSpineMesh, SILENT, false, 0.5);

    // Construct the mesh using the DMC technique
    auto reconstructedSpineMesh = DualMarchingCubes::generateMeshFromVolume(volume, SILENT);

    // Smooth the resulting surface mesh
    reconstructedSpineMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, SILENT);

    // Return a pointer to the resulting neuron
    return reconstructedSpineMesh;
}

void exportSpineMeshes(NeuronSkeletonizer* skeletonizer,
                       const Mesh* inputMesh,
                       const AppOptions* options)
{
    LOG_TITLE("Spine Segmentation");

    auto proxySpineMorphologies = skeletonizer->reconstructSpineProxyMorphologies();

    auto spineMeshes = skeletonizer->reconstructSpineMeshes(inputMesh,
                                                            options->spinesVoxelsPerMicron, 0.2);

    std::vector<Mesh*> remeshedSpines;
    remeshedSpines.resize(spineMeshes.size());

    TIMER_SET;
    LOOP_STARTS("Spine Remeshing");
    PROGRESS_SET;
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < spineMeshes.size(); ++i)
    {
        auto remeshedSpine = remeshSpine(spineMeshes[i], options->spinesVoxelsPerMicron, VERBOSE);
        remeshedSpines[i] = remeshedSpine;

        LOOP_PROGRESS(PROGRESS, spineMeshes.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Exporting the spines
    TIMER_RESET;
    LOOP_STARTS("Exporting Segmented Spines");
    for (size_t i = 0; i < spineMeshes.size(); ++i)
    {
        auto spineMesh = spineMeshes[i];
        std::stringstream stream;
        stream << options->spinesMeshPrefix << SPINE_MESH_SUFFIX << "-" << i;
        spineMesh->exportMesh(stream.str(),
                              options->exportOBJ, options->exportPLY,
                              options->exportOFF, options->exportSTL, SILENT);

        auto& remeshedSpine = remeshedSpines[i];
        if (remeshedSpine == nullptr) continue;

//        stream << REFINED_SPINE_MESH_SUFFIX;
//        remeshedSpine->exportMesh(stream.str(),
//                                  options->exportOBJ, options->exportPLY,
//                                  options->exportOFF, options->exportSTL, SILENT);

        remeshedSpine->translate(-proxySpineMorphologies[i]->getBasePoint());

        remeshedSpine->rotateTowardsTargetPoint(proxySpineMorphologies[i]->getBasePoint(),
                    proxySpineMorphologies[i]->getDirection().normalized(),
                                                proxySpineMorphologies[i]->getBasePoint() + 5 * Vector3f(0, 1, 0));
        stream << REFINED_SPINE_MESH_SUFFIX << "_rotatedy";
        remeshedSpine->exportMesh(stream.str(),
                                  options->exportOBJ, options->exportPLY,
                                  options->exportOFF, options->exportSTL, SILENT);

        LOOP_PROGRESS(i, spineMeshes.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clear the spines
    for (auto& mesh : spineMeshes) { if (mesh != nullptr) mesh->~Mesh();  }

    // Clear the re-meshed spines
    for (auto& mesh : remeshedSpines) { if (mesh != nullptr) mesh->~Mesh(); }
}

void runSpineSegmentationOperations(const AppOptions* options,
                                    NeuronSkeletonizer* skeletonizer,
                                    Mesh* neuronMesh)
{
    // Export the spine meshes
    if (options->exportSpineMeshes)
    {
        if (!options->removeSpinesFromSkeleton)
        {
            LOG_WARNING("To export the spine meshes, add the flag --remove-spines-from-skeleton "
                        "to segment the spines!");
        }
        else
        {
            exportSpineMeshes(skeletonizer, neuronMesh, options);
            skeletonizer->exportSpineExtents(options->debuggingPrefix);
        }
    }
}


}
