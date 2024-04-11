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
#include <algorithms/skeletonization/SpineSkeletonizer.h>

// Defines
#define NEURON_SMOOTHING_ITERATIONS 25

namespace Ultraliser
{

Mesh* remeshSpine(Mesh* spineMesh, const float voxelsPerMicron = 50,
                  const bool verbose = SILENT)
{
    // Compute the bounding box of the input neuron mesh
    Vector3f pMinInput, pMaxInput;
    spineMesh->computeBoundingBox(pMinInput, pMaxInput);
    const auto& bounds = pMaxInput - pMinInput;

    if (bounds.x() < 1e-16 || bounds.y() < 1e-16 || bounds.z() < 1e-16)
    {
        return nullptr;
    }

    // Compute the resolution of the volume
    const auto largestDimension = bounds.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, 0.1,
                                VOLUME_TYPE::UI8, SILENT);

    // Apply surface and solid voxelization to the input neuron mesh
    volume->surfaceVoxelization(spineMesh, SILENT, false, 1.0);
    volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ, SILENT);
    // volume->surfaceVoxelization(inputSpineMesh, SILENT, false, 0.5);

    // Construct the mesh using the DMC technique
    auto reconstructedSpineMesh = DualMarchingCubes::generateMeshFromVolume(volume, SILENT);

    // Smooth the resulting surface mesh
    reconstructedSpineMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, SILENT);

    // Return a pointer to the resulting neuron
    return reconstructedSpineMesh;
}

Meshes remeshProxySpineMeshes(const Meshes& proxySpineMeshes,
                              const AppOptions* options)
{
    // Construct a list to contain the final results of the spine meshes
    std::vector< Mesh* > spineMeshes;
    spineMeshes.resize(proxySpineMeshes.size());

    TIMER_SET;
    LOOP_STARTS("Remeshing Proxy Spine Meshes");
    for (size_t i = 0; i < proxySpineMeshes.size(); ++i)
    {
        if (proxySpineMeshes[i] == nullptr)  { spineMeshes[i] = nullptr; }
        else
        {
            spineMeshes[i] = remeshSpine(proxySpineMeshes[i],
                                         options->spinesVoxelsPerMicron, VERBOSE);
            LOOP_PROGRESS(i, proxySpineMeshes.size());
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return spineMeshes;
}

void exportDendriticProxySpineMeshes(const Meshes& proxyDendriticSpineMeshes,
                                     const AppOptions* options)
{
    TIMER_SET;
    LOOP_STARTS("Exporting Segmented Dendritic (Proxy) Spines");
    for (size_t i = 0; i < proxyDendriticSpineMeshes.size(); ++i)
    {
        const auto proxySpineMesh = proxyDendriticSpineMeshes[i];
        if (proxySpineMesh == nullptr) continue;

        std::stringstream stream;
        stream << options->dendriticSpinesProxyMeshesPrefix
               << DENDRITIC_SPINE_PROXY_SUFFIX << "-" << i;
        proxySpineMesh->exportMesh(stream.str(),
                                   options->exportOBJ, options->exportPLY,
                                   options->exportOFF, options->exportSTL, SILENT);

        LOOP_PROGRESS(i, proxyDendriticSpineMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void exportDendriticSpineMeshes(const Meshes& dendriticSpineMeshes, const AppOptions* options)
{
    TIMER_SET;
    LOOP_STARTS("Exporting Segmented Dendritic (Remeshed) Spines");
    for (size_t i = 0; i < dendriticSpineMeshes.size(); ++i)
    {
        const auto proxySpineMesh = dendriticSpineMeshes[i];
        if (proxySpineMesh == nullptr) continue;

        std::stringstream stream;
        stream << options->dendriticSpinesMeshesPrefix
               << DENDRITIC_SPINE_SUFFIX << "-" << i;
        proxySpineMesh->exportMesh(stream.str(),
                                   options->exportOBJ, options->exportPLY,
                                   options->exportOFF, options->exportSTL, SILENT);

        LOOP_PROGRESS(i, dendriticSpineMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void exportSpineMeshes(const Meshes& spinesMeshes, const AppOptions* options)
{
    TIMER_SET;
    LOOP_STARTS("Exporting Segmented Dendritic (Remeshed) Spines");
    for (size_t i = 0; i < spinesMeshes.size(); ++i)
    {
        const auto spineMesh = spinesMeshes[i];
        if (spineMesh == nullptr) continue;

        std::stringstream stream;
        stream << options->spinesMeshesPrefix << SPINE_SUFFIX << "-" << i;
        spineMesh->exportMesh(stream.str(),
                              options->exportOBJ, options->exportPLY,
                              options->exportOFF, options->exportSTL, SILENT);

        LOOP_PROGRESS(i, spinesMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void transformSpinesToOrigin(Meshes& dendriticSpineMeshes,
                             SpineMorphologies& proxySpineMorphologies)
{
    TIMER_SET;
    LOOP_STARTS("Transforming Spines to Origin");
    for (size_t i = 0; i < dendriticSpineMeshes.size(); ++i)
    {
        const auto spineMesh = dendriticSpineMeshes[i];
        if (spineMesh == nullptr) continue;

        // Get the base point to translate the spine to the origin
        const auto basePoint = proxySpineMorphologies[i]->getBasePoint();

        // Get thr orientation data needed to make the spine point up towards the Y-axis
        const auto spineDirection = proxySpineMorphologies[i]->getDirection().normalized();
        const auto targetPoint = basePoint + 5 * Vector3f(0, 1, 0);

        // Rotate the spine towards the Y-axis
        spineMesh->rotateTowardsTargetPoint(basePoint, spineDirection, targetPoint);

        // Translate the spine mesh
       spineMesh->translate(-basePoint);

        LOOP_PROGRESS(i, dendriticSpineMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

bool reconstructAndExportSpineMorphology(Mesh* spineMesh,
                                         const size_t spineIndex,
                                         const Vector3f& spineBasePoint,
                                         const AppOptions* options)
{
    // TODO: Construct the voxelization options of the spine
    Skeletonizer::VoxelizationOptions spineVoxelizationOptions;
    spineVoxelizationOptions.volumeResolution = 256;
    spineVoxelizationOptions.verbose = SILENT;
    spineVoxelizationOptions.edgeGapPrecentage = 0.5;
    std::unique_ptr< SpineSkeletonizer > skeletonizer = std::make_unique< SpineSkeletonizer >(
                spineMesh, spineIndex, spineBasePoint, spineVoxelizationOptions, true, false);

    // Run the spine skeletonization operation
    const bool success = skeletonizer->runSkeletonization(SILENT);
    if (success && options->exportSpinesSWCMorphologies)
    {
        // Export the spine morphology
        const bool resampleSkeleton = false;
        std::stringstream prefix;
        prefix << options->spinesMorphologiesPrefix << SPINE_SUFFIX << "-" << spineIndex;
        skeletonizer->exportSWCFile(prefix.str(), resampleSkeleton, SILENT);

        // Spine has been exported
        return true;
    }

    // Spine has an issue, report it
    return false;
}

void reconstructSpineMorphologies(const Meshes& spineMeshes,
                                  const SpineMorphologies& spineProxyMorphologies,
                                  const AppOptions* options)
{
    TIMER_SET;
    LOOP_STARTS("Reconstructing HQ Spine Morphologies");
    for (size_t i = 0; i < spineMeshes.size(); ++i)
    {
        auto spineMesh = spineMeshes[i];
        if (spineMesh == nullptr) continue;

        // The base point at this stage should be the origin
        // const auto basePoint = spineProxyMorphologies[i]->getBasePoint();
        const Vector3f basePoint = Vector3f(0.f);

        // Reconstruct the spine morphology
        reconstructAndExportSpineMorphology(spineMesh, i, Vector3f(0.), options);

        LOOP_PROGRESS(i, spineMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void reconstructSpines(NeuronSkeletonizer* skeletonizer,
                       const Mesh* inputNeuronMesh,
                       const AppOptions* options)
{
    // Reconstruct the dendritic proxy morphologies of the spines (not final and located on spines)
    auto proxySpinesMorphologies = skeletonizer->reconstructSpineProxyMorphologies();

    // Reconstruct dendritic proxy spine meshes (not final and located on spines)
    auto proxySpinesMeshes = skeletonizer->reconstructSpineMeshes(
                inputNeuronMesh, options->spinesVoxelsPerMicron, 0.1);

    // Export the dendritic proxy spine meshes, if needed
    if (options->exportDenditicSpinesProxyMeshes)
    {
        exportDendriticProxySpineMeshes(proxySpinesMeshes, options) ;
    }

    // Reconstruct the (final) dendritic spine meshes by remeshing the proxy spine meshes
    auto spinesMeshes = remeshProxySpineMeshes(proxySpinesMeshes, options);

    // Export the dendritic spines meshes, if needed
    if (options->exportDendriticSpinesMeshes)
    {
        exportDendriticSpineMeshes(spinesMeshes, options);
    }

    // Transform the spines to the origin
    transformSpinesToOrigin(spinesMeshes, proxySpinesMorphologies);

    // Export the spines meshes at the origin, if needed
    if (options->exportSpinesMeshes)
    {
        exportSpineMeshes(spinesMeshes, options);
    }

    // Clear the proxy spine meshes
    for (auto& mesh : proxySpinesMeshes) { if (mesh != nullptr) mesh->~Mesh();  }

    // Reconstruct the spine morphologies
    if (options->exportSpinesSWCMorphologies)
    {
        reconstructSpineMorphologies(spinesMeshes, proxySpinesMorphologies, options);
    }

    // Clear the spine proxy morphologies
    for (auto& morphology : proxySpinesMorphologies)
    {
        if (morphology != nullptr) morphology->~SpineMorphology();
    }

    // Clear the spine meshes
    for (auto& mesh : spinesMeshes) { if (mesh != nullptr) mesh->~Mesh();  }

}

void exportSpineMeshes(NeuronSkeletonizer* skeletonizer,
                       const Mesh* inputNeuronMesh,
                       const AppOptions* options)
{
    LOG_TITLE("Spine Segmentation");

    // Reconstruct the proxy morphologies of the spines (not final)
    auto proxySpineMorphologies = skeletonizer->reconstructSpineProxyMorphologies();

    // Reconstruct proxy spine meshes (not final)
    auto proxySpineMeshes = skeletonizer->reconstructSpineMeshes(
                inputNeuronMesh, options->spinesVoxelsPerMicron, 0.2);

    std::vector< Mesh* > remeshedSpines;
    remeshedSpines.resize(proxySpineMeshes.size());

    TIMER_SET;
    LOOP_STARTS("Spine Remeshing");
    for (size_t i = 0; i < proxySpineMeshes.size(); ++i)
    {
        if (proxySpineMeshes[i] == nullptr)
        {
            remeshedSpines[i] = nullptr;
            continue;
        }

        remeshedSpines[i] = remeshSpine(proxySpineMeshes[i], options->spinesVoxelsPerMicron, VERBOSE);

        LOOP_PROGRESS(i, proxySpineMeshes.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);


    // Clear the spines
    for (auto& mesh : proxySpineMeshes) { if (mesh != nullptr) mesh->~Mesh();  }

    // Clear the re-meshed spines
    for (auto& mesh : remeshedSpines) { if (mesh != nullptr) mesh->~Mesh(); }
}

void runSpineSegmentationOperations(const AppOptions* options,
                                    NeuronSkeletonizer* skeletonizer,
                                    Mesh* neuronMesh)
{
    // Ensure the presence of the spine removal from skeleton
    if (options->exportDenditicSpinesProxyMeshes    ||
        options->exportDendriticSpinesMeshes        ||
        options->exportSpinesMeshes                 ||
        options->exportSpinesSWCMorphologies)
    {
        if (!options->removeSpinesFromSkeleton)
        {
            LOG_WARNING("To perform spine segmentation and export operations, add the flag "
                        "--remove-spines-from-skeleton to segment the spines at first!");
            return;
        }

        // Reconstruct the spines
        reconstructSpines(skeletonizer, neuronMesh, options);
    }


//    // Export the spine meshes
//    if (options->exportSpineMeshes)
//    {
//        if (!options->removeSpinesFromSkeleton)
//        {
//            LOG_WARNING("To export the spine meshes, add the flag --remove-spines-from-skeleton "
//                        "to segment the spines!");
//        }
//        else
//        {
//            exportSpineMeshes(skeletonizer, neuronMesh, options);
//            skeletonizer->exportSpineExtents(options->debugPrefix);
//        }
//    }
}

}
