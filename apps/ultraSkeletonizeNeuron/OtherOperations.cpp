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
#include <data/morphologies/Utilities.h>

// Defines
#define NEURON_SMOOTHING_ITERATIONS 25

namespace Ultraliser
{

Mesh* createMeshFromSections(Sections& sections, AppOptions* options)
{
    Vector3f pMinInput, pMaxInput;
    computeSectionsBoundingBox(sections, pMinInput, pMaxInput);
    const auto& bounds = pMaxInput - pMinInput;

    // Compute the resolution of the volume
    const auto largestDimension = bounds.getLargestDimension();
    size_t resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);

    // Construct the volume from the input mesh
    auto volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                             VolumeGrid::getType(options->volumeType), SILENT);

    volume->surfaceVoxelizeSections(sections);
    volume->solidVoxelization(options->voxelizationAxis);

    // Remove the border voxels that span less than half the voxel
    auto bordeVoxels = volume->searchForBorderVoxels(true);
    for (size_t i = 0; i < bordeVoxels.size(); ++i)
    {
        for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
        {
            auto voxel = bordeVoxels[i][j]; volume->clear(voxel.x(), voxel.y(), voxel.z());
        }
        bordeVoxels[i].clear();
    }
    bordeVoxels.clear();

    // Construct the mesh using the DMC technique
    auto resultMesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Smooth the resulting surface mesh
    resultMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, true);

    // Return a pointer to the resulting neuron
    return resultMesh;
}

}
