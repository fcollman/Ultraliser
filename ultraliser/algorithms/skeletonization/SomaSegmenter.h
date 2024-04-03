/***************************************************************************************************
 * Copyright (c) 2016 - 2024
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <algorithms/skeletonization/NeuronSkeletonizer.h>
#include <data/meshes/simple/Mesh.h>
#include <data/volumes/Volume.h>

namespace Ultraliser
{

/**
 * @brief The SomaSegmenter class
 */
class SomaSegmenter
{
public:

    /**
     * @brief SomaSegmenter
     * @param neuronMesh
     * @param voxelizationVoxelsPerMicron
     */
    SomaSegmenter(Mesh *neuronMesh, const float voxelizationVoxelsPerMicron = 1);

    /**
     * @brief SomaSegmenter
     * @param neuronVolume
     */
    SomaSegmenter(Volume* neuronVolume);
    ~SomaSegmenter();

public:

    /**
     * @brief segmentSomaProxyMesh
     * @param verbose
     * @return
     */
    Mesh* segmentSomaProxyMesh(const bool verbose = SILENT);

    /**
     * @brief segmentSomaMesh
     * @param verbose
     * @return
     */
    Mesh* segmentSomaMesh(const bool verbose = SILENT);

private:

    /**
     * @brief _verifyVolume
     * @param verbose
     */
    void _verifyVolume(const bool verbose = SILENT);

private:

    /**
     * @brief _neuronMesh
     */
    Mesh* _neuronMesh;

    /**
     * @brief _neuronVolume
     */
    Volume* _neuronVolume = nullptr;

    /**
     * @brief _voxelizationVoxelsPerMicron
     */
    float _voxelizationVoxelsPerMicron;
};

}
