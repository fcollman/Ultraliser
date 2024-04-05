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
#include <algorithms/skeletonization/SomaSegmenter.h>

namespace Ultraliser
{

void runFixSomaSlicingArtifactsOperations(const AppOptions* options, Mesh* neuronMesh)
{
    if (options->fixSomaSlicingArtifacts)
    {
        LOG_TITLE("Fixing Soma Slicing Artifacts");
        TIMER_SET;

        // Use the SomaSegmenter to segment a coarse soma mesh from the neuron mesh if the mesh
        // has obvious slicing artifacts
        LOG_STATUS("Coarse Skeletonization to Segment Initial Soma Profile");
        std::unique_ptr< SomaSegmenter > somaSegmenter =
                std::make_unique< SomaSegmenter >(neuronMesh, options->somaSegmenterVPM);
        auto somaMesh = somaSegmenter->segmentSomaMesh(SILENT);
        LOG_STATS(GET_TIME_SECONDS);

        // Construct a point cloud for the neuron mesh to map the reconstructed soma to the neuron mesh
        LOG_STATUS("Mapping Initial Soma Profile to Neuron Mesh");
        std::vector< Vector3f > neuronMeshCloud;
        neuronMeshCloud.resize(neuronMesh->getNumberVertices());

        OMP_PARALLEL_FOR
                for (size_t i = 0; i < neuronMesh->getNumberVertices(); ++i)
        {
            neuronMeshCloud[i] = neuronMesh->_vertices[i];
        }

        // Map the source mesh to the destination mesh
        somaMesh->kdTreeMapping(neuronMeshCloud, SILENT);
        LOG_STATS(GET_TIME_SECONDS);

        // Append this soma mesh to the given neuron mesh to fill the slices in the soma
        LOG_STATUS("Updating Neuron Mesh with Corrected Soma");
        neuronMesh->append(somaMesh);
        LOG_STATS(GET_TIME_SECONDS);

        LOG_STATUS_IMPORTANT("Fixing Soma Slicing Artifacts Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOG_WARNING("The soma is assumed to have NO missing slices.\n"
                    "\tIf you wish to fix any slicing artifacts in the soma, or in case the\n"
                    "\tapplication fails to continue due to gaps or missing slices within the soma\n"
                    "\textent, please use the --fix-soma-slicing-artifacts command to fix the soma!");
    }
}

void runSomaExportOperations(const AppOptions* options, NeuronSkeletonizer* skeletonizer)
{
    // Export the somatic proxy mesh
    if (options->exportProxySomaMesh)
    {
        skeletonizer->exportSomaProxyMesh(options->meshPrefix,
            options->exportOBJ, options->exportPLY, options->exportOFF, options->exportSTL);
    }

    // Export the somatic mesh
    if (options->exportSomaMesh)
    {
        skeletonizer->exportSomaMesh(options->meshPrefix,
            options->exportOBJ, options->exportPLY, options->exportOFF, options->exportSTL);
    }
}

}
