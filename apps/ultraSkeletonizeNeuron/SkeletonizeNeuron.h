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

#pragma once

#include <AppOptions.h>
#include <algorithms/skeletonization/NeuronSkeletonizer.h>

namespace Ultraliser
{

/**
 * @brief createNeuronVolume
 * @param neuronMesh
 * @param options
 * @param verbose
 * @return
 */
Volume* createNeuronVolume(Mesh* neuronMesh,
                           const AppOptions* options,
                           const bool verbose = VERBOSE);

/**
 * @brief reconstructNeuronMeshFromVolume
 * @param neuronVolume
 * @param options
 * @param verbose
 * @return
 */
Mesh* reconstructNeuronMeshFromVolume(Volume* neuronVolume,
                                      AppOptions* options,
                                      const bool verbose = VERBOSE);

/**
 * @brief remeshNeuron
 * @param inputNeuronMesh
 * @param options
 * @param verbose
 * @return
 */
Mesh* remeshNeuron(Mesh* inputNeuronMesh,
                   AppOptions* options,
                   const bool verbose = VERBOSE);

/**
 * @brief projectInputNeuronVolume
 * @param options
 * @param neuronVolume
 */
void projectInputNeuronVolume(const AppOptions* options,
                              const Volume* neuronVolume);

/**
 * @brief createNeuronSkeletonizer
 * @param options
 * @param neuronVolume
 * @return
 */
NeuronSkeletonizer* createNeuronSkeletonizer(const AppOptions* options,
                                             Volume* neuronVolume);

/**
 * @brief runNeuronSkeletonizationOperations
 * @param options
 * @param skeletonizer
 */
void runNeuronSkeletonizationOperations(const AppOptions* options,
                                        NeuronSkeletonizer* skeletonizer);


/**
 * @brief runFixSomaSlicingArtifactsOperations
 * @param options
 * @param neuronMesh
 */
void runFixSomaSlicingArtifactsOperations(const AppOptions* options,
                                          Mesh* neuronMesh);

/**
 * @brief runSomaExportOperations
 * @param options
 * @param skeletonizer
 */
void runSomaExportOperations(const AppOptions* options,
                             NeuronSkeletonizer* skeletonizer);

/**
 * @brief runSpineSegmentationOperations
 * @param options
 * @param skeletonizer
 * @param neuronMesh
 */
void runSpineSegmentationOperations(const AppOptions* options,
                                    NeuronSkeletonizer* skeletonizer,
                                    Mesh* neuronMesh);

/**
 * @brief createMeshFromSections
 * @param sections
 * @param options
 * @return
 */
Mesh* createMeshFromSections(Sections& sections, const AppOptions *options);

/**
 * @brief createHighResolutionNeuronMesh
 * @param inputMesh
 * @param voxelsPerMicron
 * @return
 */
Mesh* createHighResolutionNeuronMesh(Mesh* inputMesh, float voxelsPerMicron);

/**
 * @brief runHighQualityMeshGeneration
 * @param options
 * @param skeletonizer
 * @param inputNeuronMesh
 */
void runHighQualityMeshGeneration(const AppOptions* options,
                                  NeuronSkeletonizer* skeletonizer,
                                  Mesh* inputNeuronMesh);

}
