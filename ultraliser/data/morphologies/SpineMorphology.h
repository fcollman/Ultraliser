/***************************************************************************************************
 * Copyright (c) 2016 - 2024
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include <data/morphologies/Morphology.h>
#include <data/meshes/Meshes.h>
#include <algorithms/skeletonization/SkeletonBranch.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/// Forward declaration
class Volume;

/**
 * @brief The SpineMorphology class
 * Spine morphology
 */
class SpineMorphology: public Morphology
{
public:

    /**
     * @brief SpineMorphology
     * Construct the spine morphology from a tree of branches, where the input branch is the root
     * of the tree.
     * @param root
     * The root of the tree.
     * @param spineIndex
     * The index of the spine.
     */
    SpineMorphology(SkeletonBranch* root,
                    const size_t& spineIndex);

    /**
     * @brief SpineMorphology
     * @param branches
     * @note It's recommended to avoid using this constructor and use the other one to construct
     * the tree.
     */
    SpineMorphology(SkeletonBranches branches,
                    const size_t& index);

public:

    /**
     * @brief getBasePoint
     * Returns the base point of the spine.
     * @return
     */
    Vector3f getBasePoint() const { return _basePoint; }

    /**
     * @brief getDirection
     * @return
     */
    Vector3f getDirection() const {return _direction; }

    /**
     * @brief getSpineIndex
     * @return
     */
    size_t getSpineIndex() const { return _spineIndex; }

    /**
     * @brief generateMesh
     * @param voxelsPerMicron
     * @param edgeGap
     * @param verbose
     * @return
     */
    Mesh* generateMesh(const float& voxelsPerMicron,
                       const float& edgeGap = 0.1,
                       const bool& verbose = SILENT);
    /**
     * @brief generateVolume
     * @param voxelsPerMicron
     * @param verbose
     * @return
     */
    Volume* generateVolume(const float& voxelsPerMicron,
                           const float& edgeGap = 0.1,
                           const bool& verbose = SILENT);

    /**
     * @brief generateMeshWithoutDenderiticExtent
     * @param voxelsPerMicron
     * @param edgeGap
     * @param verbose
     * @return
     */
    Mesh* generateMeshWithoutDenderiticExtent(const float& voxelsPerMicron,
                                              const float& edgeGap,
                                              const bool & verbose = SILENT);

    /**
     * @brief generateVolumeWithoutDenderiticExtent
     * @param voxelsPerMicron
     * @param edgeGap
     * @param verbose
     * @return
     */
    Volume* generateVolumeWithoutDenderiticExtent(const float& voxelsPerMicron,
                                                  const float& edgeGap,
                                                  const bool & verbose);

    /**
     * @brief exportBranches
     * @param prefix
     * @param verbose
     */
    void exportBranches(const std::string &prefix,
                        const bool &verbose = SILENT);

    /**
     * @brief exportExtents
     * @param prefix
     * @param showProgrress
     */
    void exportExtents(const std::string& prefix) const;

    /**
     * @brief getTerminals
     * @return
     */
    std::vector< Vector3f > getTerminals() const { return _terminals; }

private:

    /**
     * @brief _getNonDendrticSections
     * @return
     */
    Sections _getNonDendrticSections() const;

    /**
     * @brief _constructTreeFromLogicalBranches
     * @param root
     * @param sectionIndex
     */
    void _constructTreeFromLogicalBranches(SkeletonBranch* root,
                                           size_t& sectionIndex);

    /**
     * @brief _computeBoundingBox
     */
    void _computeBoundingBox();

    /**
     * @brief _computeOrientationVector
     */
    void _computeOrientationVector();

private:

    /**
     * @brief spineIndex
     */
    size_t _spineIndex;

    /**
     * @brief _radfiusScaleFactor
     */
    float _radfiusScaleFactor = 1.0;

    /**
     * @brief _rootSample
     */
    Sample* _rootSample;

    /**
     * @brief _basePoint
     */
    Vector3f _basePoint;

    /**
     * @brief _direction
     */
    Vector3f _direction;

    /**
     * @brief _terminals
     */
    std::vector< Vector3f > _terminals;
};

/**
 * @brief SpineMorphologies
 */
typedef std::vector< SpineMorphology* > SpineMorphologies;

}
