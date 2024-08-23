/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

#include <common/Defines.h>
#include "Skeletonizer.h"
#include "SkeletonizerUtils.h"
#include <algorithms/skeletonization/thinning/Neighbors.hh>
#include <algorithms/utilities/KdTree.h>
#include <math/Vector.h>
#include <data/meshes/simple/TriangleOperations.h>
#include <utilities/Range.h>
#include <data/volumes/voxels/NodeVoxel.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(Volume* volume,
                           const bool &useAcceleration,
                           const bool &debugSkeleton,
                           const std::string debugPrefix)
    : _volume(volume)
    , _mesh(nullptr)
    , _useAcceleration(useAcceleration)
    , _debugSkeleton(debugSkeleton)
    , _debugPrefix(debugPrefix)
    , _debug(_debugSkeleton && _debugPrefix != NONE)
{
    /// NOTE: The mesh is assigned a nullptr, until further notice

    // Mesh bounding box
    _pMinMesh = volume->getPMin();
    _pMaxMesh = volume->getPMax();
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    // TODO: Verify the volume point cloud
    // Volume bounding box
    _pMinVolume = Vector3f(0.f);
    _pMaxVolume = Vector3f((volume->getWidth() - 1) * 1.f,
                           (volume->getHeight() - 1) * 1.f,
                           (volume->getDepth() - 1) * 1.f);
    _boundsVolume = _pMaxVolume;
    _centerVolume = 0.5 * _boundsVolume;

    // Mesh to volume scale factor
    _scaleFactor = _boundsMesh / _boundsVolume;
}

Skeletonizer::Skeletonizer(Mesh* mesh,
                           const VoxelizationOptions& options,
                           const bool &useAcceleration,
                           const bool &debugSkeleton,
                           const std::string debugPrefix)
    : _volume(nullptr)
    , _mesh(mesh)
    , _voxelizationOptions(options)
    , _useAcceleration(useAcceleration)
    , _debugSkeleton(debugSkeleton)
    , _debugPrefix(debugPrefix)
    , _debug(_debugSkeleton && _debugPrefix != NONE)
{
    /// NOTE: The volume is assigned a nullptr, until further notice

    // Compute the mesh bounding box
    mesh->computeBoundingBox(_pMinMesh, _pMaxMesh);
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    /// NOTE: Compute the volume bounds after the generation of the volume
}

void Skeletonizer::_computeVolumeFromMesh()
{
    // If the mesh is a nullptr, then return, there is nothing to compute
    if (_mesh == nullptr)
    {
        LOG_ERROR("Skeletonizer::_computeVolumeFromMesh(): An empty mesh is given to voxelize!");
        return;
    }

    // The volume must be a nullptr to be able to compute it
    if (_volume == nullptr)
    {
        // Create the volume extent
        _volume = new Volume(_pMinMesh, _pMaxMesh,
                             256, // _voxelizationOptions.volumeResolution,
                             _voxelizationOptions.edgeGapPrecentage,
                             _voxelizationOptions.volumeType,
                             _voxelizationOptions.verbose);

        // Apply surface and solid voxelization to the input neuron mesh
        _volume->surfaceVoxelization(_mesh, _voxelizationOptions.verbose, false, 1.0);
        _volume->solidVoxelization(_voxelizationOptions.voxelizationAxis,
                                   _voxelizationOptions.verbose);

        // Remove the border voxels that span less than half the voxel
        auto bordeVoxels = _volume->searchForBorderVoxels(_voxelizationOptions.verbose);
        for (size_t i = 0; i < bordeVoxels.size(); ++i)
        {
            for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
            {
                _volume->clear(bordeVoxels[i][j].x(), bordeVoxels[i][j].y(), bordeVoxels[i][j].z());
            }
            bordeVoxels[i].clear();
        }
        bordeVoxels.clear();
        _volume->surfaceVoxelization(_mesh, _voxelizationOptions.verbose, false, 0.5);

        _pMinMesh = _volume->getPMin();
        _pMaxMesh = _volume->getPMax();
        _boundsMesh = _pMaxMesh - _pMinMesh;
        _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

        // Compute the scale factors
        // Volume bounding box
        _pMinVolume = Vector3f(0.f);
        _pMaxVolume = Vector3f((_volume->getWidth() - 1) * 1.f,
                               (_volume->getHeight() - 1) * 1.f,
                               (_volume->getDepth() - 1) * 1.f);

        _boundsVolume = _pMaxVolume;
        _centerVolume = 0.5 * _boundsVolume;

        // Mesh to volume scale factor
        _scaleFactor = _boundsMesh / _boundsVolume;
    }
}

void Skeletonizer::initialize(const bool verbose)
{
    TIMER_SET;
    VERBOSE_LOG(LOG_TITLE("Ultraliser Skeletonization"), verbose);
    VERBOSE_LOG(LOG_SUCCESS("Voxel Size [%f] μm", _volume->getVoxelSize()), verbose);
    VERBOSE_LOG(LOG_STATUS("Initialization - Building Structures"), verbose);

    // Compute the shell points either natively or by using the acceleration structures
    if (_useAcceleration)
    {
        // Build the ThinningVoxels acceleration structure from the input solid volume
        // NOTE: We do not rebuild the ThinningVoxels structure!
        auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

        // Compute the surface shell from the pre-built ThinningVoxels structure
        _computeShellPointsUsingAcceleration(thinningVoxels, verbose);
    }
    else { _computeShellPoints(verbose); }

    VERBOSE_LOG(LOG_STATUS_IMPORTANT("Initialization Stats."), verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

Sections Skeletonizer::getValidSections() const
{
    Sections sections;
    size_t sectionIndex = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isValid())
        {
            Section* section = new Section(sectionIndex++);
            for (size_t j = 0; j < _branches[i]->nodes.size(); ++j)
            {
                auto node = _branches[i]->nodes[j];
                if (node->isSoma) { continue; }
                section->addSample(new Sample(node->point, node->radius, j));
            }
            sections.push_back(section);
        }
    }

    return sections;
}

void Skeletonizer::_scaleShellPoints(const bool verbose)
{
    // Initialize the timer
    TIMER_SET;

    // TODO: Adjust the voxel slight shift
    // Adjust the locations of the shell points taking into consideration the mesh coordinates
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Mapping Shell Points"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _shellPoints.size(); ++i)
    {
        // Center the shell points (of the volume) at the origin
        _shellPoints[i] -= _centerVolume;

        // Scale to match the dimensions of the mesh
        _shellPoints[i].x() *= _scaleFactor.x();
        _shellPoints[i].y() *= _scaleFactor.y();
        _shellPoints[i].z() *= _scaleFactor.z();

        // Translate to the center of the mesh
        _shellPoints[i] += _centerMesh;

        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _shellPoints.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_computeShellPointsUsingAcceleration(ThinningVoxelsUI16List &thinningVoxels,
                                                        const bool verbose)
{
    // Initialize the timer
    TIMER_SET;

    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Computing Shell Points *"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        auto& voxel = thinningVoxels[i];
        if (_volume->isBorderVoxel(voxel->x, voxel->y, voxel->z))
        {
            voxel->border = true;
        }

        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, thinningVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Add all the obtained voxels in a single list
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        const auto& voxel = thinningVoxels[i];
        if (voxel->border)
        {
            _shellPoints.push_back(Vector3f(voxel->x, voxel->y, voxel->z));
        }
    }

    // Scale the shell points to match the extent of the input data
    _scaleShellPoints(verbose);
}

void Skeletonizer::_computeShellPoints(const bool verbose)
{
    // Initialize the time
    TIMER_SET;

    // Search for the border voxels (the shell voxels) of the volume
    std::vector< std::vector< Vec3ui_64 > > perSlice = _volume->searchForBorderVoxels(verbose);

    // Concatinate the points in a single list
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Computing Shell Points"), verbose);
    for (size_t i = 0; i < perSlice.size(); ++i)
    {
        for (size_t j = 0; j < perSlice[i].size(); ++j)
        {
            const auto voxel = perSlice[i][j];
            _shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSlice[i].clear();
        perSlice[i].shrink_to_fit();
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, perSlice.size()), verbose);
    }
    perSlice.clear();
    perSlice.shrink_to_fit();
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Scale the shell points to match the extent of the input data
    _scaleShellPoints(verbose);
}

void Skeletonizer::_applyVolumeThinning(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Volume Thinning"), verbose);

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Thinning Loop"), verbose);
    VERBOSE_LOG(LOOP_PROGRESS(0, 100), verbose);
    while(1)
    {
        size_t numberDeletedVoxels = _volume->deleteCandidateVoxelsParallel(thinningKernel);

        // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        VERBOSE_LOG(LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                      initialNumberVoxelsToBeDeleted), verbose);

        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_applyVolumeThinningUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Volume Thinning *"), verbose);

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Thinning Loop"), verbose);
    VERBOSE_LOG(LOOP_PROGRESS(0, 100), verbose);
    while(1)
    {
        // Delete the border voxels based on the ThinningVoxels acceleration structure
        size_t numberDeletedVoxels = _volume->deleteBorderVoxelsUsingThinningVoxels(
                    thinningKernel, thinningVoxels);

        // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        VERBOSE_LOG(LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                                  initialNumberVoxelsToBeDeleted), verbose);

        // No more voxels to be deleted
        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::skeletonizeVolumeToCenterLines(const bool verbose)
{
    if (_useAcceleration)
        _applyVolumeThinningUsingAcceleration(verbose);
    else
        _applyVolumeThinning(verbose);
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxels(const bool verbose)
{
    if (_useAcceleration)
        return _extractNodesFromVoxelsUsingAcceleration(verbose);
    else
        return _extractNodesFromVoxelsUsingSlicing(verbose);
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsNaive(const bool verbose)
{
    LOG_STATUS("Mapping Voxels to Nodes");

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;

    // A list of filled voxels to compute the elements in parallel
    std::vector< Vec4ui_64 > indicesFilledVoxels;

    // Search the filled voxels in the volume
    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Detecting Filled Voxels"), verbose);
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        for (size_t j = 0; j < _volume->getHeight(); ++j)
        {
            for (size_t k = 0; k < _volume->getDepth(); ++k)
            {
                // If the voxel is filled
                if (_volume->isFilled(i, j, k))
                {
                    // Get the 1D index of the voxel
                    size_t voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(i, j, k);

                    Vec4ui_64 index(i, j, k, voxelIndex);
                    indicesFilledVoxels.push_back(index);

                    // Mapper from voxel to node indices
                    indicesMapper.insert(std::pair< size_t, size_t >(voxelIndex, nodeIndex));

                    // New node
                    nodeIndex++;
                }
            }
        }
        VERBOSE_LOG(LOOP_PROGRESS(i, _volume->getWidth()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Resize the nodes
    _nodes.resize(indicesFilledVoxels.size());

    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Building Graph Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < indicesFilledVoxels.size(); ++n)
    {
        const size_t i = indicesFilledVoxels[n].x();
        const size_t j = indicesFilledVoxels[n].y();
        const size_t k = indicesFilledVoxels[n].z();
        const size_t voxelIndex = indicesFilledVoxels[n].w();

        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(i * 1.f, j * 1.f, k * 1.f);

        // Get a point in the same coordinate space of the mesh
        Vector3f nodePosition(voxelPosition);
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(voxelIndex, nodePosition, voxelPosition);

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, indicesFilledVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Clear the auxiliary lists
    indicesFilledVoxels.clear();
    indicesFilledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingSlicing(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Mapping Voxels to Nodes"), verbose);

    /**
     * @brief The FilledVoxel struct
     */
    struct FilledVoxel
    {
        size_t i, j, k, idx;

        FilledVoxel(size_t ii, size_t jj, size_t kk, size_t index)
        { i = ii; j = jj; k = kk; idx = index; }
    };

    /**
     * @brief FilledVoxels
     */
    typedef std::vector< FilledVoxel > FilledVoxels;

    // Make a per-slice list
    std::vector< FilledVoxels > allFilledVoxels;
    allFilledVoxels.resize(_volume->getWidth());

    TIMER_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        // Get a reference to the per-slice list
        auto& perSlice = allFilledVoxels[i];
        for (size_t j = 0; j < _volume->getHeight(); ++j)
        {
            for (size_t k = 0; k < _volume->getDepth(); ++k)
            {
                if (_volume->isFilled(i, j, k))
                {
                    perSlice.push_back(
                        FilledVoxel(i, j, k, _volume->mapTo1DIndexWithoutBoundCheck(i, j, k)));
                }
            }
        }
    }

    // Put them in a single list
    FilledVoxels filledVoxels;
    for (size_t i = 0; i < allFilledVoxels.size(); ++i)
    {
        if (allFilledVoxels[i].size() > 0)
        {
            filledVoxels.insert(filledVoxels.end(),
                                allFilledVoxels[i].begin(), allFilledVoxels[i].end());
            allFilledVoxels[i].clear();
            allFilledVoxels[i].shrink_to_fit();
        }
    }
    allFilledVoxels.clear();
    allFilledVoxels.shrink_to_fit();

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;
    for (size_t i = 0; i < filledVoxels.size(); ++i)
    {
        indicesMapper.insert(std::pair< size_t, size_t >(filledVoxels[i].idx, i));
    }

    // Resize the nodes
    _nodes.resize(filledVoxels.size());

    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < filledVoxels.size(); ++n)
    {
        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(filledVoxels[n].i * 1.f,
                               filledVoxels[n].j * 1.f,
                               filledVoxels[n].k * 1.f);

        // Get a point in the same coordinate space of the mesh
        Vector3f nodePosition(voxelPosition);
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(filledVoxels[n].idx, nodePosition, voxelPosition);

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, filledVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    filledVoxels.clear();
    filledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Mapping Voxels to Nodes *"), verbose);

    // Get all the center-line voxels from the volume
    auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

    // Construct the NodeVoxels list
    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Constructing Node Voxels"), verbose);
    PROGRESS_SET;
    NodeVoxelsUI16 nodeVoxels;
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        // Reference to the voxel
        auto& voxel = thinningVoxels[i];

        // The inactive voxels have been deactivated during the thinning
        if (voxel->active)
        {
            // Create a corresponding node to the voxel
            NodeVoxelUI16 nodeVoxel;

            // Get the location based on the voxel XYZ coordinates
            nodeVoxel.x = voxel->x; nodeVoxel.y = voxel->y; nodeVoxel.z = voxel->z;

            // Update the voxel index (used later to map the voxels to nodes)
            nodeVoxel.voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(
                        voxel->x, voxel->y, voxel->z);

            // Add th node voxel to the list
            nodeVoxels.push_back(nodeVoxel);
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, thinningVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;
    for (size_t i = 0; i < nodeVoxels.size(); ++i)
    {
        indicesMapper.insert(std::pair< size_t, size_t >(nodeVoxels[i].voxelIndex, i));
    }

    // Resize the nodes to the corresponding size of the NodeVoxels list
    _nodes.resize(nodeVoxels.size());

    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Building Graph Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < nodeVoxels.size(); ++n)
    {
        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(nodeVoxels[n].x * 1.f,
                               nodeVoxels[n].y * 1.f,
                               nodeVoxels[n].z * 1.f);

        // Get a point in the same coordinate space of the mesh
        /// TODO: Adjust the center of the node based on the actual center of the voxel
        Vector3f nodePosition(voxelPosition);

        // Adjust the location based on the dimensions of the input data
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(n, nodeVoxels[n].voxelIndex, nodePosition, voxelPosition);

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, nodeVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    return indicesMapper;
}

void Skeletonizer::_inflateNodes(const bool verbose)
{
    if (_useAcceleration)
        _inflateNodesUsingAcceleration(verbose);
    else
        _inflateNodesNatively(verbose);
}

void Skeletonizer::_inflateNodesUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Inflating Graph Nodes - Mapping to Surface"), verbose);

    auto kdtree = KdTree::from(_shellPoints);

    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Inflating Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto &node = *_nodes[i];

        auto nearestPoint = kdtree.findNearestPoint(node.point);
        auto minimumDistance = nearestPoint.distance;

        // TODO: Make some logic to detect the actual radius based on the voxel size
        if (minimumDistance > _volume->getVoxelSize())
        {
            node.radius = minimumDistance * 1.2;
        }
        else
        {
            node.radius = _volume->getVoxelSize() * 0.5;
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_inflateNodesNatively(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Inflating Graph Nodes - Mapping to Surface"), verbose);

    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Inflating Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        float minimumDistance = std::numeric_limits< float >::max();
        for (size_t j = 0; j < _shellPoints.size(); ++j)
        {
            const float distance = (_nodes[i]->point - _shellPoints[j]).abs();
            if (distance < minimumDistance) { minimumDistance = distance; }
        }

        // TODO: Make some logic to detect the actual radius based on the voxel size
        if (minimumDistance > 0.01)
        {
            _nodes[i]->radius = minimumDistance * 1.2;
        }
        else
        {
            _nodes[i]->radius = 0.1;
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_connectNodesToBuildEdges(const std::map< size_t, size_t >& indicesMapper,
                                             const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Connecting Graph Nodes"), verbose);

    size_t numberEdges;
    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Building & Linking Edges"), verbose);
    PROGRESS_SET;
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        // Check if the node has been visited before
        SkeletonNode* node = _nodes[i];

        // Count the number of the connected edges to the node
        size_t connectedEdges = 0;

        // Search for the neighbours
        for (size_t l = 0; l < 26; l++)
        {
            size_t idx = node->voxel.x() + VDX[l];
            size_t idy = node->voxel.y() + VDY[l];
            size_t idz = node->voxel.z() + VDZ[l];

            if (_volume->isFilled(idx, idy, idz))
            {
                // Increment the number of conected edges along this node
                connectedEdges++;

                // Find the index of the voxel
                const auto& vIndex = _volume->mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

                // Find the corresponding index of the node to access the node from the nodes list
                const auto& nIndex = indicesMapper.find(vIndex)->second;

                // Add the node to the edgeNodes, only to be able to access it later
                SkeletonNode* edgeNode = _nodes[nIndex];
                node->edgeNodes.push_back(edgeNode);

                // Construct the edge
                SkeletonEdge* edge = new SkeletonEdge(numberEdges, node, edgeNode);

                // Add the constructed edge to the list
                _edges.push_back(edge);

                // Increment the number of edges
                numberEdges++;
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;

        if (connectedEdges > 2)
            node->branching = true;

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_removeTriangleLoops(const bool verbose)
{
    TIMER_SET;
    VERBOSE_LOG(LOG_STATUS("Removing Triangle Loops"), verbose);

    const size_t currentNodesSize = _nodes.size();
    for (size_t i = 0; i < currentNodesSize; ++i)
    {
        if (_nodes[i]->branching)
        {
            SkeletonNodes sideNodes;
            if (isTriangleNode(_nodes[i], sideNodes))
            {
                if (_nodes[i]->visited) continue;

                auto& n1 = _nodes[i];
                auto& n2 = sideNodes[0];
                auto& n3 = sideNodes[1];

                // Collapse a triangle into a single node
                collapseTriangleIntoNode(_nodes, n1, n2, n3);

                if (n1->edgeNodes.size() > 2)
                    n1->branching = true;
                else
                    n1->branching = false;

                if (n2->edgeNodes.size() > 2)
                    n2->branching = true;
                else
                    n2->branching = false;

                if (n3->edgeNodes.size() > 2)
                    n3->branching = true;
                else
                    n3->branching = false;

                n1->visited = true;
                n2->visited = true;
                n3->visited = true;
            }
        }

        VERBOSE_LOG(LOOP_PROGRESS(i, currentNodesSize), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Reset the visiting state
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i) { _nodes[i]->visited = false; }
}

void Skeletonizer::constructGraph(const bool verbose)
{
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxels();

    // Assign accurate radii to the nodes of the graph, i.e. inflate the nodes
    _inflateNodes(verbose);

    // Connect the nodes to construct the edges of the graph
    _connectNodesToBuildEdges(indicesMapper);

    // Remove the triangular configurations
    _removeTriangleLoops(verbose);

    // Re-index the samples, for simplicity
    OMP_PARALLEL_FOR for (size_t i = 1; i <= _nodes.size(); ++i) { _nodes[i - 1]->index = i; }
}

void Skeletonizer::_buildBranchesFromNodes(const SkeletonNodes& nodes)
{
    // Used to index the branch
    size_t branchIndex = 0;

    // Construct the hierarchy to the terminal
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto& node = nodes[i];

        // The node must be branching
        if (node->branching)
        {
            // The node must be visited less number of times than its branching edges
            if (node->iVisit < node->edgeNodes.size())
            {
                // Construct the branch, starting with the edge node
                for (size_t j = 0; j < node->edgeNodes.size(); ++j)
                {
                    // Get a reference to the edge node
                    auto& edgeNode = node->edgeNodes[j];

                    if (edgeNode->iVisit >= edgeNode->edgeNodes.size()) continue;

                    // If the edge node is a terminal
                    if (edgeNode->terminal)
                    {
                        SkeletonBranch* branch = new SkeletonBranch();

                        node->iVisit += 1;
                        branch->nodes.push_back(node);

                        edgeNode->iVisit += 1;
                        branch->nodes.push_back(edgeNode);

                        branch->index = branchIndex;
                        branchIndex++;
                        _branches.push_back(branch);
                    }

                    // If the edge node is a branching node
                    else if (edgeNode->branching)
                    {
                        SkeletonBranch* branch = new SkeletonBranch();

                        node->iVisit += 1;
                        branch->nodes.push_back(node);

                        edgeNode->iVisit += 1;
                        branch->nodes.push_back(edgeNode);

                        branch->index = branchIndex;
                        branchIndex++;
                        _branches.push_back(branch);
                    }

                    // If the edge node is an intermediate node
                    else
                    {
                        // Ensure that the edge node is not visited before to make a branch
                        if (edgeNode->iVisit < 1)
                        {
                            SkeletonBranch* branch = new SkeletonBranch();

                            node->iVisit += 1;
                            branch->nodes.push_back(node);

                            edgeNode->iVisit += 1;
                            branch->nodes.push_back(edgeNode);

                            // The previous node is the first node
                            SkeletonNode *previousNode = node;

                            // The current node is the edge node
                            SkeletonNode *currentNode = edgeNode;

                            // Ensure that the current node has only two connected edges (or nodes)
                            while (true)
                            {
                                // Get a reference to the connecting nodes to the current node
                                auto edgeNode0 = currentNode->edgeNodes[0];
                                auto edgeNode1 = currentNode->edgeNodes[1];

                                // Ignore the previous node
                                if (edgeNode0->index == previousNode->index)
                                {
                                    previousNode = currentNode;
                                    currentNode = edgeNode1;
                                }
                                else
                                {
                                    previousNode = currentNode;
                                    currentNode = edgeNode0;
                                }

                                currentNode->iVisit += 1;
                                branch->nodes.push_back(currentNode);

                                if (!(currentNode->edgeNodes.size() == 2))
                                    break;
                            }

                            branch->index = branchIndex;

                            branchIndex++;
                            _branches.push_back(branch);
                        }
                    }
                }
            }
        }
    }
}

void Skeletonizer::skeletonizeVolumeBlockByBlock(const size_t& blockSize,
                                                 const size_t& numberOverlappingVoxels,
                                                 const size_t& numberZeroVoxels,
                                                 const bool verbose)
{
    thinVolumeBlockByBlock(blockSize, numberOverlappingVoxels, numberZeroVoxels);
    constructGraph(verbose);
    segmentComponents();
}

void Skeletonizer::_findClosestNodesInTwoPartitions(GraphComponent& partition1,
                                                    GraphComponent& partition2,
                                                    size_t* partition1NodeIndex,
                                                    size_t* partition2NodeIndex,
                                                    float* distance)
{
    // A parameter to store the shotrest distances between the nodes
    float shortestDistance = std::numeric_limits<float>::max();
    size_t indexNode1, indexNode2;

    // Iterate over the nodes of the first partition
    for (size_t i = 0; i < partition1.size(); ++i)
    {
        // Reference to the first node
        auto _node1 = _nodes[partition1[i]];

        // Iterate over the nodes of the second partition
        for (size_t j = 0; j < partition2.size(); ++j)
        {
            // Reference to the second node
            auto _node2 = _nodes[partition2[j]];

            // Calculate the distance
            const auto distance = _node1->point.distance(_node2->point);

            // If the calculated distance is less than the shortest distance, update
            if (distance < shortestDistance)
            {
                shortestDistance = distance;
                indexNode1 = i;
                indexNode2 = j;
            }
        }
    }

    *distance = shortestDistance;
    *partition1NodeIndex = indexNode1;
    *partition2NodeIndex = indexNode2;
}

void Skeletonizer::_connectPartition(GraphComponents& partitions,
                                     const size_t& partitionIndex,
                                     SkeletonEdges &edges)
{
    // A reference to the primary partition
    auto primaryPartition = partitions[partitionIndex];

    size_t partition1NodeIndex;
    size_t partition2NodeIndex;
    float shortestDistance = std::numeric_limits<float>::max();
    size_t secondaryPartitionIndex;

    for (size_t i = 0; i < partitions.size(); ++i)
    {
        // Skip the same partition
        if (i == partitionIndex) continue;

        // A reference to the secondary partition
        auto secondaryPartition = partitions[i];

        size_t _node1Index, _node2Index;
        float distance;

        _findClosestNodesInTwoPartitions(primaryPartition, secondaryPartition,
                                         &_node1Index, &_node2Index, &distance);

        if (distance < shortestDistance)
        {
            shortestDistance = distance;
            partition1NodeIndex = _node1Index;
            partition2NodeIndex = _node2Index;
            secondaryPartitionIndex = i;
        }
    }

    // Add the missing connectivity information
    auto primaryNode = _nodes[primaryPartition[partition1NodeIndex]];
    auto secondaryNode = _nodes[partitions[secondaryPartitionIndex][partition2NodeIndex]];

    primaryNode->edgeNodes.push_back(secondaryNode);
    secondaryNode->edgeNodes.push_back(primaryNode);

    SkeletonEdge* edge = new SkeletonEdge(edges.size(), primaryNode, secondaryNode);
    edges.push_back(edge);
}

void Skeletonizer::_verifyGraphConnectivityToClosestPartition(SkeletonEdges &edges, const bool verbose)
{
    while (true)
    {
        // Construct the graph
        auto graph = new Graph(edges, _nodes);

        // Get the number of partitions
        auto components = graph->getComponents();

        if (components.size() == 1)
        {
            VERBOSE_LOG(LOG_SUCCESS("The skeleton graph has 1 component! OK."), verbose);
            return;
        }
        else
        {
            VERBOSE_LOG(LOG_WARNING("The skeleton graph has [ %d ] components! "
                        "Running the Connectomics Algorithm", components.size()), verbose);

            // Add the partition
            _connectPartition(components, 0, edges);

            // Updating the branching and terminal nodes
            for (size_t i = 0; i < _nodes.size(); ++i)
            {
                // Check if the node has been visited before
                SkeletonNode* node = _nodes[i];

                if (node->edgeNodes.size() == 1)
                    node->terminal = true;
                else
                    node->terminal = false;

                if (node->edgeNodes.size() > 2)
                    node->branching = true;
                else
                    node->branching = false;
            }
        }
    }
}

void Skeletonizer::_verifyGraphConnectivityToMainPartition(GraphComponents &components,
                                                           SkeletonEdges &edges)
{
    // Each component in the GraphComponents is simply a list of nodes
    // Find the largest partition
    size_t primaryPartitionIndex = 0;
    size_t numberNodesPrimaryPartition = 0;
    for (size_t i = 0; i < components.size(); ++i)
    {
        if (components[i].size() > numberNodesPrimaryPartition)
        {
            primaryPartitionIndex = i;
            numberNodesPrimaryPartition = components[i].size();
        }
    }

    // Get a reference to the primary partition
    GraphComponent primaryPartition = components[primaryPartitionIndex];

    // Construct a list of the secondary partitions
    GraphComponents secondaryPartitions;
    for (size_t i = 0; i < components.size(); ++i)
    {
        if (i == primaryPartitionIndex) continue;

        secondaryPartitions.push_back(components[i]);
    }


    // Find the connections between each secondary partition and the parimary partition
    for (size_t i = 0; i < secondaryPartitions.size(); ++i)
    {
        auto secondaryPartition = secondaryPartitions[i];

        // Find the indices of the connecting nodes
        size_t closestPrimaryNodeIndex;
        size_t closesetSecondaryNodeIndex;
        float shortestDistance = 1e32;

        for (size_t j = 0; j < secondaryPartition.size(); ++j)
        {
            auto secondaryNodeIndex = secondaryPartition[j];
            auto secondaryNode = _nodes[secondaryNodeIndex];

            for (size_t k = 0; k < primaryPartition.size(); ++k)
            {
                auto primaryNodeIndex = primaryPartition[k];
                auto primaryNode = _nodes[primaryNodeIndex];

                const auto distance = primaryNode->point.distance(secondaryNode->point);
                if (distance < shortestDistance)
                {
                    shortestDistance = distance;
                    closestPrimaryNodeIndex = primaryNodeIndex;
                    closesetSecondaryNodeIndex = secondaryNodeIndex;
                }
            }
        }

        // The primary and secondary nodes are connected
        auto primaryNode = _nodes[closestPrimaryNodeIndex];
        auto secondaryNode = _nodes[closesetSecondaryNodeIndex];

        primaryNode->edgeNodes.push_back(secondaryNode);
        secondaryNode->edgeNodes.push_back(primaryNode);

        SkeletonEdge* edge = new SkeletonEdge(edges.size(), primaryNode, secondaryNode);
        edges.push_back(edge);
    }

    // Updating the branching and terminal nodes
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        // Check if the node has been visited before
        SkeletonNode* node = _nodes[i];

        if (node->edgeNodes.size() == 1)
            node->terminal = true;
        else
            node->terminal = false;

        if (node->edgeNodes.size() > 2)
            node->branching = true;
        else
            node->branching = false;
    }
}

void Skeletonizer::_verifyGraphConnectivity(SkeletonEdges& edges)
{
    // Since we have all the nodes and the edges, we can verify if the graph is conected or not
    auto graph = new Graph(edges, _nodes);

    auto components = graph->getComponents();

    if (components.size() == 1)
    {
        LOG_SUCCESS("The skeleton graph has 1 component! OK.");
    }
    else
    {
        LOG_WARNING("The skeleton graph has [ %d ] components!", components.size());

        // Verify the graph connectivity to the main partition
        _verifyGraphConnectivityToMainPartition(components, edges);

        auto newGraph = new Graph(edges, _nodes);

        auto newComponents = newGraph->getComponents();

        if (newComponents.size() == 1)
        {
            LOG_SUCCESS("The skeleton graph has now 1 component! OK.");
        }
        else
        {
            LOG_WARNING("The skeleton graph has still [ %d ] components!", newComponents.size());
        }
    }
}

void Skeletonizer::_updateParent(SkeletonBranch* branch)
{
    for(size_t j = 0; j < branch->children.size(); j++)
    {
        auto& child = branch->children[j];

        // Clear old parents if any
        child->parents.clear();
        child->parents.shrink_to_fit();

        // Add the new parent
        child->parents.push_back(branch);
        child->logicalParents.push_back(branch);

        _updateParent(child);
    }
}

void Skeletonizer::_updateParents(const bool verbose)
{
    TIMER_SET;
    VERBOSE_LOG(LOG_STATUS("Updating Parents"), verbose);
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];

        if (branch->isValid() && branch->isRoot())
        {
            for(size_t j = 0; j < branch->children.size(); j++)
            {
                auto& child = branch->children[j];

                // Clear old parents if any
                child->parents.clear();
                child->parents.shrink_to_fit();

                // Add the new parent
                child->parents.push_back(branch);
                child->logicalParents.push_back(branch);

                _updateParent(child);
            }
        }
        VERBOSE_LOG(LOOP_PROGRESS(i, _branches.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

GraphNodes Skeletonizer::_constructGraphNodesFromSkeletonNodes(
        const SkeletonNodes& skeletonNodes)
{
    // A list to conatin all the graph nodes
    GraphNodes graphNodes;

    // Resize it to be allow the parallelization of the list
    graphNodes.resize(skeletonNodes.size());

    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < skeletonNodes.size(); ++i)
    {
        const auto& skeletonNode = skeletonNodes[i];
        graphNodes[i] = (new GraphNode(skeletonNode->graphIndex,
                                       skeletonNode->point,
                                       skeletonNode->index,
                                       skeletonNode->edgeNodes.size()));
    }

    // Return the graph nodes list
    return graphNodes;
}

GraphBranches Skeletonizer::_constructGraphBranchesFromGraphNodes(
        GraphNodes &graphNodes, const int64_t& somaNodeIndex, const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Constructing Graph Branches"), verbose);

    // Use a new index to label graph branches
    size_t branchGraphIndex = 0;

    // A list of all the constructed GraphBranches
    GraphBranches graphBranches;

    // Construct the valid branches at the end
    TIMER_SET;
    for (size_t i = 0; i < graphNodes.size(); i++)
    {
        if (graphNodes[i]->children.size() > 0)
        {
            // This graph node is always the first node, becuase all the other nodes are children
            const auto& firstNodeIndex = graphNodes[i]->index;
            const auto& firstNodeSkeletonIndex = graphNodes[i]->skeletonIndex;

            for (size_t j = 0; j < graphNodes[i]->children.size(); ++j)
            {
                // This graph node is always the last node, because it is a child node
                const auto& lastNodeIndex = graphNodes[i]->children[j]->index;
                const auto& lastNodeSkeletonIndex = graphNodes[i]->children[j]->skeletonIndex;

                // Search for the branches
                for (size_t k = 0; k < _branches.size(); k++)
                {
                    // Reference to the branch
                    auto& branch = _branches[k];

                    // The branch must be valid
                    if (branch->isValid() &&
                        branch->hasTerminalNodes(firstNodeSkeletonIndex, lastNodeSkeletonIndex))
                    {
                        GraphBranch* gBranch = new GraphBranch(branchGraphIndex);
                        branchGraphIndex++;

                        gBranch->skeletonIndex = branch->index;
                        gBranch->firstNodeIndex = firstNodeIndex;
                        gBranch->firstNodeSkeletonIndex = firstNodeSkeletonIndex;
                        gBranch->lastNodeIndex = lastNodeIndex;
                        gBranch->lastNodeSkeletonIndex = lastNodeSkeletonIndex;
                        if (somaNodeIndex == firstNodeIndex) gBranch->isRoot = true;

                        graphBranches.push_back(gBranch);
                    }
                }
            }
        }
        VERBOSE_LOG(LOOP_PROGRESS(i, graphNodes.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    return graphBranches;
}

SkeletonWeightedEdges Skeletonizer::_reduceSkeletonToWeightedEdges(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Creating Simplified Weighted Graph from Skeleton"), verbose);

    /// The "thinned" skeleton has a list of edges, but indeed we would like to simplify it to
    /// avoid spending hours traversing it. Therefore, we created a "weighted skeleton", where
    /// each edge in this skeleton is a connection between two branching/terminal nodes and the
    /// weight represents the number of samples "or samples" between the two branching/terminal
    /// nodes. The "weighted edges" must be constructed only for valid branches.
    SkeletonWeightedEdges edges;
    const auto branchesCount = _branches.size();

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Constructing Weighted Edges"), verbose);
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // Reference to the current branch
        auto& branch = _branches[i];

        // The branch must be valid to be able to have a valid graph
        if (branch->isValid())
        {
            // Reset the traversal state of the "terminal" nodes of the branch
            branch->nodes.front()->visited = false;
            branch->nodes.back()->visited = false;

            // Create a weighted edge and append it to the list, where the weight is indicated by
            // the number of nodes "or samples" in the branch
            SkeletonWeightedEdge* edge = new SkeletonWeightedEdge(branch);
            edges.push_back(edge);
        }

        VERBOSE_LOG(LOOP_PROGRESS(i, branchesCount), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Return the resulting edges array that will be used for constructing the graph
    return edges;
}


SkeletonNodes Skeletonizer::_selectBranchingNodesFromWeightedEdges(
        const SkeletonWeightedEdges& edges, const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Identifying Branching Nodes"), verbose);

    // Use a new index to label the branching nodes, where the maximum value corresponds to the
    // actual number of the branching nodes in the graph
    int64_t branchingNodeIndex = 0;

    // A list to collect the branching nodes
    SkeletonNodes nodes;

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Selecting Branching Nodes for Weighted Skeleton"), verbose);
    for (size_t i = 0; i < edges.size(); ++i)
    {
        // The node must be visited once to append it to the @skeletonBranchingNodes list
        auto& edge = edges[i];
        auto node1 = edge->node1;
        auto node2 = edge->node2;

        // First node of the edge
        if (!node1->visited)
        {
            node1->graphIndex = branchingNodeIndex;
            nodes.push_back(node1);
            branchingNodeIndex++;
            node1->visited = true;
        }

        // Second node of the edge
        if (!node2->visited)
        {
            node2->graphIndex = branchingNodeIndex;
            nodes.push_back(node2);
            branchingNodeIndex++;
            node2->visited = true;
        }
        VERBOSE_LOG(LOOP_PROGRESS(i, edges.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    return nodes;
}

void Skeletonizer::_verifySkeletonNodes(const bool verbose)
{
    _totalNumberNodes = _nodes.size();
    if (_totalNumberNodes == 0)
    {
        LOG_WARNING("The skeleton has [ %ld ] centerline connected nodes!");
    }
    else
    {
        VERBOSE_LOG(LOG_SUCCESS("The skeleton has [ %ld ] centerline connected nodes. OK.",
                                _totalNumberNodes), verbose);
    }
}

void Skeletonizer::_verifySkeletonEdges(const bool verbose)
{
    _totalNumberEdges = _edges.size();
    if (_totalNumberEdges == 0)
    {
        LOG_WARNING("The skeleton has [ %ld ] centerline edges!");
    }
    else
    {
        VERBOSE_LOG(LOG_SUCCESS("The skeleton has [ %ld ] centerline edges. OK.",
                                _totalNumberEdges), verbose);
    }
}

std::vector< Vector3f > Skeletonizer::getShellPoints()
{
    return _shellPoints;
}

void Skeletonizer::_collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                                         int64_t& swcIndex, int64_t branchingNodeSWCIndex)
{
    // Get a reference to the nodes of the current branch
    auto& currentBranchNodes = branch->nodes;

    for (size_t i = 1; i < currentBranchNodes.size(); ++i)
    {
        currentBranchNodes[i]->swcIndex = swcIndex;

        if (i == 1) { currentBranchNodes[i]->prevSampleSWCIndex = branchingNodeSWCIndex;}
        else { currentBranchNodes[i]->prevSampleSWCIndex= swcIndex - 1; }

        swcIndex++;
        swcNodes.push_back(currentBranchNodes[i]);
    }

    const int64_t branchingIndex = swcIndex - 1;
    for (size_t i = 0; i < branch->children.size(); ++i)
    {
        if (branch->children[i]->isValid())
        {
            _collectSWCNodes(branch->children[i], swcNodes, swcIndex, branchingIndex);
        }
    }
}

void Skeletonizer::_exportGraphNodes(const std::string prefix, const bool verbose)
{
    // Construct the file path
    std::string filePath = prefix + NODES_EXTENSION;
    VERBOSE_LOG(LOG_STATUS("Exporting Nodes : [ %s ]", filePath.c_str()), verbose);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Writing Nodes"), verbose);
    size_t progress = 0;
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];
        if (node->terminal)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << " T " << NEW_LINE;
        }
        else
        {
        stream << node->point.x() << " "
               << node->point.y() << " "
               << node->point.z() << " "
               << node->radius << NEW_LINE;
        }

        VERBOSE_LOG(LOOP_PROGRESS(progress, _nodes.size()), verbose);
        ++progress;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Close the file
    stream.close();
}

void Skeletonizer::exportBranches(const std::string& prefix,
                                  const SkeletonBranch::BRANCH_STATE state,
                                  const bool verbose)
{
    // Construct the file path
    std::string filePath = prefix;
    std::string branchType;

    if (state == SkeletonBranch::INVALID)
    {
        filePath += INVALID_BRANCH;
        branchType = "Invalid";
    }
    else if (state == SkeletonBranch::VALID)
    {
        filePath += VALID_BRANCH;
        branchType = "Valid";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE)
    {
        filePath += TWO_SAMPLE_BRANCH;
        branchType = "Two-sample";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_VALID)
    {
        filePath += TWO_SAMPLE_VALID_BRANCH;
        branchType = "Two-sample and Valid";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_INVALID)
    {
        filePath += TWO_SAMPLE_INVALID_BRANCH;
        branchType = "Two-sample and Invalid";
    }
    else if (state == SkeletonBranch::SOMATIC)
    {
        filePath += SOMATIC_BRANCH;
        branchType = "Somatic";
    }
    else if (state == SkeletonBranch::SPINE)
    {
        filePath += SPINE_BRANCH;
        branchType = "Spine";
    }
    else
    {
        /// NOTHING
    }

    filePath += BRANCHES_EXTENSION;

    VERBOSE_LOG(LOG_STATUS("Exporting %s Branches: [ %s ]",
                           branchType.c_str(), filePath.c_str()), verbose);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    // A set of selected branches to write for debugging
    SkeletonBranches toWrite;
    if (state == SkeletonBranch::INVALID)
    {
        for (const auto& branch : _branches)
        {
            if (branch->isValid()) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE)
    {
        for (const auto& branch : _branches)
        {
            if (branch->nodes.size() == 2) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_VALID)
    {
        for (const auto& branch : _branches)
        {
            if (branch->nodes.size() == 2 && branch->isValid()) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_INVALID)
    {
        for (const auto& branch : _branches)
        {
             if (branch->nodes.size() == 2 && !branch->isValid()) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::SOMATIC)
    {
        for (const auto& branch : _branches)
        {
            if (branch->isInsideSoma()) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::SPINE)
    {
        for (const auto& branch : _branches)
        {
            if (branch->isSpine()) { toWrite.push_back(branch); }
        }
    }
    else if (state == SkeletonBranch::VALID)
    {
        for (const auto& branch : _branches)
        {
            if (branch->isValid()) { toWrite.push_back(branch); }
        }
    }
    else
    {
        for (const auto& branch : _branches)
        {
            toWrite.push_back(branch);
        }
    }

    if (toWrite.size() == 0)
    {
        LOG_WARNING("No Branches have been Collected! Aborting Export!");
        return;
    }

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Writing Branches"), verbose);
    for (size_t i = 0; i < toWrite.size(); ++i)
    {
        stream << START_BRANCH_KEYWORD << toWrite[i]->index << "\n";
        for (auto& node: toWrite[i]->nodes)
        {
            stream << node->point.x() << " " << node->point.y() << " " << node->point.z() << " "
                   << node->radius << NEW_LINE;
        }
        stream << END_BRANCH_KEYWORD << NEW_LINE;

        VERBOSE_LOG(LOOP_PROGRESS(i, toWrite.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Close the file
    stream.close();
}

}
