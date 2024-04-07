/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include "SpineMorphology.h"
#include <common/Common.h>
#include <utilities/TypeConversion.h>
#include <utilities/Utilities.h>
#include <algorithms/mcs/DualMarchingCubes.h>

namespace Ultraliser
{

void SpineMorphology::_constructTreeFromLogicalBranches(SkeletonBranch* root,
                                                        size_t& sectionIndex)
{
    // Adding the terminal nodes to visualize the terminals of the branches
    if (root->logicalChildren.size() == 0)
    {
        _terminals.push_back(root->nodes.back()->point);
    }

    for (size_t i = 0; i < root->logicalChildren.size(); ++i)
    {
        Section* section = new Section(sectionIndex++);
        for (size_t j = 0; j < root->logicalChildren[i]->nodes.size(); ++j)
        {
            auto node = root->logicalChildren[i]->nodes[j];
            section->addSample(new Sample(node->point, node->radius * _radfiusScaleFactor, j));
        }
        _sections.push_back(section);

        // Traverse the children tree
        _constructTreeFromLogicalBranches(root->logicalChildren[i], sectionIndex);
    }
}

SpineMorphology::SpineMorphology(SkeletonBranch* root, const size_t &spineIndex)
{
    // Update the base point
    _basePoint = root->nodes[0]->point;

    // Update the root sampe
    _rootSample = new Sample(root->nodes[0]->point, root->nodes[0]->radius, 0);

    // Set the spine index to be the same index of the spineRoot in the array
    _spineIndex = spineIndex;

    // The section index should be the same as the index of the section in the list
    size_t sectionIndex = 0;

    // Add the root to the list of the sections
    Section* section = new Section(sectionIndex++);
    for (size_t i = 0; i < root->nodes.size(); ++i)
    {
        auto node = root->nodes[i];
        section->addSample(new Sample(node->point, node->radius * _radfiusScaleFactor, i));
    }
    _sections.push_back(section);

    // Traverse the children tree
    _constructTreeFromLogicalBranches(root, sectionIndex);

    // Compute the bounding box of the entire morphology
    _computeBoundingBox();

    // Compute the orientation vector
    _computeOrientationVector();
}

Sections SpineMorphology::_getNonDendrticSections() const
{
    Sections nonDendriticSections;

    auto dendriticCenter = _rootSample[0].getPosition();
    auto dendriticExtent = _rootSample[0].getRadius();

    size_t validSectionIndex = 0;

    for (size_t i = 0; i < _sections.size(); ++i)
    {
        auto section = _sections[i];
        size_t numberValidSamples = 0;
        for (size_t j = 0; j < section->getSamples().size(); ++j)
        {
            auto sample = section->getSamples()[j];
            auto sampleCenter = sample->getPosition();
            if (!Utilities::isPointInsideSphere(sampleCenter, dendriticCenter, dendriticExtent))
            {
                numberValidSamples++;
            }
        }

        // If more than a sample, then it is a valid section
        if (numberValidSamples > 1)
        {
            // Construct the section
            Section* validSection = new Section(validSectionIndex);
            validSectionIndex++;

            for (size_t j = 0; j < section->getSamples().size(); ++j)
            {
                auto sample = section->getSamples()[j];
                auto sampleCenter = sample->getPosition();
                if (!Utilities::isPointInsideSphere(sampleCenter, dendriticCenter, dendriticExtent))
                {
                    Sample* newSample = new Sample(sample->getPosition(), sample->getRadius() * 1, j);
                    validSection->addSample(newSample);
                }
            }
            nonDendriticSections.push_back(validSection);
        }
    }

    return nonDendriticSections;
}

SpineMorphology::SpineMorphology(SkeletonBranches branches, const size_t &index)
{
    // Update the spine index
    _spineIndex = index;

    for (size_t i = 0; i < branches.size(); ++i)
    {
        const auto branch = branches[i];
        Section* section = new Section(i);
        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            auto node = branch->nodes[j];
            section->addSample(new Sample(node->point, node->radius, j));
        }
        _sections.push_back(section);
    }

    // Compute the bounding box of the entire morphology
    _computeBoundingBox();

    // Since we do not know where the root branch is, use the center of the bounds as the base point
    _basePoint = _pMin + 0.5 * (_pMax - _pMin);

    // Compute the orientation vector
    _computeOrientationVector();
}

void SpineMorphology::_computeBoundingBox()
{
    // Bounding box data
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(-1 * std::numeric_limits<float>::max());

    for (const auto& section: _sections)
    {
        for (const auto& sample: section->getSamples())
        {
            const auto position = sample->getPosition();
            const auto radius = sample->getRadius();

            Vector3f pMaxSample = position + Vector3f(radius);
            Vector3f pMinSample = position - Vector3f(radius);

            if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
            if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
            if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

            if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
            if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
            if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
        }
    }
}

void SpineMorphology::_computeOrientationVector()
{
    Vector3f orientation(0.);
    size_t numberValidSegments = 0;
    for (const auto& section : _sections)
    {
        // If the section has under any conditions less than two samples, continue
        if (section->getSamples().size() < 2) continue;

        // Compute the orientation of each edge (or segment) in the section
        for (size_t i = 0; i < section->getSamples().size() - 1; ++i)
        {
            const auto& p1 = section->getSamples()[i]->getPosition();
            const auto& p2 = section->getSamples()[i + 1]->getPosition();

            const auto segmentOrientation = p1.directionTo(p2);
            orientation += segmentOrientation;

            // One more valid segment
            numberValidSegments++;
        }
    }

    // Normalize
    orientation /= numberValidSegments;
    _direction = orientation;
}

Volume* SpineMorphology::generateVolumeWithoutDenderiticExtent(const float& voxelsPerMicron,
                                                               const float& edgeGap,
                                                               const bool & verbose)
{
    auto nonDendriticSections = _getNonDendrticSections();

    if (nonDendriticSections.size() == 0)
        return nullptr;

    // Bounding box data
    Vector3f pMinInput = Vector3f(std::numeric_limits<float>::max());
    Vector3f pMaxInput = Vector3f(-1 * std::numeric_limits<float>::max());

    for (const auto& section: nonDendriticSections)
    {
        for (const auto& sample: section->getSamples())
        {
            const auto position = sample->getPosition();
            const auto radius = sample->getRadius();

            Vector3f pMaxSample = position + Vector3f(radius);
            Vector3f pMinSample = position - Vector3f(radius);

            if (pMaxSample.x() > pMaxInput.x()) pMaxInput.x() = pMaxSample.x();
            if (pMaxSample.y() > pMaxInput.y()) pMaxInput.y() = pMaxSample.y();
            if (pMaxSample.z() > pMaxInput.z()) pMaxInput.z() = pMaxSample.z();

            if (pMinSample.x() < pMinInput.x()) pMinInput.x() = pMinSample.x();
            if (pMinSample.y() < pMinInput.y()) pMinInput.y() = pMinSample.y();
            if (pMinSample.z() < pMinInput.z()) pMinInput.z() = pMinSample.z();
        }
    }

    Vector3f inputBB = pMaxInput - pMinInput;
    Vector3f inputCenter = pMinInput + (0.5f * inputBB);

    // Expand the bounding box to be able to caprture all the details missed
    pMinInput -= edgeGap * inputBB;
    pMaxInput += edgeGap * inputBB;
    inputBB = pMaxInput - pMinInput;
    inputCenter = pMinInput + (0.5f * inputBB);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, edgeGap, VOLUME_TYPE::UI8, verbose);

    // Rasterize the morphologies into the volume
    volume->surfaceVoxelizeSections(nonDendriticSections, verbose);
    volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ, verbose);

    // Return the volume
    return volume;
}
Volume* SpineMorphology::generateVolume(const float& voxelsPerMicron,
                                        const float& edgeGap,
                                        const bool & verbose)
{
    // Get the bounding box of the morphology
    Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
    getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

    // Expand the bounding box to be able to caprture all the details missed
    _pMin -= edgeGap * inputBB;
    _pMax += edgeGap * inputBB;
    inputBB = _pMax - _pMin;
    inputCenter = _pMin + (0.5f * inputBB);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, edgeGap, VOLUME_TYPE::UI8, verbose);

    // Rasterize the morphologies into the volume
    volume->surfaceVoxelizeSpineMorphology(this, POLYLINE_SPHERE_PACKING);
    volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ, verbose);

    // Return the volume
    return volume;
}

Mesh* SpineMorphology::generateMesh(const float &voxelsPerMicron,
                                       const float& edgeGap,
                                       const bool &verbose)
{
    // Reconstruct the volume
    auto volume = generateVolume(voxelsPerMicron, edgeGap, verbose);

    // Use the DMC algorithm to reconstruct a mesh
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume, verbose);

    // Smooth the mesh to be able to have correct mapping
    mesh->smoothSurface(10, verbose);

    // Return the mesh
    return mesh;
}

Mesh* SpineMorphology::generateMeshWithoutDenderiticExtent(const float& voxelsPerMicron,
                                                           const float& edgeGap,
                                                           const bool& verbose)
{
    // Reconstruct the volume
    auto volume = generateVolumeWithoutDenderiticExtent(voxelsPerMicron, edgeGap, verbose);

    // Instead of returning a nullptr, generate a spine with the dendritic extent included
    if (volume == nullptr)
    {
        // TODO: FIX ME!
        return nullptr;
        volume = generateVolume(voxelsPerMicron, edgeGap, verbose);
    }

    // Use the DMC algorithm to reconstruct a mesh
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume, verbose);

    // Smooth the mesh to be able to have correct mapping
    mesh->smoothSurface(DEFAULT_SMOOTHING_ITERATIONS, verbose);
    mesh->subdivideTrianglseAtCentroid();

    // Return the mesh
    return mesh;
}

void SpineMorphology::exportExtents(const std::string& prefix) const
{
    // Construct the file path
    std::stringstream sstream;
    sstream << prefix << SPINE_SUFFIX << "-" <<_spineIndex << EXTENTS_EXTENSION;
    std::fstream stream;
    stream.open(sstream.str(), std::ios::out);

    // Compute the BB
    const auto bounds = _pMax - _pMin;
    const auto center = _pMin + bounds * 0.5f;

    // Export the data
    stream << center.x() << " " << center.y() << " " << center.z() << " "
           << bounds.x() << " " << bounds.y() << " " << bounds.z() << NEW_LINE;

    // Close the file
    stream.close();
}

void SpineMorphology::exportBranches(const std::string& prefix,
                                     const bool& verbose)
{
    // Start the timer
    TIMER_SET;

    // Prefix
    std::stringstream sstream;
    sstream << prefix << _spineIndex << BRANCHES_EXTENSION;
    std::string filePath = sstream.str();
    VERBOSE_LOG(LOG_STATUS("Exporting Spine Branches: [ %s ]", filePath.c_str()), verbose);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    VERBOSE_LOG(LOOP_STARTS("Writing Spine Branches"), verbose);
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        const auto& section = _sections[i];

        // The @start marks a new branch in the file
        stream << START_BRANCH_KEYWORD << section->getIndex() << NEW_LINE;
        for (auto& sample: section->getSamples())
        {
            stream << sample->getPosition().x() << " "
                   << sample->getPosition().y() << " "
                   << sample->getPosition().z() << " "
                   << sample->getRadius() << NEW_LINE;
        }
        // The @end marks the terminal sample of a branch
        stream << END_BRANCH_KEYWORD << NEW_LINE;

        VERBOSE_LOG(LOOP_PROGRESS(i, _sections.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Close the file
    stream.close();
}
SpineMorphology::~SpineMorphology()
{
    /// EMPTY DESTRUCTOR
}

}
