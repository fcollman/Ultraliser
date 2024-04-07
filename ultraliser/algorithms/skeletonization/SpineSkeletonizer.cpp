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


#include <algorithms/skeletonization/SpineSkeletonizer.h>
#include <algorithms/skeletonization/SkeletonizerUtils.h>

namespace Ultraliser
{

SpineSkeletonizer::SpineSkeletonizer(Volume* spineVolume,
                                     const bool useAcceleration,
                                     const bool debugSkeleton,
                                     const std::string debugPrefix)
    : Skeletonizer(spineVolume, useAcceleration, debugSkeleton, debugPrefix)
{
    /// EMPTY CONSTRUCTOR
}

SpineSkeletonizer::SpineSkeletonizer(Mesh* spineMesh,
                                     const size_t& index,
                                     const Vector3f basePoint,
                                     const VoxelizationOptions& options,
                                     const bool useAcceleration,
                                     const bool debugSkeleton,
                                     const std::string debugPrefix)
    : Skeletonizer(spineMesh, options, useAcceleration, debugSkeleton, debugPrefix)
{
    _index = index;
    _basePoint = basePoint;
    _computeVolumeFromMesh();
}

void adjustRootBranchOrientation(SkeletonBranch* rootBranch, const SkeletonNode* rootNode)
{
    // We can rely on the front node only to verify
    auto firstNode = rootBranch->nodes.front();

    // If the first node of the branch is the root node, then return as nothing to be fixed
    if (firstNode->index == rootNode->index)
        return;

    // Otherwise, reverse the order of the nodes
    std::reverse(std::begin(rootBranch->nodes), std::end(rootBranch->nodes));
}

void adjustChildrenBranchOrientation(SkeletonBranch* rootBranch, SkeletonBranches& branches)
{
    // The root branch should be okay
    // Get the last sample of the root branch
    auto lastNodeRootBranch = rootBranch->nodes.back();

    // If the last node has no edge nodes, return
    if (lastNodeRootBranch->edgeNodes.size() == 0)
        return;

    // Find the branches that have the same terminal node
    for (size_t i = 0; i < branches.size(); ++i)
    {
        auto branch = branches[i];

        // If the same branch, then continue
        if (branch->index == rootBranch->index) continue;

        // If the first node of the branch matches the last node of the root branch, this is a direct child
        auto firstNode = branch->nodes.front();
        if (lastNodeRootBranch->index == firstNode->index)
        {
            rootBranch->children.push_back(branch);
            branch->parents.push_back(rootBranch);

            // Continue
            continue;
        }

        // Otherwise, check if the last node
        auto lastNode = branch->nodes.back();
        if (lastNodeRootBranch->index == lastNode->index)
        {
            // We must fix the branch first
            std::reverse(std::begin(branch->nodes), std::end(branch->nodes));

            // Then we can append
            rootBranch->children.push_back(branch);
            branch->parents.push_back(rootBranch);
        }
    }

    // Do it recursively
    for (size_t i = 0; i < rootBranch->children.size(); ++i)
    {
        adjustChildrenBranchOrientation(rootBranch->children[i], branches);
    }
}


void SpineSkeletonizer::_buildSpineBranchesFromNodes()
{
    // Used to index the branch
    size_t branchIndex = 0;

    // Construct the hierarchy to the terminal
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];

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

        // This case will handle single branches that have terminals only
        else if (node->terminal)
        {
            if (node->visited) continue;

            // Create a new branch and set its index to zero
            SkeletonBranch* branch = new SkeletonBranch();
            branch->index = branchIndex++;

            SkeletonNode *currentNode = node;
            branch->nodes.push_back(currentNode);
            currentNode->visited = true;

            SkeletonNode* nextNode = currentNode->edgeNodes[0];

            while(true)
            {
                if (nextNode->visited) break;

                // Add the current node to branch
                branch->nodes.push_back(nextNode);
                nextNode->visited = true;

                // The other terminal node
                if (nextNode->edgeNodes.size() == 1)
                {
                    nextNode = nextNode->edgeNodes[0];
                    branch->nodes.push_back(nextNode);
                    nextNode->visited = true;
                    break;
                }
                else
                {
                    auto edgeNode0 = nextNode->edgeNodes[0];
                    auto edgeNode1 = nextNode->edgeNodes[1];

                    if (edgeNode0->visited)
                        nextNode = edgeNode1;
                    else
                        nextNode = edgeNode0;

                }
            }
            _branches.push_back(branch);
        }

    }
}


void SpineSkeletonizer::run(const bool verbose)
{
    // Initialize
    initialize(verbose);

    std::stringstream prefix;
    prefix << "/ssd2/skeletonization-project/spine-extraction/output/refacotr-1/spines/864691134832191490_" << _index << "_volume_";
    _volume->project(prefix.str(), true);

    // Skeletonize the volume to center-lines
    skeletonizeVolumeToCenterLines(verbose);

    prefix << SKELETON_SUFFIX;
    //_volume->project(prefix.str(), true);

    /// Extract the nodes of the skeleton from the center-line "thinned" voxels and return a
    /// mapper that maps the indices of the voxels in the volume and the nodes in the skeleton
    auto indicesMapper = _extractNodesFromVoxels(verbose);


    /// Connect the nodes of the skeleton to construct its edges. This operation will not connect
    /// any gaps, it will just connect the nodes extracted from the voxels.
    _connectNodesToBuildEdges(indicesMapper, verbose);


    // export nodes
    _exportGraphNodes(prefix.str(), false);



    /// Inflate the nodes, i.e. adjust their radii
    _inflateNodes(verbose);

    /// Reconstruct the sections "or branches" from the nodes using the edges data
    _buildSpineBranchesFromNodes();

    /// todo: remove triangles edges





     exportBranches(prefix.str(), verbose);


    // Identify the connections at the terminals of each branch
    // identifyTerminalConnections(_branches);

    //
    // identifyTerminalBranchesForSpine(_branches);

    // Identify the root branch
//    auto rootBranch = identifyRootBranchForSpine(_branches, _basePoint);
//    if (rootBranch == nullptr)
//        std::cout << "NULL ROOT BRANCH \n";

//    // Identify the root node
//    auto rootNode = identifyRootNodeForSpine(rootBranch, _basePoint);
//    if (rootNode == nullptr)
//        std::cout << "NULL ROOT NODE\n";

//    // Verify the orientation of the root branch
//    // mainly for the root branch
//    adjustRootBranchOrientation(rootBranch, rootNode);


//    // If the spine has more than one branch
//    adjustChildrenBranchOrientation(rootBranch, _branches);




    // Identify the paths of the morphology

    // Construct the tree of the spine

    // Compute the orientation of the spine

    // exportSWCFile(_debugPrefix, false, SILENT);



}

void SpineSkeletonizer::segmentComponents(const bool verbose)
{

}

void SpineSkeletonizer::exportBranches(const std::string& prefix, const bool)
{

    // Construct the file path
    std::string filePath = prefix + BRANCHES_EXTENSION;

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // The @start marks a new branch in the file
        stream << "SB " << _branches[i]->index << "\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " " << node->point.y() << " " << node->point.z() << " "
                   << node->radius << NEW_LINE;
        }
        // The @end marks the terminal sample of a branch
        stream << "EB\n";
    }

    // Close the file
    stream.close();
}

void SpineSkeletonizer::_collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                                         int64_t& swcIndex, int64_t branchingNodeSWCIndex, const bool)
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

SkeletonNodes SpineSkeletonizer::_constructSWCTable(const bool& resampleSkeleton, const bool)
{
    // A table, or list that contains all the nodes in order
    SkeletonNodes swcNodes;

    // A global index that will be used to correctly index the nodes
    int64_t swcIndex = 1;

    // Fake soma node
    SkeletonNode* fakeSomaNode = new SkeletonNode();

    // Append the somatic mode
    fakeSomaNode->swcIndex = swcIndex;
    fakeSomaNode->prevSampleSWCIndex = -1;
    fakeSomaNode->point.x() = 0;
    fakeSomaNode->point.y() = 0;
    fakeSomaNode->point.z() = 0;
    fakeSomaNode->radius = 0;

    swcIndex++;
    swcNodes.push_back(fakeSomaNode);

    // Resample the skeleton
    TIMER_SET;
    if (resampleSkeleton)
    {
        LOOP_STARTS("Resampling Skeleton");
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            auto& branch = _branches[i];

            // Do not resample the root sections
            if (branch->isRoot()) continue;

            // Resample only valid branches
            if (branch->isValid()) { branch->resampleAdaptively(); }
            LOOP_PROGRESS(i, _branches.size());
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Get all the root branches
    TIMER_RESET;
    LOOP_STARTS("Constructing SWC Table");
    const size_t numberBranches = _branches.size();
    for (size_t i = 0; i < numberBranches ; ++i)
    {

        auto& branch = _branches[i];
        if (branch->isRoot() && branch->isValid())
        {
            // The branching index is that of the soma
            _collectSWCNodes(branch, swcNodes, swcIndex, 1);
        }
        LOOP_PROGRESS(i, numberBranches);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return swcNodes;
}

void SpineSkeletonizer::exportSWCFile(const std::string& prefix, const bool& resampleSkeleton, const bool verbose)
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + SWC_EXTENSION;
    LOG_STATUS("Exporting Spine to SWC file: [ %s ]", filePath.c_str());

    auto swcNodes = _constructSWCTable(resampleSkeleton);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    auto somaNode = swcNodes[0];
    stream << somaNode->swcIndex << " "
           << SWC_SOMA_STRUCT_IDENTIFIER << " "
           << somaNode->point.x() << " "
           << somaNode->point.y() << " "
           << somaNode->point.z() << " "
           << somaNode->radius << " "
           << "-1" << "\n";

    LOOP_STARTS("Writing SWC Table");
    const size_t numberSWCNodes = swcNodes.size();
    for (size_t i = 1; i < numberSWCNodes; ++i)
    {
        // TODO: Export all the branches as basal dendrites UFN
        auto swcNode = swcNodes[i];
        stream << swcNode->swcIndex << " "
               << SWC_BASAL_DENDRITE_STRUCT_IDENTIFIER << " "
               << swcNode->point.x() << " "
               << swcNode->point.y() << " "
               << swcNode->point.z() << " "
               << swcNode->radius << " "
               << swcNode->prevSampleSWCIndex << "\n";
    LOOP_PROGRESS(i, numberSWCNodes);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

}
