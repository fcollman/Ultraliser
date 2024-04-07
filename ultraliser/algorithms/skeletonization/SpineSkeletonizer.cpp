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

    // Count the number of terminal nodes
    size_t numberTerminalNodes = 0;
    for (const auto node : _nodes) { if (node->terminal) numberTerminalNodes++; }

    // This is a single branch
    if (numberTerminalNodes == 2)
    {
        // Construct the hierarchy to the terminal
        for (size_t i = 0; i < _nodes.size(); ++i)
        {
            auto& node = _nodes[i];

            // This case will handle single branches that have terminals only
            if (node->terminal)
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

                // The branch is valid
                branch->setValid();
                _branches.push_back(branch);
            }
        }
    }
    else // Branching spine
    {
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

                            // The branch is valid
                            branch->setValid();
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

                            // The branch is valid
                            branch->setValid();
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

                                // The branch is valid
                                branch->setValid();
                                _branches.push_back(branch);
                            }
                        }
                    }
                }
            }
        }
    }

    // Reset the nodes
    for (auto node : _nodes) { node->visited = false; }
}

void SpineSkeletonizer::_connectSpineBranches()
{

}


EdgesIndices SpineSkeletonizer::_findShortestPathsFromTerminalNodesToRoot(
        SkeletonWeightedEdges& edges,
        SkeletonNodes &skeletonBranchingNodes, GraphNodes &graphNodes,
        const int64_t& somaNodeIndex, const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Identifying Short Paths from Terminals To Root"), verbose);

    std::cout << "cx 1\n";
    // Identify the terminal nodes to process the paths in parallel
    SkeletonNodes terminalNodes;
    for (size_t i = 0; i < skeletonBranchingNodes.size(); i++)
    {
        if (skeletonBranchingNodes[i]->terminal)
            terminalNodes.push_back(skeletonBranchingNodes[i]);
    }
    std::cout << "cx 2\n";

    std::vector< EdgesIndices > edgesIndicesList;
    edgesIndicesList.resize(skeletonBranchingNodes.size());

    std::cout << "cx 3\n";
    // Generate the ShortestPathFinder only once for all the path retrival functions
    const ShortestPathFinder pathFinder(edges, skeletonBranchingNodes.size());

    std::cout << "cx 4\n";
    // Search for all the terminal nodes
    size_t numberTerminalNodes = terminalNodes.size();
    std::vector< PathIndices > terminalsToSomaPaths;
    terminalsToSomaPaths.resize(numberTerminalNodes);

    std::cout << "cx 5\n";
    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Searching Terminals-to-soma Paths *"), verbose);
    OMP_PARALLEL_FOR
    for (size_t iNode = 0; iNode < numberTerminalNodes; ++iNode)
    {
#ifdef REVERSE
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = pathFinder->findPath(terminalNodes[iNode]->graphIndex,
                                                       somaNodeIndex);

        // Reverse the terminal to soma path to have the correct order
        std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());
#else
        // Find the path between the terminal node and the soma node
        terminalsToSomaPaths[iNode] = pathFinder.findPath(
                    somaNodeIndex, terminalNodes[iNode]->graphIndex);
#endif
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, numberTerminalNodes), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    std::cout << "cx 6\n";
    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Updating Graph"), verbose);
    for (size_t iNode = 0; iNode < numberTerminalNodes; ++iNode)
    {
        // Get a reference to the EdgesIndices list
        EdgesIndices& perTerminalEdgesIndices = edgesIndicesList[iNode];

#ifdef REVERSE
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = pathFinder->findPath(terminalNodes[iNode]->graphIndex,
                                                       somaNodeIndex);

        // Reverse the terminal to soma path to have the correct order
        std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());
#else
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = terminalsToSomaPaths[iNode];

        // Find the edges
        for (size_t j = 0; j < terminalToSomaPath.size() - 1; ++j)
        {
            auto currentNodeIndex = terminalToSomaPath[j];
            auto nextNodeIndex = terminalToSomaPath[j + 1];

            // Add the edge indices to the list
            perTerminalEdgesIndices.push_back(EdgeIndex(currentNodeIndex, nextNodeIndex));

            // If the next node index is not in the current node index, then add it
            if (!graphNodes[currentNodeIndex]->isNodeInChildren(nextNodeIndex))
            {
                graphNodes[currentNodeIndex]->children.push_back(graphNodes[nextNodeIndex]);
            }
        }
#endif
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, numberTerminalNodes), verbose);
        VERBOSE_LOG(PROGRESS_UPDATE, verbose);
    }

    std::cout << "cx 7\n";
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Clear the terminal nodes list
    terminalNodes.clear();
    terminalNodes.shrink_to_fit();

    // Clear the terminal nodes list
    terminalsToSomaPaths.clear();
    terminalsToSomaPaths.shrink_to_fit();

    std::cout << "cx 8\n";
    // The indices of all the edges that have been traversed
    EdgesIndices edgesIndices;
    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Composing Path Edges"), verbose);
    for (size_t i = 0; i < edgesIndicesList.size(); ++i)
    {
        // Get the terminal edge identified per terminal
        auto perTerminalEdgesIndices = edgesIndicesList[i];
        if (perTerminalEdgesIndices.size() > 0)
        {
            // Append them to the edgeIndices list
            edgesIndices.insert(edgesIndices.end(),
                                perTerminalEdgesIndices.begin(), perTerminalEdgesIndices.end());

            // Clean the list per terminal
            perTerminalEdgesIndices.clear();
            perTerminalEdgesIndices.shrink_to_fit();
        }
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
    std::cout << "cx 9\n";
    // Clean the list used to collect the edges in parallel
    edgesIndicesList.clear();
    edgesIndicesList.shrink_to_fit();

    // Return the EdgesIndices list
    return edgesIndices;
}

void SpineSkeletonizer::_identifySpineBranchConnections()
{
    // If there is only a single branch, then return
    if (_branches.size() < 2) return;

    // Identify the branch connectivity
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& iBranch = _branches[i];

        // Invalid branch, ignore
        if (!iBranch->isValid()) continue;

        auto& iBranchT1 = iBranch->nodes.front();
        auto& iBranchT2 = iBranch->nodes.back();

        for (size_t j = 0; j < _branches.size(); ++j)
        {
            auto& jBranch = _branches[j];

            // Invalid branch, ignore
            if (!jBranch->isValid()) continue;

            // Same branch, next branch
            if (iBranch->index == jBranch->index) continue;

            auto& jBranchT1 = jBranch->nodes.front();
            auto& jBranchT2 = jBranch->nodes.back();

            if (iBranchT1->index == jBranchT1->index)
            {
                iBranch->t1Connections.push_back(jBranch);
            }

            if (iBranchT1->index == jBranchT2->index)
            {
                iBranch->t1Connections.push_back(jBranch);
            }

            if (iBranchT2->index == jBranchT1->index)
            {
                iBranch->t2Connections.push_back(jBranch);
            }

            if (iBranchT2->index == jBranchT2->index)
            {
                iBranch->t2Connections.push_back(jBranch);
            }
        }
    }
}

void SpineSkeletonizer::_identifyRootBranch()
{
    // _basePoint = Vector3f(0.f);

    // The shortest distance between the base point and the branches
    float shortestDistanceToBranch = std::numeric_limits< float >::max();

    // The closest branch to the base point
    _root = nullptr;

    // If a single branch is reconstructed, then make it the root
    if (_branches.size() == 0) { return; }
    if (_branches.size() == 1) { _root = _branches[0]; }
    else
    {
        // Identify the root branch based on the given base point.
        for (const auto& branch : _branches)
        {
            // std::cout << branch->nodes.front()->index << ", " << branch->nodes.back()->index << "\n";
            for (const auto& node : branch->nodes)
            {
                const auto pNode = node->point;
                const auto distance = pNode.distance(_basePoint);
                if (distance < shortestDistanceToBranch)
                {
                    shortestDistanceToBranch = distance;
                    _root = branch;
                }
            }
        }
    }

    // Identify the direction of the root branch based on the given base point
    // Compute the distance from the front and back nodes to the base point
    const auto distanceToFrontNode = _root->nodes.front()->point.distance(_basePoint);
    const auto distanceToBackNode = _root->nodes.back()->point.distance(_basePoint);

    // Reverse the orientation of the root branch only if back node is closer to the base point
    if (distanceToBackNode < distanceToFrontNode)
    {
         std::reverse(_root->nodes.begin(), _root->nodes.end());
    }

    _root->setRoot();


    _rootNode = _root->nodes.front(); // new SkeletonNode();
    // _rootNode->index = _nodes.back()->index + 1;

    // By default, the soma node is actually the soma
    _rootNode->isSoma = true;

    // The somatic node is considered inside the soma in the processing
    _rootNode->insideSoma = true;

    // Initially, we set the soma node to some value that does not make any conflict
    // _rootNode->radius = 0.1;
    // _nodes.push_back(_rootNode);

    // Update the root node
    //_rootNode->point = _root->nodes.front()->point;
    //_rootNode->radius = _root->nodes.front()->radius;
}

void constructTreeS(SkeletonBranch* branch)
{
    if (branch == nullptr)
        return;

    // std::cout << "Index " << branch->index << "\n";

    // The given branch to this function is assumed to be oriented in the right direction
    if (branch->t1Connections.size() > 0)
    {
        for (auto t1 : branch->t1Connections)
        {
            //std::cout << t1->index << " \n";
            // Ensure the correct order of the nodes in the branch

//                if (t1->nodes.front()->index != branch->nodes.back()->index)
//                {
//                    reverse(t1->nodes.begin(), t1->nodes.end());
//                }

//                 // Append the t1 to the children of the branch
//                  branch->children.push_back(t1);

//                 // Append the branch to be the parent of t1
//                  t1->parents.push_back(branch);

                // Traverse the rest of the tree
                constructTreeS(t1);

        }
    }
    else if (branch->t2Connections.size() > 0)
    {
        for (auto t2 : branch->t2Connections)
        {
            // std::cout << t2->index << " \n";
            // Ensure the correct order of the nodes in the branch

//                if (t2->nodes.front()->index != branch->nodes.back()->index)
//                {
//                     reverse(t2->nodes.begin(), t2->nodes.end());
//                }

//                // Append the t2 to the children of the branch
//                 branch->children.push_back(t2);

//                // Append the branch to be the parent of t2
//                 t2->parents.push_back(branch);

                // Traverse the rest of the tree
                constructTreeS(t2);

        }
    }
    else
    {
        /// NOTHING TO BE DONE
    }
}

void SpineSkeletonizer::_constructGraphHierarchy()
{
    if (_root == nullptr)
    {
        LOG_WARNING("Invalid Spine! No Root Branch!");
        _validSpine = false;
        return;
    }

    /// NOTE: The _root branch must be identified to
    if (_root != nullptr)
    {
        const auto t1Connections = _root->t1Connections;
        const auto t2Connections = _root->t2Connections;

        // If the root node has connecting branches along both terminals, then that's an issue
        if (t1Connections.size() > 0 && t2Connections.size() > 0)
        {
            // Report and Issue
            LOG_WARNING("Spine root connected from both terminals!");
            _validSpine = false;
        }
        else
        {

            constructTreeS(_root);
        }
    }
}


int64_t SpineSkeletonizer::_getRootNodeIndexFromGraphNodes(const SkeletonNodes& nodes) const
{
    int64_t rootNodeIndex = -1;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        const auto& node = nodes[i];
        if (node->isSoma)
        {
            rootNodeIndex = node->graphIndex;
            break;
        }
    }

    // Return the soma node index
    return rootNodeIndex;
}


void SpineSkeletonizer::_addRootNode()
{
    _rootNode = new SkeletonNode();
    _rootNode->index = _nodes.back()->index + 1;

    // By default, the soma node is actually the soma
    _rootNode->isSoma = true;

    // The somatic node is considered inside the soma in the processing
    _rootNode->insideSoma = true;

    // Initially, we set the soma node to some value that does not make any conflict
    _rootNode->radius = 0.1;
    _nodes.push_back(_rootNode);
}

void SpineSkeletonizer::_constructSkeletonHierarchy(GraphBranches& graphBranches,
                                                     const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Constructing Skeleton Hierarchy"), verbose);

    TIMER_SET;
    for(size_t i = 0; i < graphBranches.size(); ++i)
    {
        // Reference to the GraphBranch
        const auto& graphBranch = graphBranches[i];

        // Reference to the corresponding SkeletonBranch
        auto& skeletonBranch = _branches[graphBranch->skeletonIndex];

        // Adjust the direction, i.e. the samples list
        skeletonBranch->adjustDirection(graphBranch->firstNodeSkeletonIndex,
                                        graphBranch->lastNodeSkeletonIndex);

        // Update the tree
        for(size_t j = 0; j < graphBranch->children.size(); ++j)
        {
            skeletonBranch->children.push_back(
                        _branches[graphBranch->children[j]->skeletonIndex]);
            skeletonBranch->logicalChildren.push_back(
                        _branches[graphBranch->children[j]->skeletonIndex]);
        }

        VERBOSE_LOG(LOOP_PROGRESS(i, graphBranches.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void SpineSkeletonizer::_constructGraphHierarchy(GraphBranches& graphBranches, const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Constructing Graph Hierarchy"), verbose);
    TIMER_SET;
    for (size_t i = 0; i < graphBranches.size(); ++i)
    {
        auto& iBranch = graphBranches[i];

        // Get the last node of the iBranch
        const auto& iLastNodeIndex = iBranch->lastNodeIndex;

        for (size_t j = 0; j < graphBranches.size(); ++j)
        {
            auto& jBranch = graphBranches[j];

            // If the same branch, next branch
            if (iBranch->index == jBranch->index) continue;

            // Get the first node of the jBranch
            const auto& jFirstNodeIndex = jBranch->firstNodeIndex;

            // jBranch is a child of the iBranch if the nodes are the same
            if (iLastNodeIndex == jFirstNodeIndex)
            {
                iBranch->children.push_back(jBranch);
            }
        }

        VERBOSE_LOG(LOOP_PROGRESS(i, graphBranches.size()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}


void SpineSkeletonizer::_updateParent(SkeletonBranch* branch)
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

void SpineSkeletonizer::_updateParents(const bool verbose)
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


void SpineSkeletonizer::run(const bool verbose)
{

    if (_index == 0) return;
    // Initialize
    initialize(verbose);

    std::stringstream prefix;
    prefix << "/ssd2/skeletonization-project/spine-extraction/output/refacotr-1/spines/864691134832191490_" << _index << "_volume_";
    // _volume->project(prefix.str(), true);

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

    _removeTriangleLoops(false);

    // export nodes
    _exportGraphNodes(prefix.str(), false);



    /// Inflate the nodes, i.e. adjust their radii
    _inflateNodes(verbose);

    // _addRootNode();




    /// Reconstruct the sections "or branches" from the nodes using the edges data
    _buildSpineBranchesFromNodes();


    // Connect the spine branches
    _identifySpineBranchConnections();

    identifyTerminalBranchesForSpine(_branches);


    // Identify the root node
    _identifyRootBranch();

    if (_root == nullptr || !_validSpine) return;

    if (_branches.size() > 1)
    {
        /// Reduce the skeleton into a list of SkeletonWeightedEdge's
        SkeletonWeightedEdges weighteEdges = _reduceSkeletonToWeightedEdges(verbose);

        /// Get a list of all the branching/terminal nodes within the skeleton from the
        /// SkeletonWeightedEdges list
        SkeletonNodes skeletonBranchingNodes = _selectBranchingNodesFromWeightedEdges(weighteEdges,
                                                                                      verbose);

        /// Get the soma node index within the weighted graph
        int64_t rootNodeIndex = _getRootNodeIndexFromGraphNodes(skeletonBranchingNodes);

        // Construct the graph nodes list
        GraphNodes graphNodes = _constructGraphNodesFromSkeletonNodes(skeletonBranchingNodes);

        /// After having the weighted edges and the nodes computed, compute the number of components in
        /// the graph, if the result is more than 1 then then re-connect them to be in a single graph
        auto graph = new Graph(weighteEdges, graphNodes);
        auto components = graph->getComponents();
//        if (components.size() == 1)
//        {
//            LOG_SUCCESS("The skeleton graph has 1 component! OK.");
//        }
//        else
//        {
//            LOG_WARNING("The skeleton graph has [ %d ] components!", components.size());
//        }

        std::cout << "1\n";
        // Find the shortest paths of all the terminals and get a list of the indices of the active edges
        EdgesIndices edgeIndices = _findShortestPathsFromTerminalNodesToRoot(
                    weighteEdges, skeletonBranchingNodes, graphNodes, rootNodeIndex, verbose);

        std::cout << "2\n";

        // Construct the GraphBranches from the GraphNodes
        GraphBranches graphBranches = _constructGraphBranchesFromGraphNodes(
                    graphNodes, rootNodeIndex, verbose);

        std::cout << "3\n";

        // Construct the hierarchy of the graph
        _constructGraphHierarchy(graphBranches, verbose);

        std::cout << "4\n";

        // Construct the hierarchy of the skeleton
        _constructSkeletonHierarchy(graphBranches, verbose);


        std::cout << "5\n";

        // Update all the parents
        _updateParents(verbose);

    }

    /// todo: remove triangles edges

//    std::cout << "\nIndex: " <<  _index << ": \n";
//    for (auto branch : _branches)
//    {
//        std::cout << branch->t1Connections.size() << " " << branch->t2Connections.size();// << "\n";
//        if (branch->isTerminal())
//            std::cout << " T \n";
//        else
//            std::cout << "\n";
//    }


    // _constructGraphHierarchy();



     exportBranches(prefix.str(), verbose);

      exportSWCFile(prefix.str(), false, SILENT);


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
    fakeSomaNode->point.x() = _root->nodes.front()->point.x();
    fakeSomaNode->point.y() = _root->nodes.front()->point.y();
    fakeSomaNode->point.z() = _root->nodes.front()->point.z();
    fakeSomaNode->radius = _root->nodes.front()->radius;

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
    if (!_validSpine)
        return;

    // Construct the file path
    std::string filePath = prefix + SWC_EXTENSION;
    VERBOSE_LOG(LOG_STATUS("Exporting Spine to SWC file: [ %s ]", filePath.c_str()), verbose);

    TIMER_SET;

    auto swcNodes = _constructSWCTable(resampleSkeleton);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    // Create a fake node for the soma, simply from the root branch
    auto somaNode = swcNodes[0];
    stream << somaNode->swcIndex << " "
           << SWC_SOMA_STRUCT_IDENTIFIER << " "
           << somaNode->point.x() << " "
           << somaNode->point.y() << " "
           << somaNode->point.z() << " "
           << somaNode->radius << " "
           << "-1" << "\n";

    VERBOSE_LOG(LOOP_STARTS("Writing SWC Table"), verbose);
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
        VERBOSE_LOG(LOOP_PROGRESS(i, numberSWCNodes), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Close the file
    stream.close();
}

}
