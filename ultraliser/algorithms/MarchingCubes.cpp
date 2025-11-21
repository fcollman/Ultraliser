/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include "MarchingCubes.h"
#include "MarchingCubes.hh"
#include <data/volumes/voxels/Voxels.h>
#include <data/meshes/simple/Mesh.h>
#include <set>

namespace Ultraliser
{

MarchingCubes::MarchingCubes(Volume* volume,
                             const size_t _isoValue,
                             const bool keepOpenBoundaries)
    : _volume(volume)
    , _isoValue(_isoValue)
    , _keepOpenBoundaries(keepOpenBoundaries)
{
    /// EMPTY CONSTRUCTOR
}

double interpolateIsoValue(double _isoValue, double f1, double f2, double t1, double t2)
{
    if(f2 == f1)
        return 0.5 * (t2 + t1);

    // Return interpolated
    return (t2 - t1) * (_isoValue - f1) / (f2 - f1) + t1;
}

size_t addSharedVertex(double x1, double y1, double z1,
                       double c2,
                       int axis, double f1, double f2,
                       double _isoValue,
                       Vertices& vertices)
{
    size_t vertexIndex = vertices.size();

    if(axis == 0)
    {
        double x = interpolateIsoValue(_isoValue, f1, f2, x1, c2);
        vertices.push_back(Vector3f(x, y1, z1));
        return vertexIndex;
    }

    if(axis == 1)
    {
        double y = interpolateIsoValue(_isoValue, f1, f2, y1, c2);
        vertices.push_back(Vector3f(x1, y, z1));
        return vertexIndex;
    }

    if(axis == 2)
    {
        double z = interpolateIsoValue(_isoValue, f1, f2, z1, c2);
        vertices.push_back(Vector3f(x1, y1, z));
        return vertexIndex;
    }

    // Just for complication
    return 0;
}

Mesh* MarchingCubes::generateMesh()
{
    LOG_TITLE("Mesh Reconstruction with Marching Cubes");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    
    // Track which vertices are on the volume boundary
    std::set<size_t> borderVertexIndices;
    _buildSharedVertices(vertices, triangles, borderVertexIndices);

    // Reconstruct the mesh
    Mesh* mesh = new Mesh(vertices, triangles);
    
    // Mark border vertices in the mesh
    mesh->markBorderVertices(borderVertexIndices);

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

AdvancedMesh* MarchingCubes::generateAdvancedMesh()
{
    LOG_TITLE("Mesh Reconstruction with Marching Cubes");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    
    // Build the mesh with border vertex detection
    _borderVertexIndices.clear();
    _buildSharedVertices(vertices, triangles, _borderVertexIndices);

    // Reconstruct the mesh
    AdvancedMesh* mesh = new AdvancedMesh(vertices, triangles);

    // Mark border vertices if any were detected
    if (!_borderVertexIndices.empty())
    {
        mesh->markBorderVertices(_borderVertexIndices);
    }

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

void MarchingCubes::_buildSharedVertices(Vertices& vertices, Triangles &triangles, std::set<size_t>& borderVertexIndices)
{
    // Start the timer
    TIMER_SET;

    // Polygons
    std::vector< size_t > polygons;
    
    // Clear border vertex indices
    borderVertexIndices.clear();
    
    // Get volume dimensions for boundary checking
    const int64_t volumeWidth = static_cast<int64_t>(_volume->getWidth());
    const int64_t volumeHeight = static_cast<int64_t>(_volume->getHeight());
    const int64_t volumeDepth = static_cast<int64_t>(_volume->getDepth());
    const int64_t width = volumeWidth;
    const int64_t height = volumeHeight;
    const int64_t depth = volumeDepth;

    // Adding a little bit of extra voxels
    const int64_t extraVoxels = 2;
    const int64_t minValue = -1 * extraVoxels;
    const int64_t maxValue = extraVoxels;

    const int64_t maxX = _volume->getWidth() + maxValue;
    const int64_t maxY = _volume->getHeight() + maxValue;
    const int64_t maxZ = _volume->getDepth() + maxValue;

    const size_t sizeX = maxX + extraVoxels;
    const size_t sizeY = maxY + extraVoxels;
    const size_t sizeZ = maxZ + extraVoxels;

    // Total size
    const size_t size = (sizeX * sizeY * sizeZ * 3);

    // Mesh shared indices list
    size_t* sharedIndices;
    try {
        sharedIndices = new size_t[size];

    }  catch (...) {
        LOG_ERROR("Cannot allocate required memory for the Marching Cubes implementation!");
    }

    // A list of lists of DMCVoxel's
    // This list will have reducedX entries to get filled in parallel
    MCVoxelsList volumeMCVoxels;
    volumeMCVoxels.resize(sizeX);

    // For computations and indexing ...
    const size_t z3 = sizeZ * 3;
    const size_t yz3 = sizeY * z3;

    LOOP_STARTS("Searching Filled Voxels");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = minValue; i < maxX; ++i)
    {
        // Get a reference to the slice
        MCVoxels& sliceMCVoxels = volumeMCVoxels[static_cast< size_t >(i + extraVoxels)];

        for (int64_t j = minValue; j < maxY; ++j)
        {
            for (int64_t k = minValue; k < maxZ; ++k)
            {
                double v[8];
                v[0] = _volume->getValueUI64(i, j, k);
                v[1] = _volume->getValueUI64(i + 1, j, k);
                v[2] = _volume->getValueUI64(i + 1, j + 1, k);
                v[3] = _volume->getValueUI64(i, j + 1, k);
                v[4] = _volume->getValueUI64(i, j, k + 1);
                v[5] = _volume->getValueUI64(i + 1, j, k + 1);
                v[6] = _volume->getValueUI64(i + 1, j + 1, k + 1);
                v[7] = _volume->getValueUI64(i, j + 1, k + 1);

                // Get the cube index
                size_t cubeIndex = 0;
                for (uint8_t l = 0; l < 8; ++l)
                    if (v[l] <= _isoValue)
                        cubeIndex |= 1 << l;

                if (cubeIndex == 0 || cubeIndex == 255)
                    continue;

                sliceMCVoxels.push_back(new MCVoxel(i, j, k, cubeIndex, v));
            }
        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, sizeX);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Delta between the voxels
    const float voxelSize = _volume->getVoxelSize();

    // Building the shared vertices
    TIMER_RESET;
    LOOP_STARTS("Building Shared Vertices");
    for (size_t ii = 0; ii < volumeMCVoxels.size(); ii++)
    {
        LOOP_PROGRESS(ii, volumeMCVoxels.size());
        for (size_t jj = 0; jj < volumeMCVoxels[ii].size(); jj++)
        {
            // Reference to the voxel
            MCVoxel* mcVoxel = volumeMCVoxels[ii][jj];

            // Voxel config
            double* v = mcVoxel->vConfig;

            const int64_t i = mcVoxel->x;
            const int64_t j = mcVoxel->y;
            const int64_t k = mcVoxel->z;

            // Geth the voxel bounding box
            Vector3f pMin, pMax;
            _volume->getVoxelBoundingBox(i, j, k, pMin, pMax);

            const float x = pMax.x() - 0.5 * voxelSize;       // dx * i;
            const float x_dx = pMax.x() + 0.5 * voxelSize;    // dx * (i + 1);

            const float y = pMax.y() - 0.5 * voxelSize;       // dy * j;
            const float y_dy = pMax.y() + 0.5 * voxelSize;    // dy * (j + 1);

            const float z = pMax.z() - 0.5 * voxelSize;       // dz * k;
            const float z_dz = pMax.z() + 0.5 * voxelSize;    // dz * (k + 1);

            // Get the edges from the edge table
            int edges = MC_EDGE_TABLE[mcVoxel->cubeIndex];

            // An array of the unique indices per case
            size_t uniqueIndices[12];
            for (size_t idx = 0; idx < 12; ++idx)
                uniqueIndices[idx] = 0;

            // Helper lambda to check if a vertex should be marked as border based on edge position
            // Vertices are marked as border only if the edge they're on is exactly on the volume boundary
            auto checkAndMarkBorder = [&](size_t vertexIndex, bool isEdgeOnBoundary) {
                if (isEdgeOnBoundary && vertexIndex < vertices.size())
                {
                    borderVertexIndices.insert(vertexIndex);
                }
            };

            // Generate vertices and avoid DUPLICATE vertices!
            // Note: For edges shared between cubes, we check the edge position, not just cube position
            if(edges & 0x040)
            {
                uniqueIndices[6] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 0;
                sharedIndices[idx] = uniqueIndices[6];
                addSharedVertex(x_dx, y_dy, z_dz, x, 0, v[6], v[7], _isoValue, vertices);
                // Edge 6 is between voxels at (i+1, j+1, k+1) and (i, j+1, k+1)
                // Edge is on boundary if it touches a volume face
                checkAndMarkBorder(uniqueIndices[6], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1));
            }

            if(edges & 0x020)
            {
                uniqueIndices[5] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 1;
                sharedIndices[idx] = uniqueIndices[5];
                addSharedVertex(x_dx, y, z_dz, y_dy, 1, v[5], v[6], _isoValue, vertices);
                // Edge 5 is between voxels at (i+1, j, k+1) and (i+1, j+1, k+1)
                // Edge is on boundary if it touches a volume face
                checkAndMarkBorder(uniqueIndices[5], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1));
            }

            if(edges & 0x400)
            {
                uniqueIndices[10] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 2;
                sharedIndices[idx] = uniqueIndices[10];
                addSharedVertex(x_dx, y_dy, z, z_dz, 2, v[2], v[6], _isoValue, vertices);
                // Edge 10 is between voxels at (i+1, j+1, k) and (i+1, j+1, k+1)
                // Edge is on boundary if it touches a volume face
                checkAndMarkBorder(uniqueIndices[10], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1));
            }

            if(edges & 0x001)
            {
                // Edge 0 is between voxels at positions (i, j, k) and (i+1, j, k)
                // Edge is on boundary if the edge touches a volume face:
                // - i == -1 (left boundary: edge before voxel 0) OR i == width-1 (right boundary: edge after voxel width-1)
                // - j == -1 (bottom boundary: edge before voxel 0) OR j == height-1 (top boundary: edge after voxel height-1)
                // - k == -1 (back boundary: edge before voxel 0) OR k == depth-1 (front boundary: edge after voxel depth-1)
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(j == 0 || k == 0)
                {
                    // This is the first time creating this edge
                    uniqueIndices[0] = vertices.size();
                    addSharedVertex(x, y, z, x_dx, 0, v[0], v[1], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[0], isEdgeOnBoundary);
                }
                else
                {
                    // This edge was already created by a previous cube
                    int64_t idx = i * yz3 + (j - 1) * z3 + (k - 1) * 3 + 0;
                    uniqueIndices[0] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    // Need to check if the edge between voxels at (i, j, k) and (i+1, j, k) is on boundary
                    checkAndMarkBorder(uniqueIndices[0], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1));
                }
            }

            if(edges & 0x002)
            {
                // Edge 1 is between voxels at (i+1, j, k) and (i+1, j+1, k)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(k == 0)
                {
                    uniqueIndices[1] = vertices.size();
                    addSharedVertex(x_dx, y, z, y_dy, 1, v[1], v[2], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[1], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 1;
                    uniqueIndices[1] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[1], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1));
                }
            }

            if(edges & 0x004)
            {
                // Edge 2 is between voxels at (i+1, j+1, k) and (i, j+1, k)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(k == 0)
                {
                    uniqueIndices[2] = vertices.size();
                    addSharedVertex(x_dx, y_dy, z, x, 0, v[2], v[3], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[2], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 0;
                    uniqueIndices[2] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[2], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1));
                }
            }

            if(edges & 0x008)
            {
                // Edge 3 is between voxels at (i, j+1, k) and (i, j, k)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(i == 0 || k == 0)
                {
                    uniqueIndices[3] = vertices.size();
                    addSharedVertex(x, y_dy, z, y, 1, v[3], v[0], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[3], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + (k - 1) * 3 + 1;
                    uniqueIndices[3] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[3], (i < 0 || i >= width - 1 || j < 0 || j >= height - 1));
                }
            }

            if(edges & 0x010)
            {
                // Edge 4 is between voxels at (i, j, k) and (i, j, k+1)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(j == 0)
                {
                    uniqueIndices[4] = vertices.size();
                    addSharedVertex(x, y, z_dz, x_dx, 0, v[4], v[5], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[4], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 0;
                    uniqueIndices[4] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[4], (i < 0 || i >= width - 1 || k < 0 || k >= depth - 1));
                }
            }

            if(edges & 0x080)
            {
                // Edge 7 is between voxels at (i, j+1, k+1) and (i, j, k+1)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(i == 0)
                {
                    uniqueIndices[7] = vertices.size();
                    addSharedVertex(x, y_dy, z_dz, y, 1, v[7], v[4], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[7], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 1;
                    uniqueIndices[7] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[7], (j < 0 || j >= height - 1 || k < 0 || k >= depth - 1));
                }
            }

            if(edges & 0x100)
            {
                // Edge 8 is between voxels at (i, j, k) and (i, j, k+1)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(i == 0 || j == 0)
                {
                    uniqueIndices[8] = vertices.size();
                    addSharedVertex(x, y, z, z_dz, 2, v[0], v[4], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[8], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + (j - 1) * z3 + k * 3 + 2;
                    uniqueIndices[8] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[8], (k < 0 || k >= depth - 1));
                }
            }

            if(edges & 0x200)
            {
                // Edge 9 is between voxels at (i+1, j, k) and (i+1, j, k+1)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(j == 0)
                {
                    uniqueIndices[9] = vertices.size();
                    addSharedVertex(x_dx, y, z, z_dz, 2, v[1], v[5], _isoValue,vertices);
                    checkAndMarkBorder(uniqueIndices[9], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 2;
                    uniqueIndices[9] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[9], (i < 0 || i >= width - 1 || k < 0 || k >= depth - 1));
                }
            }

            if(edges & 0x800)
            {
                // Edge 11 is between voxels at (i, j+1, k) and (i, j+1, k+1)
                // Edge is on boundary if it touches a volume face
                bool isEdgeOnBoundary = (i < 0 || i >= width - 1 || j < 0 || j >= height - 1 || k < 0 || k >= depth - 1);
                if(i == 0)
                {
                    uniqueIndices[11] = vertices.size();
                    addSharedVertex(x, y_dy, z, z_dz, 2, v[3], v[7], _isoValue, vertices);
                    checkAndMarkBorder(uniqueIndices[11], isEdgeOnBoundary);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 2;
                    uniqueIndices[11] = sharedIndices[idx];
                    // For shared edges, check if THIS edge position is on boundary
                    checkAndMarkBorder(uniqueIndices[11], (j < 0 || j >= height - 1 || k < 0 || k >= depth - 1));
                }
            }

            // Generate all faces (including boundary faces)
            // Boundary faces will be removed later using border_vertices information
            int64_t vertexIndex;
            int32_t* triangleTablePtr = MC_TRIANGLE_TABLE[mcVoxel->cubeIndex];
            for (size_t tIndex = 0;
                vertexIndex = triangleTablePtr[tIndex], vertexIndex != -1; ++tIndex)
            {
                polygons.push_back(uniqueIndices[vertexIndex]);
            }
        }
    }

    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Construct the triangles from the indices
    for (size_t i = 0; i < polygons.size() / 3; i++)
    {
        Triangle t;
        t[0] = polygons[i * 3 + 0];
        t[1] = polygons[i * 3 + 1];
        t[2] = polygons[i * 3 + 2];
        triangles.push_back(t);
    }

    // Clear
    delete [] sharedIndices;
    polygons.clear();
    polygons.shrink_to_fit();
}

bool MarchingCubes::_isOnBoundary(const int64_t i, const int64_t j, const int64_t k) const
{
    const int64_t width = static_cast<int64_t>(_volume->getWidth());
    const int64_t height = static_cast<int64_t>(_volume->getHeight());
    const int64_t depth = static_cast<int64_t>(_volume->getDepth());

    // Check if cube is on or outside the volume boundary
    // A cube at (i, j, k) samples voxels at (i, j, k) to (i+1, j+1, k+1)
    // The cube is on the boundary if it touches a volume face
    // Volume has voxels at positions 0 to width-1, so cube is on boundary if:
    // - i < 0 (left boundary) OR i >= width-1 (right boundary)
    // - j < 0 (bottom boundary) OR j >= height-1 (top boundary)
    // - k < 0 (back boundary) OR k >= depth-1 (front boundary)
    return (i < 0 || i >= width - 1) ||
           (j < 0 || j >= height - 1) ||
           (k < 0 || k >= depth - 1);
}

Mesh* MarchingCubes::generateMeshFromVolume(Volume* volume)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< MarchingCubes > workflow = std::make_unique< MarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateMesh();
}

AdvancedMesh* MarchingCubes::generateAdvancedMeshFromVolume(Volume* volume)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< MarchingCubes > workflow = std::make_unique< MarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateAdvancedMesh();
}

}
