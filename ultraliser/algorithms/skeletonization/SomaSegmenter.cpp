#include "SomaSegmenter.h"

namespace Ultraliser
{

SomaSegmenter::SomaSegmenter(Mesh *neuronMesh, const float voxelizationVoxelsPerMicron)
    : _neuronMesh(neuronMesh)
    , _voxelizationVoxelsPerMicron(voxelizationVoxelsPerMicron)
{
    /// EMPTY CONSTRUCTOR
}

SomaSegmenter::SomaSegmenter(Volume* neuronVolume)
    : _neuronVolume(neuronVolume)
{
    /// EMPTY CONSTRUCTOR
}

void SomaSegmenter::_verifyVolume(const bool verbose)
{
    // If the _neuronVolume does not exist, create it from the _neuronMesh
    if (_neuronVolume == nullptr)
    {
        // If the _neuronMesh is null, exit
        if (_neuronMesh == nullptr)
        {
            LOG_ERROR("SomaSegmenter: The given mesh of neuron is INVALID!");
        }

        // Compute the bounding box of the mesh
        Vector3f pMinInput, pMaxInput;
        _neuronMesh->computeBoundingBox(pMinInput, pMaxInput);
        const auto& meshBoundingBox = pMaxInput - pMinInput;

        // Compute the resolution of the volume
        const auto largestDimension = meshBoundingBox.getLargestDimension();

        size_t resolution = static_cast< size_t >(_voxelizationVoxelsPerMicron* largestDimension);

        // Create the _neuronVolume from the _neuronMesh;
        _neuronVolume = new Volume(pMinInput, pMaxInput, resolution, 0.1, VOLUME_TYPE::BIT, SILENT);

        // Adaptive and conservative Voxelization
        _neuronVolume->surfaceVoxelization(_neuronMesh, SILENT, false, 1.0);
        _neuronVolume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ, SILENT);
    }
}

Mesh* SomaSegmenter::segmentSomaProxyMesh(const bool verbose)
{
    // Verify the presence of the neuron volume before we proceed
    _verifyVolume(verbose);

    // Construct a NeuronSkeletonizer
    const bool ignoreSpines = false;
    const bool useAcceleration = true;
    const bool ignoreDebugging = false;
    std::unique_ptr< NeuronSkeletonizer > skeletonizer = std::make_unique< NeuronSkeletonizer >(
                _neuronVolume, ignoreSpines, useAcceleration, ignoreDebugging);

    // Initialize the skeletonizer
    skeletonizer->initialize(verbose);

    // Skeletonize the volume to obtain the centerlines
    skeletonizer->skeletonizeVolumeToCenterLines(verbose);

    // Construct the neuron graph from the volume
    return skeletonizer->constructSomaProxyMeshFromGraph(verbose);
}

Mesh* SomaSegmenter::segmentSomaMesh(const bool verbose)
{
    // Verify the presence of the neuron volume before we proceed
    _verifyVolume(verbose);

    // Construct a NeuronSkeletonizer
    const bool ignoreSpines = false;
    const bool useAcceleration = true;
    const bool ignoreDebugging = false;
    NeuronSkeletonizer*  skeletonizer = new NeuronSkeletonizer (
                _neuronVolume, ignoreSpines, useAcceleration, ignoreDebugging);

    // Initialize the skeletonizer
    skeletonizer->initialize(verbose);

    // Skeletonize the volume to obtain the centerlines
    skeletonizer->skeletonizeVolumeToCenterLines(verbose);

    return skeletonizer->constructSomaMeshFromGraph(verbose);
}

SomaSegmenter::~SomaSegmenter()
{
    // Clear the _neuronVolume
    if (_neuronVolume != nullptr) { _neuronVolume->~Volume(); }
}

}
