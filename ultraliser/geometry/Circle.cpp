#include "Circle.h"
#include <math/Vector4f.h>
#include <math/Quat4f.h>
#include <math/Matrix4f.h>

namespace Ultraliser
{

Circle::Circle(const Vector3f center, const float radius, const size_t numberSegments)
    : _center(center)
    , _radius(radius)
    , _numberSegments(numberSegments)
{
    // Construct the vertices and edges directly
    _generateVerticesAndEdges();
}

void Circle::_generateVerticesAndEdges()
{
    const float angleIncrement = 2 * ULTRALISER_PI / _numberSegments;
    for (size_t i = 0; i < _numberSegments; ++i)
    {
        // Construct the vertex and add it to the list
        const float angle = i * angleIncrement;
        const float x = _center.x() + _radius * cos(angle);
        const float y = _center.y() + _radius * sin(angle);
        const float z = _center.z();
        _vertices.push_back(Vector3f(x, y, z));

        // Connect vertices in circular order to form the edges
        size_t nextEdgeIndex = (i + 1) % _numberSegments;
        _edges.push_back(CircleEdge(i, nextEdgeIndex));
    }
}

const Vector3f Circle::_constructNormal() const
{
    // Construct a triangle from the center and first edge (since we know that vertices are in order)
    const auto& v1 = _center;
    const auto& v2 = _vertices[0];
    const auto& v3 = _vertices[1];

    // Construct the verctors
    const auto vec1 = (v2 - v1).normalized();
    const auto vec2 = (v3 - v1).normalized();
    return Vector3f::cross(vec1, vec2).normalized();
}

void Circle::_rotate(const Matrix4f& matrix)
{
    for (size_t i = 0; i < _vertices.size(); ++i)
    {
        Vector4f result = matrix * Vector4f(_vertices[i]);
        _vertices[i] = Vector3f(result.x(), result.y(), result.z());
    }
}

void Circle::mapToPointCloud( const std::vector< Vector3f >& pointCloud)
{

}

void Circle::getNormals(Vector3f& normal1, Vector3f& normal2) const
{
    // Construct a triangle from the center and first edge (since we know that vertices are in order)
    const auto& v1 = _center;
    const auto& v2 = _vertices[0];
    const auto& v3 = _vertices[1];

    // Construct the verctors
    const auto vec1 = (v2 - v1).normalized();
    const auto vec2 = (v3 - v1).normalized();
    normal1 = Vector3f::cross(vec1, vec2).normalized();
    normal2 = -normal1;
}

void Circle::rotateTowardsTargetPoint(const Vector3f& targetPoint)
{
    // Compute the direction, normal and rotation difference quatrenion
    const auto direction = (targetPoint - _center).normalized();
    const auto normal = Vector3f::normal(_center, _vertices[0], _vertices[1]);
    auto q = Vector3f::rotationDifference(normal, direction);

    // Construct the rotation matrix
    auto rotationMatrix = Matrix4f::rotation(q.normalized());

    // Translate all the vertices to the center at first
    for (auto& vertex : _vertices) { vertex = vertex - _center; vertex.print(); }

    // Rotate the vertices using the computed rotation matrix
    _rotate(rotationMatrix);

    // Translate back to the center
    for (auto& vertex : _vertices) { vertex = vertex + _center; }
}

void Circle::printData() const
{
    std::cout << "Vertices:" << NEW_LINE;
    for (const auto& vertex : _vertices)
    {
        std::cout << "("
                  << vertex.x() << ", " << vertex.y() << ", " << vertex.z()
                  << ")" << NEW_LINE;
    }

    std::cout << "Edges:" << NEW_LINE;
    for (const auto& edge : _edges)
    {
        std::cout << edge.v1 << " -> " << edge.v2 << NEW_LINE;
    }
}

void Circle::exportOBJ(const std::string prefix) const
{
    // Open the file
    std::string fileName = prefix + OBJ_EXTENSION;
    std::ofstream stream(fileName.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot write circle mesh file [ %s ]", fileName.c_str());
    }

    // Write vertices
    for (const auto& vertex : _vertices)
    {
        stream << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << NEW_LINE;
    }

    // Write faces
    for (const auto& edge : _edges)
    {
        stream << "l " << edge.v1 + 1 << " " << edge.v2 + 1 << "\n";
    }

    // Close the file stream
    stream.close();
}

}
