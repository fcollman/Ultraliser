#pragma once

#include <iostream>
#include <fstream>

#include <vector>
#include <math/Vector3f.h>
#include <math/Constants.h>
#include <math/Functions.h>
#include <common/Common.h>
#include <math/Matrix4f.h>

namespace Ultraliser
{

/**
 * @brief The Circle class
 */
class Circle {

    /**
     * @brief The CircleEdge class
     */
    class CircleEdge
    {
    public:
        size_t v1, v2;

        CircleEdge(size_t v1, size_t v2)
            : v1(v1), v2(v2)
        {
            /// EMPTY CONSTRUCTOR
        }
    };

    /**
     * @brief CircleEdges
     */
    typedef std::vector< CircleEdge > CircleEdges;

public:

    /**
     * @brief Circle
     * @param center
     * @param radius
     * @param numberSegments
     */
    Circle(const Vector3f center, const float radius, const size_t numberSegments);

    /**
     * @brief mapToPointCloud
     * @param pointCloud
     */
    void mapToPointCloud( const std::vector< Vector3f >& pointCloud);

    /**
     * @brief printData
     */
    void printData() const;

    /**
     * @brief exportOBJ
     * @param prefix
     */
    void exportOBJ(const std::string prefix) const;

    /**
     * @brief getCenter
     * @return
     */
    const Vector3f getCenter() const { return _center; }

    /**
     * @brief getRadius
     * @return
     */
    const float getRadius() const { return _radius; }

    /**
     * @brief getVertices
     * @return
     */
    const std::vector< Vector3f > getVertices() const { return _vertices; }

    /**
     * @brief getNormals
     * @param normal1
     * @param normal2
     */
    void getNormals(Vector3f& normal1, Vector3f& normal2) const;

    /**
     * @brief rotateTowardsTargetPoint
     * @param targetPoint
     */
    void rotateTowardsTargetPoint(const Vector3f& targetPoint);

private:

    /**
     * @brief _generateVerticesAndEdges
     */
    void _generateVerticesAndEdges();

    /**
     * @brief _constructNormal
     * @return
     */
    const Vector3f _constructNormal() const;

    /**
     * @brief _rotate
     * @param matrix
     */
    void _rotate(const Matrix4f& matrix);

private:

    /**
     * @brief _center
     */
    Vector3f _center;

    /**
     * @brief _radius
     */
    float _radius;

    /**
     * @brief _numberSegments
     */
    size_t  _numberSegments;

    /**
     * @brief _vertices
     */
    std::vector< Vector3f > _vertices;

    /**
     * @brief _edges
     */
    CircleEdges _edges;
};

}
