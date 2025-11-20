import math

import numpy as np
import pytest

from ultraliser import volume_to_mesh


@pytest.fixture(scope="module")
def binary_sphere():
    size = 64
    radius = 20.0
    center = (size - 1) / 2.0

    z, y, x = np.ogrid[:size, :size, :size]
    distance_sq = (x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2
    volume = (distance_sq <= radius ** 2).astype(np.uint8)

    return volume, radius


def _edges_from_faces(faces):
    edges = {}
    for tri in faces:
        for i in range(3):
            a = int(tri[i])
            b = int(tri[(i + 1) % 3])
            edge = (a, b) if a < b else (b, a)
            edges[edge] = edges.get(edge, 0) + 1
    return edges


def _mesh_is_watertight(faces):
    edge_usage = _edges_from_faces(faces)
    return all(count == 2 for count in edge_usage.values())


def _triangles_have_positive_area(vertices, faces, min_area=1e-8):
    tri_vertices = vertices[faces]
    cross = np.cross(tri_vertices[:, 1] - tri_vertices[:, 0],
                     tri_vertices[:, 2] - tri_vertices[:, 0])
    areas = np.linalg.norm(cross, axis=1) * 0.5
    return np.all(areas > min_area)


def _sphere_shape_is_preserved(vertices, expected_radius, tol_ratio=0.1, tol_std=0.8):
    centroid = vertices.mean(axis=0)
    distances = np.linalg.norm(vertices - centroid, axis=1)
    mean_radius = distances.mean()
    std_radius = distances.std()

    return (
        math.isclose(mean_radius, expected_radius, rel_tol=tol_ratio, abs_tol=0.5)
        and std_radius < tol_std
    )


def test_volume_to_mesh_generates_watertight_sphere(binary_sphere):
    volume, radius = binary_sphere

    mesh_data = volume_to_mesh(
        volume,
        algorithm="dmc",
        solid_voxels=False,
        iso_value=1,
        optimize=True,
        optimization_iterations=10,
        smoothing_iterations=5,
        dense_factor=0.4,
        laplacian_iterations=2,
        laplacian_lambda=0.2,
        laplacian_mu=0.1,
        smooth_iterations=3,
        smooth_normals=True,
    )

    vertices = np.asarray(mesh_data["vertices"], dtype=np.float32)
    faces = np.asarray(mesh_data["faces"], dtype=np.int64)

    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert faces.ndim == 2 and faces.shape[1] == 3

    assert _mesh_is_watertight(faces), "Mesh should be watertight with closed surfaces"
    assert _triangles_have_positive_area(vertices, faces), "Mesh triangles must have positive area"
    assert _sphere_shape_is_preserved(vertices, radius), "Mesh should approximate the input sphere"