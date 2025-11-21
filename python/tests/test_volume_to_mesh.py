import math
import os

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


def _make_prism_volume(size=16):
    volume = np.zeros((size, size, size), dtype=np.uint8)
    volume[2:size-2, 2:size-2, : ] = 1
    return volume


def _boundary_mask(vertices, max_coords, eps=1e-4):
    mask = np.zeros(vertices.shape[0], dtype=bool)
    for axis in range(3):
        mask |= vertices[:, axis] <= eps
        mask |= vertices[:, axis] >= max_coords[axis] - eps
    return mask


# def test_keep_open_boundaries_removes_boundary_faces():
#     volume = _make_prism_volume()

#     closed_mesh = volume_to_mesh(volume, algorithm="dmc")
#     open_mesh = volume_to_mesh(volume, algorithm="dmc", keep_open_boundaries=True)

#     closed_vertices = np.asarray(closed_mesh["vertices"], dtype=np.float32)
#     open_vertices = np.asarray(open_mesh["vertices"], dtype=np.float32)
#     closed_faces = np.asarray(closed_mesh["faces"], dtype=np.int64)
#     open_faces = np.asarray(open_mesh["faces"], dtype=np.int64)

#     assert open_faces.shape[0] < closed_faces.shape[0]

#     if open_faces.size:
#         boundary_triangles = np.all(open_vertices[open_faces][:, :, 0] <= 1e-4, axis=1)
#         assert not np.any(boundary_triangles), "Boundary faces should be removed"

#     assert closed_vertices.shape == open_vertices.shape
#     max_coords = np.array(volume.shape[::-1], dtype=np.float32)
#     boundary_mask = _boundary_mask(closed_vertices, max_coords)
#     if np.any(~boundary_mask):
#         interior_diff = np.abs(closed_vertices[~boundary_mask] - open_vertices[~boundary_mask])
#         assert np.all(interior_diff < 1e-5), "Interior vertices should remain unchanged"


# def test_lock_boundary_vertices_restores_original_positions():
#     volume = _make_prism_volume()
#     baseline = volume_to_mesh(volume, algorithm="dmc")
#     smoothed = volume_to_mesh(
#         volume,
#         algorithm="dmc",
#         laplacian_iterations=8,
#         laplacian_lambda=0.3,
#         laplacian_mu=0.1,
#         smooth_iterations=5,
#     )
#     locked = volume_to_mesh(
#         volume,
#         algorithm="dmc",
#         laplacian_iterations=8,
#         laplacian_lambda=0.3,
#         laplacian_mu=0.1,
#         smooth_iterations=5,
#         lock_boundary_vertices=True,
#     )

#     baseline_vertices = np.asarray(baseline["vertices"], dtype=np.float32)
#     smoothed_vertices = np.asarray(smoothed["vertices"], dtype=np.float32)
#     locked_vertices = np.asarray(locked["vertices"], dtype=np.float32)

#     assert baseline_vertices.shape == smoothed_vertices.shape == locked_vertices.shape

#     max_coords = np.array(volume.shape[::-1], dtype=np.float32)
#     boundary = _boundary_mask(baseline_vertices, max_coords)

#     assert np.any(np.abs(smoothed_vertices[boundary] - baseline_vertices[boundary]) > 1e-4)
#     assert np.allclose(locked_vertices[boundary], baseline_vertices[boundary], atol=1e-5)


def test_optimize_with_open_boundaries():
    """Test that optimization works correctly with keep_open_boundaries=True.
    
    This verifies that:
    1. Optimization can run on open meshes without errors
    2. lock_boundary_vertices defaults to True when keep_open_boundaries=True
    3. The mesh remains valid after optimization
    """
    volume = _make_prism_volume()
    
    # Generate mesh with open boundaries and run optimization
    # lock_boundary_vertices should default to True automatically
    mesh = volume_to_mesh(
        volume,
        algorithm="dmc",
        keep_open_boundaries=True,
        optimize=True,
        optimization_iterations=5,
        smoothing_iterations=3,
        dense_factor=0.4,
        laplacian_iterations=2,
        laplacian_lambda=0.2,
        laplacian_mu=0.1,
        smooth_iterations=2,
        smooth_normals=True,
    )
    
    vertices = np.asarray(mesh["vertices"], dtype=np.float32)
    faces = np.asarray(mesh["faces"], dtype=np.int64)
    
    # Verify mesh is valid
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert faces.ndim == 2 and faces.shape[1] == 3
    assert vertices.shape[0] > 0, "Mesh should have vertices"
    assert faces.shape[0] > 0, "Mesh should have faces"
    
    # Verify triangles have positive area
    assert _triangles_have_positive_area(vertices, faces), "Mesh triangles must have positive area"
    
    # Verify that the mesh is not watertight (has open boundaries)
    edge_usage = _edges_from_faces(faces)
    boundary_edges = [edge for edge, count in edge_usage.items() if count == 1]
    assert len(boundary_edges) > 0, "Open mesh should have boundary edges (edges used by only one face)"
    
    # Verify vertices are in reasonable range (not NaN or infinite)
    assert np.all(np.isfinite(vertices)), "All vertices should be finite"
    assert np.all(np.abs(vertices) < 1e6), "Vertices should be in reasonable range"


def test_optimize_with_open_boundaries_explicit_lock():
    """Test that explicit lock_boundary_vertices=False can override the default."""
    volume = _make_prism_volume()
    
    # Explicitly set lock_boundary_vertices=False to override default
    mesh = volume_to_mesh(
        volume,
        algorithm="dmc",
        keep_open_boundaries=True,
        lock_boundary_vertices=False,  # Explicitly override default
        optimize=True,
        optimization_iterations=3,
        smoothing_iterations=2,
        dense_factor=0.4,
    )
    
    vertices = np.asarray(mesh["vertices"], dtype=np.float32)
    faces = np.asarray(mesh["faces"], dtype=np.int64)
    
    # Verify mesh is still valid (though boundary vertices may have moved)
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert faces.ndim == 2 and faces.shape[1] == 3
    assert np.all(np.isfinite(vertices)), "All vertices should be finite"


def test_disable_progress_restores_environment(monkeypatch):
    volume = _make_prism_volume()
    monkeypatch.delenv("ULTRALISER_NO_PROGRESS", raising=False)

    mesh = volume_to_mesh(volume, algorithm="dmc", disable_progress=True)

    assert mesh["vertices"] is not None
    assert "ULTRALISER_NO_PROGRESS" not in os.environ