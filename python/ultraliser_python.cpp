#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Ultraliser.h>
#include <algorithms/DualMarchingCubes.h>
#include <algorithms/MarchingCubes.h>
#include <data/meshes/simple/Mesh.h>
#include <data/volumes/Volume.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {
constexpr float kEpsilon = 1e-6f;

struct PostprocessConfig {
  bool optimize = false;
  std::optional<size_t> optimizationIterations;
  std::optional<int64_t> smoothingIterations;
  std::optional<float> denseFactor;
  uint32_t laplacianIterations = 0;
  float laplacianLambda = 0.2f;
  float laplacianMu = 0.1f;
  size_t surfaceSmoothIterations = 0;
  bool smoothNormals = false;
};

PostprocessConfig make_postprocess_config(
    bool optimize, std::optional<size_t> optimizationIterations,
    std::optional<int64_t> smoothingIterations,
    std::optional<float> denseFactor, uint32_t laplacianIterations,
    float laplacianLambda, float laplacianMu, size_t surfaceSmoothIterations,
    bool smoothNormals) {
  PostprocessConfig cfg;
  cfg.optimize = optimize;
  cfg.optimizationIterations = optimizationIterations;
  cfg.smoothingIterations = smoothingIterations;
  cfg.denseFactor = denseFactor;
  cfg.laplacianIterations = laplacianIterations;
  cfg.laplacianLambda = laplacianLambda;
  cfg.laplacianMu = laplacianMu;
  cfg.surfaceSmoothIterations = surfaceSmoothIterations;
  cfg.smoothNormals = smoothNormals;

  const bool hasPartialOptimization = cfg.optimizationIterations.has_value() ||
                                      cfg.smoothingIterations.has_value() ||
                                      cfg.denseFactor.has_value();

  if (hasPartialOptimization) {
    if (!(cfg.optimizationIterations && cfg.smoothingIterations &&
          cfg.denseFactor)) {
      throw py::value_error("optimization_iterations, smoothing_iterations, "
                            "and dense_factor must all be provided together");
    }
    cfg.optimize = true;
  }

  return cfg;
}

void apply_postprocess(Ultraliser::Mesh *mesh, const PostprocessConfig &cfg) {
  if (!mesh)
    return;

  if (cfg.optimize) {
    if (cfg.optimizationIterations && cfg.smoothingIterations &&
        cfg.denseFactor) {
      mesh->optimize(*cfg.optimizationIterations, *cfg.smoothingIterations,
                     *cfg.denseFactor);
    } else {
      mesh->optimizeUsingDefaultParameters();
    }
  }

  if (cfg.laplacianIterations > 0) {
    mesh->applyLaplacianSmooth(cfg.laplacianIterations, cfg.laplacianLambda,
                               cfg.laplacianMu);
  }

  if (cfg.surfaceSmoothIterations > 0)
    mesh->smoothSurface(cfg.surfaceSmoothIterations);

  if (cfg.smoothNormals)
    mesh->smoothNormals();
}

std::array<float, 3> parse_spacing(const py::object &obj) {
  std::array<float, 3> spacing{1.f, 1.f, 1.f};

  if (obj.is_none())
    return spacing;

  if (py::isinstance<py::float_>(obj) || py::isinstance<py::int_>(obj)) {
    const float value = obj.cast<float>();
    spacing.fill(value);
    return spacing;
  }

  const py::sequence seq = obj.cast<py::sequence>();
  if (seq.size() != 3)
    throw py::value_error(
        "voxel_size must be a scalar or an iterable of length 3");

  for (size_t i = 0; i < 3; ++i)
    spacing[i] = seq[i].cast<float>();

  return spacing;
}

Ultraliser::Mesh *run_mesher(Ultraliser::Volume *volume,
                             const std::string &algorithm, size_t iso_value) {
  std::string lower(algorithm.size(), '\0');
  std::transform(
      algorithm.begin(), algorithm.end(), lower.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (lower == "dmc" || lower == "dual_marching_cubes") {
    auto workflow =
        std::make_unique<Ultraliser::DualMarchingCubes>(volume, iso_value);
    return workflow->generateMesh();
  }

  if (lower == "mc" || lower == "marching_cubes") {
    auto workflow =
        std::make_unique<Ultraliser::MarchingCubes>(volume, iso_value);
    return workflow->generateMesh();
  }

  throw py::value_error("algorithm must be either 'dmc' or 'mc'");
}

bool needs_scaling(const std::array<float, 3> &spacing) {
  return (std::fabs(spacing[0] - 1.f) > kEpsilon) ||
         (std::fabs(spacing[1] - 1.f) > kEpsilon) ||
         (std::fabs(spacing[2] - 1.f) > kEpsilon);
}

py::tuple mesh_to_numpy(const Ultraliser::Mesh &mesh) {
  const auto vertex_count = static_cast<py::ssize_t>(mesh.getNumberVertices());
  const auto face_count = static_cast<py::ssize_t>(mesh.getNumberTriangles());

  py::array_t<float> vertices({vertex_count, static_cast<py::ssize_t>(3)});
  auto vertices_mut = vertices.mutable_unchecked<2>();

  const auto *vertex_data = mesh.getVertices();
  for (py::ssize_t i = 0; i < vertex_count; ++i) {
    const auto &v = vertex_data[i];
    vertices_mut(i, 0) = v.x();
    vertices_mut(i, 1) = v.y();
    vertices_mut(i, 2) = v.z();
  }

  py::array_t<uint32_t> faces({face_count, static_cast<py::ssize_t>(3)});
  auto faces_mut = faces.mutable_unchecked<2>();

  const auto *face_data = mesh.getTriangles();
  for (py::ssize_t i = 0; i < face_count; ++i) {
    const auto &f = face_data[i];
    faces_mut(i, 0) = static_cast<uint32_t>(f[0]);
    faces_mut(i, 1) = static_cast<uint32_t>(f[1]);
    faces_mut(i, 2) = static_cast<uint32_t>(f[2]);
  }

  return py::make_tuple(std::move(vertices), std::move(faces));
}

} // namespace

py::dict volume_to_mesh(py::array volume, const std::string &algorithm,
                        const py::object &voxel_size, bool solid_voxels,
                        size_t iso_value, bool optimize,
                        std::optional<size_t> optimization_iterations,
                        std::optional<int64_t> smoothing_iterations,
                        std::optional<float> dense_factor,
                        uint32_t laplacian_iterations, float laplacian_lambda,
                        float laplacian_mu, size_t smooth_iterations,
                        bool smooth_normals) {
  py::array_t<uint8_t, py::array::c_style | py::array::forcecast> volume_bytes(
      volume);
  const auto info = volume_bytes.request();

  if (info.ndim != 3)
    throw py::value_error("volume must be a 3-D array");

  const auto depth = static_cast<int64_t>(info.shape[0]);
  const auto height = static_cast<int64_t>(info.shape[1]);
  const auto width = static_cast<int64_t>(info.shape[2]);
  const auto voxel_count = static_cast<size_t>(width) *
                           static_cast<size_t>(height) *
                           static_cast<size_t>(depth);

  auto volume_ptr = std::make_unique<Ultraliser::Volume>(
      width, height, depth, Ultraliser::Vector3f::ZERO,
      Ultraliser::Vector3f::ZERO, Ultraliser::VOLUME_TYPE::UI8);

  volume_ptr->copyFromBuffer(static_cast<uint8_t *>(info.ptr), voxel_count);

  if (solid_voxels) {
    volume_ptr->solidVoxelization(
        Ultraliser::Volume::SOLID_VOXELIZATION_AXIS::XYZ);
  }

  const auto spacing = parse_spacing(voxel_size);

  std::unique_ptr<Ultraliser::Mesh> mesh(
      run_mesher(volume_ptr.get(), algorithm, iso_value));
  volume_ptr.reset();

  const PostprocessConfig postprocess_cfg = make_postprocess_config(
      optimize, optimization_iterations, smoothing_iterations, dense_factor,
      laplacian_iterations, laplacian_lambda, laplacian_mu, smooth_iterations,
      smooth_normals);
  apply_postprocess(mesh.get(), postprocess_cfg);

  if (needs_scaling(spacing)) {
    mesh->scale(spacing[0], spacing[1], spacing[2]);
  }

  auto arrays = mesh_to_numpy(*mesh);
  mesh.reset();

  py::dict result;
  result["vertices"] = arrays[0];
  result["faces"] = arrays[1];
  result["voxel_size"] = py::make_tuple(spacing[0], spacing[1], spacing[2]);

  return result;
}

PYBIND11_MODULE(_core, m) {
  m.doc() = "Python bindings for Ultraliser volume-to-mesh reconstruction";

  m.def("volume_to_mesh", &volume_to_mesh, py::arg("volume"),
        py::arg("algorithm") = "dmc", py::arg("voxel_size") = py::none(),
        py::arg("solid_voxels") = false, py::arg("iso_value") = 1,
        py::arg("optimize") = false,
        py::arg("optimization_iterations") = py::none(),
        py::arg("smoothing_iterations") = py::none(),
        py::arg("dense_factor") = py::none(),
        py::arg("laplacian_iterations") = 0, py::arg("laplacian_lambda") = 0.2f,
        py::arg("laplacian_mu") = 0.1f, py::arg("smooth_iterations") = 0,
        py::arg("smooth_normals") = false,
        R"pbdoc(
Convert a binary volume into a surface mesh using Ultraliser.

Parameters
----------
volume : numpy.ndarray
    3-D array (depth, height, width) of type bool or uint8 with values {0, 1}.
algorithm : str, optional
    Either "dmc" (default) for Dual Marching Cubes or "mc" for classic Marching Cubes.
voxel_size : float or sequence of three floats, optional
    Per-axis scaling factors applied to the reconstructed mesh vertices.
solid_voxels : bool, optional
    Apply Ultraliser solid voxelization before meshing.
iso_value : int, optional
    Threshold used by the marching cubes algorithm; defaults to 1 for binary masks.
optimize : bool, optional
    Run Ultraliser mesh optimization. Defaults to False.
optimization_iterations : int, optional
    Number of optimization iterations (requires optimize=True). Provide together with
    smoothing_iterations and dense_factor.
smoothing_iterations : int, optional
    Number of smoothing iterations used during optimization.
dense_factor : float, optional
    Dense factor used during optimization to collapse short edges. Controls the aggressiveness
    of decimation. Higher values result in coarser meshes (more vertices removed).
    Typical values range from 1.0 to 2.0.
    - 1.0: Delete vertices where all connected edges are shorter than average.
    - 2.0: More aggressive, allows deleting vertices with slightly longer edges.
    - > 2.0: Very aggressive decimation.
    - < 1.0: Less aggressive, preserves more detail.
laplacian_iterations : int, optional
    Number of Laplacian smoothing passes to apply after optimization.
laplacian_lambda : float, optional
    Lambda parameter for Laplacian smoothing (default 0.2).
laplacian_mu : float, optional
    Mu parameter for Laplacian smoothing (default 0.1).
smooth_iterations : int, optional
    Number of additional surface smoothing passes.
smooth_normals : bool, optional
    Smooth vertex normals after all mesh operations.

Returns
-------
Dict with numpy arrays:
    vertices : (N, 3) float32 array of vertex positions.
    faces : (M, 3) uint32 array of triangle indices.
    voxel_size : tuple of the applied scaling factors.
)pbdoc");
}
