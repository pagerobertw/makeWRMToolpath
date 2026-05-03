// makeWRMToolpath.cpp
// Reads a binary STL (output of makeSTL), generates a ball-mill toolpath
// following steepest gradient, outputs Fanuc-style G-code.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <map>
#include <iomanip>

struct Point {
    float x, y, z;
};

struct Triangle {
    Point normal;
    Point vertices[3];
};

bool readSTL(const std::string& filename, std::vector<Triangle>& triangles) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file: " << filename << std::endl;
        return false;
    }

    char header[80];
    file.read(header, sizeof(header));

    uint32_t numTriangles;
    file.read(reinterpret_cast<char*>(&numTriangles), sizeof(numTriangles));
    std::cout << "Triangles: " << numTriangles << std::endl;

    triangles.resize(numTriangles);
    for (auto& tri : triangles) {
        file.read(reinterpret_cast<char*>(&tri.normal),      sizeof(Point));
        file.read(reinterpret_cast<char*>(&tri.vertices[0]), sizeof(Point));
        file.read(reinterpret_cast<char*>(&tri.vertices[1]), sizeof(Point));
        file.read(reinterpret_cast<char*>(&tri.vertices[2]), sizeof(Point));
        uint16_t attr;
        file.read(reinterpret_cast<char*>(&attr), sizeof(attr));
    }

    if (!file) {
        std::cerr << "Error: file read failed or truncated." << std::endl;
        return false;
    }
    return true;
}

struct Grid {
    std::vector<float> xs;          // unique X values, sorted
    std::vector<float> ys;          // unique Y values, sorted
    std::vector<std::vector<float>> z; // z[row][col], row=Y index, col=X index
    int nrows() const { return (int)ys.size(); }
    int ncols() const { return (int)xs.size(); }
};

bool reconstructGrid(const std::vector<Triangle>& triangles, Grid& grid) {
    std::map<float, int> xmap, ymap;

    for (const auto& tri : triangles)
        for (const auto& v : tri.vertices) {
            xmap[v.x] = 0;
            ymap[v.y] = 0;
        }

    for (auto& kv : xmap) grid.xs.push_back(kv.first);
    for (auto& kv : ymap) grid.ys.push_back(kv.first);

    // index maps for fast lookup
    for (int i = 0; i < (int)grid.xs.size(); ++i) xmap[grid.xs[i]] = i;
    for (int i = 0; i < (int)grid.ys.size(); ++i) ymap[grid.ys[i]] = i;

    grid.z.assign(grid.nrows(), std::vector<float>(grid.ncols(), 0.0f));

    for (const auto& tri : triangles)
        for (const auto& v : tri.vertices)
            grid.z[ymap[v.y]][xmap[v.x]] = v.z;

    std::cout << "Grid reconstructed: " << grid.ncols() << " cols x "
              << grid.nrows() << " rows" << std::endl;
    return true;
}

// --- Offset surface -------------------------------------------------------

// Finite-difference helpers for gradient at grid point (r, c)
static float dzdxAt(const Grid& g, int r, int c) {
    if (c == 0)            return (g.z[r][1]     - g.z[r][0])     / (g.xs[1]   - g.xs[0]);
    if (c == g.ncols()-1)  return (g.z[r][c]     - g.z[r][c-1])   / (g.xs[c]   - g.xs[c-1]);
    return                        (g.z[r][c+1]   - g.z[r][c-1])   / (g.xs[c+1] - g.xs[c-1]);
}

static float dzdyAt(const Grid& g, int r, int c) {
    if (r == 0)            return (g.z[1][c]     - g.z[0][c])     / (g.ys[1]   - g.ys[0]);
    if (r == g.nrows()-1)  return (g.z[r][c]     - g.z[r-1][c])   / (g.ys[r]   - g.ys[r-1]);
    return                        (g.z[r+1][c]   - g.z[r-1][c])   / (g.ys[r+1] - g.ys[r-1]);
}

// Offset surface: each grid point shifted along its surface normal by ball_radius.
// Result is stored as a 2D array of 3D points (X and Y shift slightly on steep slopes).
struct Surface {
    int rows, cols;
    std::vector<std::vector<Point>> pts;
};

Surface computeOffsetSurface(const Grid& grid, float ball_radius) {
    int nr = grid.nrows(), nc = grid.ncols();
    Surface s;
    s.rows = nr;  s.cols = nc;
    s.pts.assign(nr, std::vector<Point>(nc));

    for (int r = 0; r < nr; ++r) {
        for (int c = 0; c < nc; ++c) {
            float gx = dzdxAt(grid, r, c);
            float gy = dzdyAt(grid, r, c);
            // Surface normal: (-dz/dx, -dz/dy, 1), normalized
            float nx = -gx, ny = -gy, nz = 1.0f;
            float len = std::sqrt(nx*nx + ny*ny + nz*nz);
            nx /= len;  ny /= len;  nz /= len;

            s.pts[r][c] = { grid.xs[c] + nx * ball_radius,
                            grid.ys[r] + ny * ball_radius,
                            grid.z[r][c] + nz * ball_radius };
        }
    }
    return s;
}

void printSurfaceBounds(const Surface& s) {
    float zmin =  std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    for (const auto& row : s.pts)
        for (const auto& p : row) {
            zmin = std::min(zmin, p.z);
            zmax = std::max(zmax, p.z);
        }
    std::cout << "Offset surface Z: " << zmin << " to " << zmax
              << "  (" << (zmax - zmin) << " in)" << std::endl;
}

// --------------------------------------------------------------------------

void printBoundingBox(const std::vector<Triangle>& triangles) {
    float xmin =  std::numeric_limits<float>::max();
    float ymin =  std::numeric_limits<float>::max();
    float zmin =  std::numeric_limits<float>::max();
    float xmax = -std::numeric_limits<float>::max();
    float ymax = -std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();

    for (const auto& tri : triangles) {
        for (const auto& v : tri.vertices) {
            xmin = std::min(xmin, v.x);  xmax = std::max(xmax, v.x);
            ymin = std::min(ymin, v.y);  ymax = std::max(ymax, v.y);
            zmin = std::min(zmin, v.z);  zmax = std::max(zmax, v.z);
        }
    }

    std::cout << "Bounding box:" << std::endl;
    std::cout << "  X: " << xmin << " to " << xmax << "  (" << (xmax-xmin) << " in)" << std::endl;
    std::cout << "  Y: " << ymin << " to " << ymax << "  (" << (ymax-ymin) << " in)" << std::endl;
    std::cout << "  Z: " << zmin << " to " << zmax << "  (" << (zmax-zmin) << " in)" << std::endl;
}

// --- Toolpath tracing -----------------------------------------------------

struct Toolpath {
    std::vector<Point> pts;
};

// 2D occupancy grid: cells sized at step_size so detection is reliable.
// A cell is marked when a completed path passes through it.
// New paths stop when they enter an already-occupied cell.
struct OccupancyGrid {
    float xmin, ymin, cell;
    int   cols, rows;
    std::vector<bool> cells;

    OccupancyGrid(float xmin, float ymin, float xmax, float ymax, float cell_size)
        : xmin(xmin), ymin(ymin), cell(cell_size) {
        cols = (int)((xmax - xmin) / cell_size) + 2;
        rows = (int)((ymax - ymin) / cell_size) + 2;
        cells.assign(rows * cols, false);
    }

    int toCol(float px) const { return (int)((px - xmin) / cell); }
    int toRow(float py) const { return (int)((py - ymin) / cell); }

    bool isOccupied(float px, float py) const {
        int c = toCol(px), r = toRow(py);
        if (c < 0 || c >= cols || r < 0 || r >= rows) return false;
        return cells[r * cols + c];
    }

    void mark(float px, float py) {
        int c = toCol(px), r = toRow(py);
        if (c >= 0 && c < cols && r >= 0 && r < rows)
            cells[r * cols + c] = true;
    }
};

// Precomputed interpolation data for the offset surface
struct TraceData {
    const Grid* grid;
    std::vector<std::vector<float>> surfZ;  // offset surface Z at each grid point
    std::vector<std::vector<float>> gx;     // dZ/dX on offset surface
    std::vector<std::vector<float>> gy;     // dZ/dY on offset surface
};

TraceData buildTraceData(const Grid& grid, const Surface& surf) {
    int nr = surf.rows, nc = surf.cols;
    TraceData td;
    td.grid = &grid;
    td.surfZ.assign(nr, std::vector<float>(nc));
    td.gx.assign(nr, std::vector<float>(nc, 0.0f));
    td.gy.assign(nr, std::vector<float>(nc, 0.0f));

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            td.surfZ[r][c] = surf.pts[r][c].z;

    for (int r = 0; r < nr; ++r) {
        for (int c = 0; c < nc; ++c) {
            if (c == 0)
                td.gx[r][c] = (td.surfZ[r][1]   - td.surfZ[r][0])   / (grid.xs[1]   - grid.xs[0]);
            else if (c == nc-1)
                td.gx[r][c] = (td.surfZ[r][c]   - td.surfZ[r][c-1]) / (grid.xs[c]   - grid.xs[c-1]);
            else
                td.gx[r][c] = (td.surfZ[r][c+1] - td.surfZ[r][c-1]) / (grid.xs[c+1] - grid.xs[c-1]);

            if (r == 0)
                td.gy[r][c] = (td.surfZ[1][c]   - td.surfZ[0][c])   / (grid.ys[1]   - grid.ys[0]);
            else if (r == nr-1)
                td.gy[r][c] = (td.surfZ[r][c]   - td.surfZ[r-1][c]) / (grid.ys[r]   - grid.ys[r-1]);
            else
                td.gy[r][c] = (td.surfZ[r+1][c] - td.surfZ[r-1][c]) / (grid.ys[r+1] - grid.ys[r-1]);
        }
    }
    return td;
}

// Map physical (px, py) to fractional grid indices (fc, fr).
// If clamp=true, positions outside the boundary are clamped to the edge
// (so the caller gets the edge Z value) and the return value indicates
// whether the position was inside (true) or outside/clamped (false).
bool physToFrac(const Grid& grid, float px, float py, float& fc, float& fr,
                bool clamp = false) {
    bool inside = !(px < grid.xs.front() || px > grid.xs.back() ||
                    py < grid.ys.front() || py > grid.ys.back());
    if (!inside && !clamp) return false;

    float cpx = std::max(grid.xs.front(), std::min(grid.xs.back(),  px));
    float cpy = std::max(grid.ys.front(), std::min(grid.ys.back(), py));

    auto xIt = std::lower_bound(grid.xs.begin(), grid.xs.end(), cpx);
    int c1 = (int)(xIt - grid.xs.begin());
    if (c1 == 0) c1 = 1;
    if (c1 >= grid.ncols()) c1 = grid.ncols() - 1;
    int c0 = c1 - 1;
    fc = c0 + (cpx - grid.xs[c0]) / (grid.xs[c1] - grid.xs[c0]);

    auto yIt = std::lower_bound(grid.ys.begin(), grid.ys.end(), cpy);
    int r1 = (int)(yIt - grid.ys.begin());
    if (r1 == 0) r1 = 1;
    if (r1 >= grid.nrows()) r1 = grid.nrows() - 1;
    int r0 = r1 - 1;
    fr = r0 + (cpy - grid.ys[r0]) / (grid.ys[r1] - grid.ys[r0]);

    return inside;
}

float bilerp(const std::vector<std::vector<float>>& f, float fc, float fr) {
    int nc = (int)f[0].size(), nr = (int)f.size();
    int c0 = std::max(0, std::min((int)fc, nc - 2));
    int r0 = std::max(0, std::min((int)fr, nr - 2));
    float tx = fc - c0, ty = fr - r0;
    return (1-ty)*((1-tx)*f[r0][c0]   + tx*f[r0][c0+1])
         +    ty *((1-tx)*f[r0+1][c0] + tx*f[r0+1][c0+1]);
}

// Trace one flow line from seed (px0, py0).
// uphill=true follows the gradient (toward peak); false follows negative gradient.
Toolpath traceFlowLine(const TraceData& td, float px0, float py0,
                        float step_size, bool uphill, int max_steps,
                        float ball_radius, OccupancyGrid* occ = nullptr) {
    Toolpath path;
    float px = px0, py = py0;
    const float min_grad = 1e-3f;
    float dir = uphill ? 1.0f : -1.0f;
    float prev_z = uphill ? -1e10f : 1e10f;

    const Grid& grid = *td.grid;
    float xmin = grid.xs.front(), xmax = grid.xs.back();
    float ymin = grid.ys.front(), ymax = grid.ys.back();
    float last_dx = 0.0f, last_dy = 0.0f;  // last normalized step direction

    for (int i = 0; i < max_steps; ++i) {
        float fc, fr;
        bool inside = physToFrac(grid, px, py, fc, fr, /*clamp=*/true);

        if (!inside) {
            float ox = std::max(0.0f, std::max(xmin - px, px - xmax));
            float oy = std::max(0.0f, std::max(ymin - py, py - ymax));
            if (std::sqrt(ox*ox + oy*oy) >= ball_radius) break;
        }

        // Stop if this cell was already covered by an earlier path
        if (inside && occ && occ->isOccupied(px, py)) break;

        float gx = bilerp(td.gx, fc, fr);
        float gy = bilerp(td.gy, fc, fr);
        float gmag = std::sqrt(gx*gx + gy*gy);
        if (inside && gmag < min_grad) break;

        float z = bilerp(td.surfZ, fc, fr);
        if (uphill  && z < prev_z - 1e-5f) break;
        if (!uphill && z > prev_z + 1e-5f) break;
        prev_z = z;

        path.pts.push_back({px, py, z});

        if (inside && gmag >= min_grad) {
            last_dx = dir * gx / gmag;
            last_dy = dir * gy / gmag;
        }
        px += last_dx * step_size;
        py += last_dy * step_size;
    }
    return path;
}

// Seed points on a regular step_over grid across the entire model surface,
// then trace each one downhill (or uphill). This "rain" strategy gives uniform
// seed coverage regardless of terrain shape.
std::vector<Toolpath> generateToolpaths(const TraceData& td,
                                         float step_over, float step_size,
                                         bool uphill, int max_steps,
                                         float ball_radius) {
    const Grid& grid = *td.grid;
    std::vector<Toolpath> paths;
    float xmin = grid.xs.front(), xmax = grid.xs.back();
    float ymin = grid.ys.front(), ymax = grid.ys.back();

    OccupancyGrid occ(xmin, ymin, xmax, ymax, step_size);

    for (float x = xmin; x <= xmax + 1e-6f; x += step_over) {
        for (float y = ymin; y <= ymax + 1e-6f; y += step_over) {
            Toolpath path = traceFlowLine(td, x, y, step_size, uphill, max_steps,
                                          ball_radius, &occ);
            for (const auto& pt : path.pts)
                occ.mark(pt.x, pt.y);
            paths.push_back(std::move(path));
        }
    }

    paths.erase(std::remove_if(paths.begin(), paths.end(),
        [](const Toolpath& p){ return p.pts.size() < 2; }), paths.end());

    std::cout << "Generated " << paths.size() << " toolpaths" << std::endl;
    return paths;
}

// --- G-code output --------------------------------------------------------

void writeGCode(const std::vector<Toolpath>& paths, const std::string& filename,
                float feedrate, float safe_z) {
    std::ofstream f(filename);
    f << std::fixed << std::setprecision(4);
    f << "( makeWRMToolpath )\n";
    f << "G90 G94\n";
    f << "F" << std::setprecision(0) << feedrate << "\n";
    f << std::setprecision(4);
    f << "G0 Z" << safe_z << "\n";

    int total_pts = 0;
    for (const auto& path : paths) {
        if (path.pts.empty()) continue;
        f << "G0 X" << path.pts[0].x << " Y" << path.pts[0].y << "\n";
        f << "G1 Z" << path.pts[0].z << "\n";
        for (size_t i = 1; i < path.pts.size(); ++i)
            f << "G1 X" << path.pts[i].x
              << " Y"   << path.pts[i].y
              << " Z"   << path.pts[i].z << "\n";
        f << "G0 Z" << safe_z << "\n";
        total_pts += (int)path.pts.size();
    }

    f << "G0 Z" << safe_z << "\n";
    f << "M30\n";
    std::cout << "Wrote " << filename << " (" << paths.size() << " paths, "
              << total_pts << " points)" << std::endl;
}

// --------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    std::string inputFile = "RainierPeakReduced.stl";
    if (argc > 1)
        inputFile = argv[1];

    std::vector<Triangle> triangles;
    if (!readSTL(inputFile, triangles))
        return 1;

    printBoundingBox(triangles);

    Grid grid;
    if (!reconstructGrid(triangles, grid))
        return 1;

    const float ball_radius = 0.125f;  // 1/4" ball mill
    Surface offset = computeOffsetSurface(grid, ball_radius);
    printSurfaceBounds(offset);

    const float step_over = 0.1f;     // inches between seed points
    const float step_size = 0.010f;   // inches per integration step along path
    const float feedrate  = 60.0f;    // ipm -- edit at top of output .nc file
    const bool  uphill    = false;    // false = trace downhill (rain model)

    TraceData td = buildTraceData(grid, offset);
    auto paths   = generateToolpaths(td, step_over, step_size, uphill, 1500, ball_radius);

    float safe_z = 0.0f;
    for (const auto& row : offset.pts)
        for (const auto& p : row)
            safe_z = std::max(safe_z, p.z);
    safe_z += 0.10f;

    writeGCode(paths, "output.nc", feedrate, safe_z);
    return 0;
}
