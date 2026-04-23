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

    return 0;
}
