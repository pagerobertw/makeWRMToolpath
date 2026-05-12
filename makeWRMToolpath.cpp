// makeWRMToolpath.cpp
// Reads a binary STL (output of makeSTL), generates a ball-mill toolpath,
// outputs Fanuc-style LinuxCNC-compatible G-code.

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

    for (int i = 0; i < (int)grid.xs.size(); ++i) xmap[grid.xs[i]] = i;
    for (int i = 0; i < (int)grid.ys.size(); ++i) ymap[grid.ys[i]] = i;

    grid.z.assign(grid.nrows(), std::vector<float>(grid.ncols(), 0.0f));

    for (const auto& tri : triangles)
        for (const auto& v : tri.vertices)
            grid.z[ymap[v.y]][xmap[v.x]] = std::max(grid.z[ymap[v.y]][xmap[v.x]], v.z);

    std::cout << "Grid reconstructed: " << grid.ncols() << " cols x "
              << grid.nrows() << " rows" << std::endl;
    return true;
}

// --- Offset surface -------------------------------------------------------

static float dzdxAt(const Grid& g, int r, int c) {
    if (c == 0)            return (g.z[r][1]   - g.z[r][0])   / (g.xs[1]   - g.xs[0]);
    if (c == g.ncols()-1)  return (g.z[r][c]   - g.z[r][c-1]) / (g.xs[c]   - g.xs[c-1]);
    return                        (g.z[r][c+1] - g.z[r][c-1]) / (g.xs[c+1] - g.xs[c-1]);
}

static float dzdyAt(const Grid& g, int r, int c) {
    if (r == 0)            return (g.z[1][c]   - g.z[0][c])   / (g.ys[1]   - g.ys[0]);
    if (r == g.nrows()-1)  return (g.z[r][c]   - g.z[r-1][c]) / (g.ys[r]   - g.ys[r-1]);
    return                        (g.z[r+1][c] - g.z[r-1][c]) / (g.ys[r+1] - g.ys[r-1]);
}

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

void printSurfaceBounds(const Surface& s, const Grid& grid) {
    float zmin =  std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    int rmin=-1, cmin=-1, rmax=-1, cmax=-1;
    for (int r=0; r<s.rows; ++r)
        for (int c=0; c<s.cols; ++c) {
            float z = s.pts[r][c].z;
            if (z < zmin) { zmin=z; rmin=r; cmin=c; }
            if (z > zmax) { zmax=z; rmax=r; cmax=c; }
        }
    std::cout << "Offset surface Z: " << zmin << " to " << zmax
              << "  (" << (zmax - zmin) << " in)" << std::endl;
    std::cout << "  min at r=" << rmin << " c=" << cmin
              << " (x=" << grid.xs[cmin] << " y=" << grid.ys[rmin]
              << " terrain_z=" << grid.z[rmin][cmin] << ")" << std::endl;
    std::cout << "  max at r=" << rmax << " c=" << cmax
              << " (x=" << grid.xs[cmax] << " y=" << grid.ys[rmax]
              << " terrain_z=" << grid.z[rmax][cmax] << ")" << std::endl;
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

// --- Toolpath types -------------------------------------------------------

struct Toolpath {
    std::vector<Point> pts;
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

// --- Lawnmower toolpath generation ----------------------------------------
// Unidirectional passes: always X+ to X- (climb cutting).
// Each pass spans from xmax+ball_radius down to xmin-ball_radius.
// Y steps from ymin-ball_radius to ymax+ball_radius in step_over increments.
// Z at each point comes from bilinear interpolation of the offset surface;
// outside the grid the terrain is flat base so Z = ball_radius.
std::vector<Toolpath> generateLawnmower(const TraceData& td,
                                         float ball_radius, float step_over,
                                         float step_size) {
    const Grid& grid = *td.grid;
    float xmin = grid.xs.front(), xmax = grid.xs.back();
    float ymin = grid.ys.front(), ymax = grid.ys.back();
    float gdx = grid.xs[1] - grid.xs[0];
    float gdy = grid.ys[1] - grid.ys[0];
    int nr = grid.nrows(), nc = grid.ncols();

    auto surfZ = [&](float x, float y) -> float {
        if (x < xmin || x > xmax || y < ymin || y > ymax)
            return ball_radius;
        int c = std::clamp((int)((x - xmin) / gdx), 0, nc - 2);
        int r = std::clamp((int)((y - ymin) / gdy), 0, nr - 2);
        float tx = (x - grid.xs[c]) / gdx;
        float ty = (y - grid.ys[r]) / gdy;
        return (1-tx)*(1-ty)*td.surfZ[r  ][c  ]
             +    tx *(1-ty)*td.surfZ[r  ][c+1]
             + (1-tx)*   ty *td.surfZ[r+1][c  ]
             +    tx *   ty *td.surfZ[r+1][c+1];
    };

    std::vector<Toolpath> paths;
    for (float y = ymin - ball_radius;
         y <= ymax + ball_radius + 1e-5f; y += step_over) {
        Toolpath p;
        for (float x = xmax + ball_radius;
             x >= xmin - ball_radius - 1e-5f; x -= step_size)
            p.pts.push_back({x, y, surfZ(x, y)});
        if ((int)p.pts.size() >= 2) paths.push_back(std::move(p));
    }
    std::cout << "Lawnmower: " << paths.size() << " passes\n";
    return paths;
}

// --- G-code output --------------------------------------------------------
// Coordinate transforms:
//   X  : unchanged
//   Y  : y - y_max   (Y0 at back wall; workpiece interior is negative Y)
//   Z  : z - z_top   (Z0 at top of part; all cuts are negative Z)
// link_dist: max XY distance for chaining consecutive paths without a lift.
void writeGCode(const std::vector<Toolpath>& paths, const std::string& filename,
                float feedrate, float safe_z, float link_dist,
                float y_max, float z_top) {
    auto fy = [&](float y){ return y - y_max; };
    auto fz = [&](float z){ return z - z_top; };

    int n = (int)paths.size();

    std::vector<bool> chains(n, false);
    for (int i = 1; i < n; ++i) {
        if (paths[i-1].pts.empty() || paths[i].pts.empty()) continue;
        const Point& prev = paths[i-1].pts.back();
        const Point& next = paths[i].pts[0];
        float dx = next.x - prev.x, dy = next.y - prev.y;
        if (std::sqrt(dx*dx + dy*dy) <= link_dist) chains[i] = true;
    }

    std::ofstream f(filename);
    f << std::fixed << std::setprecision(4);
    f << "( makeWRMToolpath )\n";
    f << "G90 G94\n";
    f << "F" << std::setprecision(0) << feedrate << "\n";
    f << std::setprecision(4);
    f << "G0 Z" << fz(safe_z) << "\n";

    int total_pts = 0, nchained = 0;
    for (int pi = 0; pi < n; ++pi) {
        const auto& path = paths[pi];
        if (path.pts.empty()) continue;

        if (chains[pi]) {
            f << "G1 X" << path.pts[0].x << " Y" << fy(path.pts[0].y)
              << " Z"   << fz(path.pts[0].z) << "\n";
            ++nchained;
        } else {
            f << "G0 X" << path.pts[0].x << " Y" << fy(path.pts[0].y) << "\n";
            f << "G1 Z" << fz(path.pts[0].z) << "\n";
        }

        for (size_t i = 1; i < path.pts.size(); ++i)
            f << "G1 X" << path.pts[i].x
              << " Y"   << fy(path.pts[i].y)
              << " Z"   << fz(path.pts[i].z) << "\n";

        bool next_chains = (pi + 1 < n) && chains[pi + 1];
        if (!next_chains)
            f << "G0 Z" << fz(safe_z) << "\n";

        total_pts += (int)path.pts.size();
    }

    f << "G0 Z" << fz(safe_z) << "\n";
    f << "M30\n";
    std::cout << "Wrote " << filename << " (" << paths.size() << " paths, "
              << total_pts << " points, " << nchained << " chained)\n";
}

// --------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    std::cout << "build: " << __DATE__ << " " << __TIME__ << std::endl;
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

    const float ball_radius = 0.09375f;  // 3/16" dia ball mill
    const float step_over   = 0.018f;    // Y increment per pass
    const float step_size   = 0.010f;    // X resolution per point
    const float feedrate    = 60.0f;     // ipm

    Surface offset = computeOffsetSurface(grid, ball_radius);
    printSurfaceBounds(offset, grid);

    TraceData td = buildTraceData(grid, offset);
    auto paths = generateLawnmower(td, ball_radius, step_over, step_size);

    // Z0 = top of part (max offset surface Z); safe_z clears the highest point
    float z_top = 0.0f;
    for (const auto& row : offset.pts)
        for (const auto& p : row)
            z_top = std::max(z_top, p.z);
    float safe_z = z_top + 0.10f + ball_radius;

    // Y0 = back wall of workpiece
    float y_max = grid.ys.back();

    const float link_dist = 0.5f;
    writeGCode(paths, "output.nc", feedrate, safe_z, link_dist, y_max, z_top);
    return 0;
}
