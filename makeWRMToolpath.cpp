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

    // Only include upward-facing triangle vertices.
    // Use the computed (geometric) normal, not the stored STL normal — many exporters
    // write zero stored normals for flat faces such as the flat base.
    // Walls have all vertices at the same x or y, so their computed nz = 0 exactly.
    // The bottom face has nz < 0 (downward-facing).
    // Only terrain and flat-border triangles have nz > 0.
    auto computedNz = [](const Triangle& t) -> float {
        float ax = t.vertices[1].x - t.vertices[0].x;
        float ay = t.vertices[1].y - t.vertices[0].y;
        float bx = t.vertices[2].x - t.vertices[0].x;
        float by = t.vertices[2].y - t.vertices[0].y;
        return ax * by - ay * bx;  // z component of cross product (unnormalized)
    };

    for (const auto& tri : triangles) {
        if (computedNz(tri) <= 0.0f) continue;
        for (const auto& v : tri.vertices) {
            xmap[v.x] = 0;
            ymap[v.y] = 0;
        }
    }

    for (auto& kv : xmap) grid.xs.push_back(kv.first);
    for (auto& kv : ymap) grid.ys.push_back(kv.first);

    for (int i = 0; i < (int)grid.xs.size(); ++i) xmap[grid.xs[i]] = i;
    for (int i = 0; i < (int)grid.ys.size(); ++i) ymap[grid.ys[i]] = i;

    grid.z.assign(grid.nrows(), std::vector<float>(grid.ncols(), 0.0f));

    for (const auto& tri : triangles) {
        if (computedNz(tri) <= 0.0f) continue;
        for (const auto& v : tri.vertices)
            grid.z[ymap[v.y]][xmap[v.x]] = std::max(grid.z[ymap[v.y]][xmap[v.x]], v.z);
    }

    std::cout << "Grid reconstructed: " << grid.ncols() << " cols x "
              << grid.nrows() << " rows" << std::endl;
    return true;
}

// --- Offset surface -------------------------------------------------------

struct Surface {
    int rows, cols;
    std::vector<std::vector<Point>> pts;
};

// Correct ball-mill offset surface: for each grid point (xi,yi), find the
// maximum safe ball-center Z by checking all terrain points within ball_radius.
// Z_offset = max over neighbors j of: z_j + sqrt(r^2 - dx^2 - dy^2)
// This is a 2D morphological dilation with a spherical structuring element.
// The local-normal approach is incorrect in concave terrain (craters, valleys)
// because it ignores neighboring terrain that would be gouged.
Surface computeOffsetSurface(const Grid& grid, float ball_radius) {
    int nr = grid.nrows(), nc = grid.ncols();
    float gdx = grid.xs[1] - grid.xs[0];
    float gdy = grid.ys[1] - grid.ys[0];
    int dc = (int)std::ceil(ball_radius / gdx) + 1;
    int dr = (int)std::ceil(ball_radius / gdy) + 1;
    float r2 = ball_radius * ball_radius;

    Surface s;
    s.rows = nr; s.cols = nc;
    s.pts.assign(nr, std::vector<Point>(nc));

    for (int ri = 0; ri < nr; ++ri) {
        for (int ci = 0; ci < nc; ++ci) {
            float xi = grid.xs[ci], yi = grid.ys[ri];
            float z_max = -std::numeric_limits<float>::max();
            for (int rj = std::max(0, ri-dr); rj <= std::min(nr-1, ri+dr); ++rj) {
                float ddy = yi - grid.ys[rj]; ddy *= ddy;
                if (ddy > r2) continue;
                for (int cj = std::max(0, ci-dc); cj <= std::min(nc-1, ci+dc); ++cj) {
                    float ddx = xi - grid.xs[cj];
                    float d2 = ddx*ddx + ddy;
                    if (d2 > r2) continue;
                    float zb = grid.z[rj][cj] + std::sqrt(r2 - d2);
                    if (zb > z_max) z_max = zb;
                }
            }
            s.pts[ri][ci] = {xi, yi, z_max};
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
};

TraceData buildTraceData(const Grid& grid, const Surface& surf) {
    int nr = surf.rows, nc = surf.cols;
    TraceData td;
    td.grid = &grid;
    td.surfZ.assign(nr, std::vector<float>(nc));
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            td.surfZ[r][c] = surf.pts[r][c].z;
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
// Unidirectional passes at an arbitrary angle (climb cutting).
//   angle_deg: direction of each cut pass in degrees, standard math convention:
//     0   = cut toward +X (east)
//     90  = cut toward +Y (north)
//     180 = cut toward -X (west)  ← original behavior
// The stepover direction is 90° clockwise from the cut direction, which
// maintains climb cutting for any angle.
// stock_*: physical workpiece extents (all triangles, walls included).
//   Pass endpoints extend ball_radius beyond the stock edge in both axes.
//   Z lookup is clamped to the terrain grid boundary.
std::vector<Toolpath> generateLawnmower(const TraceData& td,
                                         float ball_radius, float step_over,
                                         float step_size, float angle_deg,
                                         float stock_xmin, float stock_xmax,
                                         float stock_ymin, float stock_ymax) {
    const float deg2rad = 3.14159265358979f / 180.0f;
    float theta = angle_deg * deg2rad;
    float cdx = std::cos(theta),  cdy = std::sin(theta);   // cut direction unit vector
    float sdx = std::sin(theta),  sdy = -std::cos(theta);  // step direction (90° CW from cut)

    // Project all four stock corners onto cut and step axes to find sweep extents.
    float cut_min =  std::numeric_limits<float>::max();
    float cut_max = -std::numeric_limits<float>::max();
    float step_min =  std::numeric_limits<float>::max();
    float step_max = -std::numeric_limits<float>::max();
    for (float cx : {stock_xmin, stock_xmax}) {
        for (float cy : {stock_ymin, stock_ymax}) {
            float cp = cx*cdx + cy*cdy;
            float sp = cx*sdx + cy*sdy;
            cut_min  = std::min(cut_min,  cp);  cut_max  = std::max(cut_max,  cp);
            step_min = std::min(step_min, sp);  step_max = std::max(step_max, sp);
        }
    }

    const Grid& grid = *td.grid;
    float gxmin = grid.xs.front(), gxmax = grid.xs.back();
    float gymin = grid.ys.front(), gymax = grid.ys.back();
    float gdx = grid.xs[1] - grid.xs[0];
    float gdy = grid.ys[1] - grid.ys[0];
    int nr = grid.nrows(), nc = grid.ncols();

    // Z from offset surface.  Outside the terrain grid, clamp to the nearest
    // grid boundary and hold that terrain Z.
    auto surfZ = [&](float x, float y) -> float {
        float cx = std::clamp(x, gxmin, gxmax);
        float cy = std::clamp(y, gymin, gymax);
        int c = std::clamp((int)((cx - gxmin) / gdx), 0, nc - 2);
        int r = std::clamp((int)((cy - gymin) / gdy), 0, nr - 2);
        float tx = (cx - grid.xs[c]) / gdx;
        float ty = (cy - grid.ys[r]) / gdy;
        return (1-tx)*(1-ty)*td.surfZ[r  ][c  ]
             +    tx *(1-ty)*td.surfZ[r  ][c+1]
             + (1-tx)*   ty *td.surfZ[r+1][c  ]
             +    tx *   ty *td.surfZ[r+1][c+1];
    };

    // Each pass: fixed step position s, cut parameter t sweeps the full cut range.
    // World position: (x,y) = s*(sdx,sdy) + t*(cdx,cdy)
    std::vector<Toolpath> paths;
    for (float s = step_min - ball_radius;
         s <= step_max + ball_radius + 1e-5f; s += step_over) {
        Toolpath p;
        for (float t = cut_min - ball_radius;
             t <= cut_max + ball_radius + 1e-5f; t += step_size) {
            float x = s*sdx + t*cdx;
            float y = s*sdy + t*cdy;
            p.pts.push_back({x, y, surfZ(x, y)});
        }
        if ((int)p.pts.size() >= 2) paths.push_back(std::move(p));
    }
    std::cout << "Lawnmower(" << angle_deg << " deg): " << paths.size() << " passes\n";
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

    // Skip interior points that are collinear with their neighbors.
    // Cross product of the two direction vectors is zero iff the three points
    // are collinear; use a relative tolerance to handle floating-point noise.
    auto isCollinear = [](const Point& a, const Point& b, const Point& c) -> bool {
        float d1x = b.x-a.x, d1y = b.y-a.y, d1z = b.z-a.z;
        float d2x = c.x-b.x, d2y = c.y-b.y, d2z = c.z-b.z;
        float cx = d1y*d2z - d1z*d2y;
        float cy = d1z*d2x - d1x*d2z;
        float cz = d1x*d2y - d1y*d2x;
        float cross2 = cx*cx + cy*cy + cz*cz;
        float len2 = (d1x*d1x+d1y*d1y+d1z*d1z) * (d2x*d2x+d2y*d2y+d2z*d2z);
        return cross2 < 1e-10f * len2;
    };

    int total_pts = 0, total_suppressed = 0, nchained = 0;
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

        for (size_t i = 1; i < path.pts.size(); ++i) {
            if (i < path.pts.size()-1 &&
                isCollinear(path.pts[i-1], path.pts[i], path.pts[i+1])) {
                ++total_suppressed;
                continue;
            }
            f << "G1 X" << path.pts[i].x
              << " Y"   << fy(path.pts[i].y)
              << " Z"   << fz(path.pts[i].z) << "\n";
        }

        bool next_chains = (pi + 1 < n) && chains[pi + 1];
        if (!next_chains)
            f << "G0 Z" << fz(safe_z) << "\n";

        total_pts += (int)path.pts.size();
    }

    f << "G0 Z" << fz(safe_z) << "\n";
    f << "M30\n";
    std::cout << "Wrote " << filename << " (" << paths.size() << " paths, "
              << total_pts << " points, " << total_suppressed << " suppressed, "
              << (total_pts - total_suppressed) << " emitted, "
              << nchained << " chained)\n";
}

// --------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    std::cout << "build: " << __DATE__ << " " << __TIME__ << std::endl;
    std::string inputFile = "RainierPeakReduced.stl";
    if (argc > 1)
        inputFile = argv[1];
    float angle_deg = 180.0f;
    if (argc > 2)
        angle_deg = std::stof(argv[2]);

    std::vector<Triangle> triangles;
    if (!readSTL(inputFile, triangles))
        return 1;

    printBoundingBox(triangles);

    // Stock XY bounds from ALL triangles (walls define the true workpiece edge).
    // These are used for G-code coordinate origin and pass extents.
    // The terrain grid (built after wall filtering) may be slightly smaller.
    float stock_xmin =  std::numeric_limits<float>::max();
    float stock_xmax = -std::numeric_limits<float>::max();
    float stock_ymin =  std::numeric_limits<float>::max();
    float stock_ymax = -std::numeric_limits<float>::max();
    for (const auto& tri : triangles)
        for (const auto& v : tri.vertices) {
            stock_xmin = std::min(stock_xmin, v.x);  stock_xmax = std::max(stock_xmax, v.x);
            stock_ymin = std::min(stock_ymin, v.y);  stock_ymax = std::max(stock_ymax, v.y);
        }
    std::cout << "Stock bounds: X " << stock_xmin << " to " << stock_xmax
              << "  Y " << stock_ymin << " to " << stock_ymax << std::endl;

    Grid grid;
    if (!reconstructGrid(triangles, grid))
        return 1;

    const float ball_radius   = 0.09375f;  // 3/16" dia ball mill
    const float step_over     = 0.018f;    // stepover between passes
    const float step_size     = 0.010f;    // sample spacing along each pass
    const float feedrate      = 60.0f;     // ipm

    Surface offset = computeOffsetSurface(grid, ball_radius);
    printSurfaceBounds(offset, grid);

    TraceData td = buildTraceData(grid, offset);
    auto paths = generateLawnmower(td, ball_radius, step_over, step_size, angle_deg,
                                   stock_xmin, stock_xmax, stock_ymin, stock_ymax);

    // Z0 = top of part (max offset surface Z); safe_z clears the highest point
    float z_top = 0.0f;
    for (const auto& row : offset.pts)
        for (const auto& p : row)
            z_top = std::max(z_top, p.z);
    float safe_z = z_top + 0.10f + ball_radius;

    // Y0 = back wall of workpiece (stock_ymax from full bounding box, not terrain grid edge)
    float y_max = stock_ymax;

    const float link_dist = 0.5f;
    writeGCode(paths, "output.nc", feedrate, safe_z, link_dist, y_max, z_top);
    return 0;
}
