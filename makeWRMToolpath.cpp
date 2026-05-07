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
#include <queue>
#include <tuple>

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
            grid.z[ymap[v.y]][xmap[v.x]] = std::max(grid.z[ymap[v.y]][xmap[v.x]], v.z);

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
    float z_extreme = uphill ? -1e10f : 1e10f;  // running min (downhill) or max (uphill)

    const Grid& grid = *td.grid;
    float xmin = grid.xs.front(), xmax = grid.xs.back();
    float ymin = grid.ys.front(), ymax = grid.ys.back();
    float last_dx = 0.0f, last_dy = 0.0f;  // last normalized step direction

    // Stuck detection: every check_every steps, verify Z has dropped (downhill)
    // or risen (uphill) by at least min_z_drop. Allows spiraling around the summit
    // as long as the path keeps descending.
    const int   check_every = 100;
    const float min_z_drop  = step_size;
    float fc0_init, fr0_init;
    physToFrac(grid, px0, py0, fc0_init, fr0_init, true);
    float chk_z = bilerp(td.surfZ, fc0_init, fr0_init);

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

        // Gradient of the bilinear patch directly — avoids spurious circular
        // flows that arise from bilinear-interpolating a pre-computed gradient.
        int nc = (int)td.surfZ[0].size(), nr = (int)td.surfZ.size();
        int c0 = std::max(0, std::min((int)fc, nc - 2));
        int r0 = std::max(0, std::min((int)fr, nr - 2));
        float tx = fc - c0, ty = fr - r0;
        float Z00 = td.surfZ[r0][c0],   Z10 = td.surfZ[r0][c0+1];
        float Z01 = td.surfZ[r0+1][c0], Z11 = td.surfZ[r0+1][c0+1];
        float gx = ((1.0f-ty)*(Z10-Z00) + ty*(Z11-Z01)) / (grid.xs[c0+1] - grid.xs[c0]);
        float gy = ((1.0f-tx)*(Z01-Z00) + tx*(Z11-Z10)) / (grid.ys[r0+1] - grid.ys[r0]);
        float gmag = std::sqrt(gx*gx + gy*gy);
        if (inside && gmag < min_grad) {
            // Coast in last known direction through flat areas (glacial benches, etc.).
            // If we have no direction yet, give up — nothing to coast on.
            if (last_dx == 0.0f && last_dy == 0.0f) break;
        }

        float z = bilerp(td.surfZ, fc, fr);
        if (!inside) z = std::max(z, ball_radius);  // don't plunge below stock base
        if (uphill  && z < z_extreme - step_size * 5.0f) break;
        if (!uphill && z > z_extreme + step_size * 5.0f) break;
        z_extreme = uphill ? std::max(z_extreme, z) : std::min(z_extreme, z);

        // Z-progress stuck detector
        if (inside && i > 0 && i % check_every == 0) {
            if (!uphill && z > chk_z - min_z_drop) break;
            if ( uphill && z < chk_z + min_z_drop) break;
            chk_z = z;
        }

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

// Priority-Flood sink filling: raises every local depression to its outlet level
// so every interior cell has at least one lower neighbor → all D8 paths reach the boundary.
std::vector<std::vector<float>> fillSinks(const Grid& grid) {
    int nr = grid.nrows(), nc = grid.ncols();
    const float eps = 1e-4f;  // tiny slope added at each fill step

    std::vector<std::vector<float>> filled = grid.z;
    std::vector<std::vector<bool>> processed(nr, std::vector<bool>(nc, false));

    using Cell = std::tuple<float, int, int>;
    std::priority_queue<Cell, std::vector<Cell>, std::greater<Cell>> pq;

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            if (r == 0 || r == nr-1 || c == 0 || c == nc-1) {
                pq.push({filled[r][c], r, c});
                processed[r][c] = true;
            }

    while (!pq.empty()) {
        Cell top = pq.top(); pq.pop();
        float z = std::get<0>(top);
        int r   = std::get<1>(top);
        int c   = std::get<2>(top);
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                if (dr == 0 && dc == 0) continue;
                int rr = r+dr, cc = c+dc;
                if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                if (processed[rr][cc]) continue;
                processed[rr][cc] = true;
                float nz = z + eps;
                if (filled[rr][cc] > nz) nz = filled[rr][cc];
                filled[rr][cc] = nz;
                pq.push(Cell(nz, rr, cc));
            }
        }
    }
    return filled;
}

// D8 flow direction + accumulation computed on an elevation array sz.
struct FlowData {
    std::vector<std::vector<int>>   dir_r;  // D8 row increment per cell
    std::vector<std::vector<int>>   dir_c;  // D8 col increment per cell
    std::vector<std::vector<float>> accum;  // flow accumulation count
};

FlowData computeFlowData(const Grid& grid,
                          const std::vector<std::vector<float>>& sz) {
    int nr = grid.nrows(), nc = grid.ncols();
    FlowData fd;
    fd.dir_r.assign(nr, std::vector<int>(nc, 0));
    fd.dir_c.assign(nr, std::vector<int>(nc, 0));
    fd.accum.assign(nr, std::vector<float>(nc, 1.0f));

    for (int r = 0; r < nr; ++r) {
        for (int c = 0; c < nc; ++c) {
            float z = sz[r][c];
            float best = 0.0f;
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    if (dr == 0 && dc == 0) continue;
                    int rr = r+dr, cc = c+dc;
                    if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                    float dz = z - sz[rr][cc];
                    if (dz <= 0.0f) continue;
                    float ddx = grid.xs[cc] - grid.xs[c];
                    float ddy = grid.ys[rr] - grid.ys[r];
                    float slope = dz / std::sqrt(ddx*ddx + ddy*ddy);
                    if (slope > best) { best = slope; fd.dir_r[r][c] = dr; fd.dir_c[r][c] = dc; }
                }
            }
        }
    }

    std::vector<std::pair<int,int>> order;
    order.reserve(nr * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            order.push_back({r, c});
    std::sort(order.begin(), order.end(), [&](const std::pair<int,int>& a,
                                              const std::pair<int,int>& b) {
        return sz[a.first][a.second] > sz[b.first][b.second];
    });

    int oob = 0;
    for (const auto& rc : order) {
        int r = rc.first, c = rc.second;
        int dr = fd.dir_r[r][c], dc = fd.dir_c[r][c];
        if (dr == 0 && dc == 0) continue;
        int tr = r+dr, tc = c+dc;
        if (tr < 0 || tr >= nr || tc < 0 || tc >= nc) { ++oob; continue; }
        fd.accum[tr][tc] += fd.accum[r][c];
    }
    if (oob > 0) std::cout << "D8 accum: " << oob << " out-of-bounds skipped" << std::endl;
    return fd;
}

// Trace a toolpath by following D8 flow direction cell-by-cell on the offset surface.
// Cannot spiral: each step moves to a strictly lower neighbor. Continues ball_radius
// past the model boundary then stops. Stops early if occupancy blocks the path.
Toolpath traceD8Path(const Grid& grid, const TraceData& td,
                     const FlowData& flow,
                     int r0, int c0, float ball_radius,
                     OccupancyGrid* occ) {
    Toolpath path;
    int nr = grid.nrows(), nc = grid.ncols();
    float gdx = grid.xs[1] - grid.xs[0];
    float gdy = grid.ys[1] - grid.ys[0];
    float xmin = grid.xs.front(), xmax = grid.xs.back();
    float ymin = grid.ys.front(), ymax = grid.ys.back();

    int r = r0, c = c0;
    int last_dr = 0, last_dc = 0;
    float px = grid.xs[c0], py = grid.ys[r0];
    float last_z = td.surfZ[r0][c0];  // hold last on-grid Z for boundary extension

    for (int steps = 0; steps < 4000; ++steps) {
        bool inside = (r >= 0 && r < nr && c >= 0 && c < nc);

        float x, y, z;
        if (inside) {
            bool onBoundary = (r <= 1 || r >= nr-2 || c <= 1 || c >= nc-2);
            x = grid.xs[c]; y = grid.ys[r];
            z = td.surfZ[r][c];
            if (!onBoundary) {
                last_z = z;
            } else {
                // Boundary row/col: offset surface drops to ~0 at model edge.
                // Hold last interior Z so the tool doesn't plunge into the back wall.
                z = last_z;
            }
            px = x; py = y;
        } else {
            x = px; y = py;
            z = last_z;
            float ox = std::max(0.0f, std::max(xmin - x, x - xmax));
            float oy = std::max(0.0f, std::max(ymin - y, y - ymax));
            if (std::sqrt(ox*ox + oy*oy) >= ball_radius) break;
        }

        if (inside && occ && occ->isOccupied(x, y)) break;

        path.pts.push_back({x, y, z});

        if (inside) {
            int dr = flow.dir_r[r][c], dc = flow.dir_c[r][c];
            if (dr == 0 && dc == 0) break;
            last_dr = dr; last_dc = dc;
            r += dr; c += dc;
        } else {
            px += last_dc * gdx;
            py += last_dr * gdy;
            r += last_dr; c += last_dc;
        }
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

    float gdx = grid.xs[1] - grid.xs[0];
    float gdy = grid.ys[1] - grid.ys[0];

    // Compute D8 flow direction + accumulation on the offset surface.
    // The tool follows the offset surface, so flow on that surface is what matters.
    auto filledZ = fillSinks(grid);
    FlowData flow = computeFlowData(grid, filledZ);

    float max_accum = 0.0f;
    for (int r = 0; r < grid.nrows(); ++r)
        for (int c = 0; c < grid.ncols(); ++c)
            if (flow.accum[r][c] > max_accum) max_accum = flow.accum[r][c];
    std::cout << "Flow accum (terrain, sink-filled): max=" << (int)max_accum << std::endl;

    // Seed on a regular step_over grid. Each seed snaps to the nearest offset-surface
    // grid cell and traces via D8 direction — guaranteed to reach the model boundary.
    // Sort by offset-surface elevation (highest first) so summit paths run before
    // lower paths can mark their territory.
    struct Seed { float z; int r, c; };
    std::vector<Seed> seeds;

    int nx_seeds = (int)std::floor((xmax - xmin) / step_over + 1e-6f) + 1;
    int ny_seeds = (int)std::floor((ymax - ymin) / step_over + 1e-6f) + 1;
    for (int ix = 0; ix < nx_seeds; ++ix) {
        float x = xmin + ix * step_over;
        for (int iy = 0; iy < ny_seeds; ++iy) {
            float y = ymin + iy * step_over;
            int col = std::max(0, std::min((int)std::round((x - xmin) / gdx), grid.ncols() - 1));
            int row = std::max(0, std::min((int)std::round((y - ymin) / gdy), grid.nrows() - 1));
            if (grid.z[row][col] < ball_radius) continue;
            seeds.push_back({td.surfZ[row][col], row, col});
        }
    }

    std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b){ return a.z > b.z; });

    OccupancyGrid occ(xmin, ymin, xmax, ymax, step_over * 0.5f);

    for (const auto& s : seeds) {
        Toolpath path = traceD8Path(grid, td, flow, s.r, s.c, ball_radius, &occ);
        for (const auto& pt : path.pts)
            occ.mark(pt.x, pt.y);
        paths.push_back(std::move(path));
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

    const float ball_radius = 0.125f;  // 1/4" ball mill
    Surface offset = computeOffsetSurface(grid, ball_radius);
    printSurfaceBounds(offset, grid);

    const float step_over = 0.05f;     // test value; production = 0.020"
    const float step_size = 0.010f;   // inches per integration step along path
    const float feedrate  = 60.0f;    // ipm -- edit at top of output .nc file
    const bool  uphill    = false;    // false = trace downhill (rain model)

    TraceData td = buildTraceData(grid, offset);
    auto paths   = generateToolpaths(td, step_over, step_size, uphill, 1200, ball_radius);

    float safe_z = 0.0f;
    for (const auto& row : offset.pts)
        for (const auto& p : row)
            safe_z = std::max(safe_z, p.z);
    safe_z += 0.10f;

    writeGCode(paths, "output.nc", feedrate, safe_z);
    return 0;
}
