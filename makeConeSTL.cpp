// makeConeSTL.cpp
// Generates a binary STL of a circular cone with side walls and base.
// The cone peak is at the center; Z=0 at all four edges and beyond.
// Output is suitable as test input to makeWRMToolpath.

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <string>

struct Point { float x, y, z; };
struct Triangle { Point normal; Point vertices[3]; };

static void setNormal(Triangle& tri) {
    float ux = tri.vertices[1].x - tri.vertices[0].x;
    float uy = tri.vertices[1].y - tri.vertices[0].y;
    float uz = tri.vertices[1].z - tri.vertices[0].z;
    float vx = tri.vertices[2].x - tri.vertices[0].x;
    float vy = tri.vertices[2].y - tri.vertices[0].y;
    float vz = tri.vertices[2].z - tri.vertices[0].z;
    float nx = uy*vz - uz*vy;
    float ny = uz*vx - ux*vz;
    float nz = ux*vy - uy*vx;
    float len = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (len > 0) { nx /= len; ny /= len; nz /= len; }
    tri.normal = {nx, ny, nz};
}

static void writeSTL(const std::vector<Triangle>& tris, const std::string& filename) {
    std::ofstream f(filename, std::ios::binary);
    char header[80] = {};
    f.write(header, sizeof(header));
    uint32_t n = (uint32_t)tris.size();
    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for (const auto& tri : tris) {
        f.write(reinterpret_cast<const char*>(&tri.normal), sizeof(Point));
        for (const auto& v : tri.vertices)
            f.write(reinterpret_cast<const char*>(&v), sizeof(Point));
        uint16_t attr = 0;
        f.write(reinterpret_cast<const char*>(&attr), sizeof(attr));
    }
}

static void addTri(std::vector<Triangle>& tris,
                   Point a, Point b, Point c) {
    Triangle t;
    t.vertices[0] = a; t.vertices[1] = b; t.vertices[2] = c;
    setNormal(t);
    tris.push_back(t);
}

int main() {
    const int   N       = 10;    // grid points per side (9 segments per edge)
    const float W       = 2.0f;  // model width  (inches, X axis)
    const float H       = 2.0f;  // model height (inches, Y axis)
    const float peakZ   = 0.9f;  // peak height  (inches)
    const float baseZ   = 0.3f;  // terrain height at the boundary edges (inches)
                                 // walls run from Z=0 to baseZ, making them visible

    // Terrain Z: circular cone with radius scaled to the corner distance so the
    // cone reaches baseZ at all four corners and has nonzero gradient everywhere.
    // No ridges, no flat triangles.  Edge midpoints end up slightly above baseZ.
    float R = std::sqrt(W * W / 4.0f + H * H / 4.0f);  // distance center→corner
    std::vector<std::vector<float>> z(N, std::vector<float>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float x  = j * W / (N - 1);
            float y  = i * H / (N - 1);
            float dx = x - W / 2.0f;
            float dy = y - H / 2.0f;
            float r  = std::sqrt(dx*dx + dy*dy);
            z[i][j]  = baseZ + (peakZ - baseZ) * std::max(0.0f, 1.0f - r / R);
        }
    }

    auto P = [&](int j, int i) -> Point {
        return { j * W / (N-1), i * H / (N-1), z[i][j] };
    };

    std::vector<Triangle> tris;

    // ---- Terrain surface (top) ------------------------------------------
    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N-1; ++j) {
            addTri(tris, P(j,i),   P(j+1,i), P(j,  i+1));
            addTri(tris, P(j+1,i), P(j+1,i+1), P(j,i+1));
        }
    }

    // ---- Side walls (terrain edge → Z=0) --------------------------------
    // Left wall  (X=0)
    for (int i = 0; i < N-1; ++i) {
        float y0 = i*(H/(N-1)), y1 = (i+1)*(H/(N-1));
        float z0 = z[i][0],     z1 = z[i+1][0];
        addTri(tris, {0,y0,z0}, {0,y0,0},  {0,y1,z1});
        addTri(tris, {0,y0,0},  {0,y1,0},  {0,y1,z1});
    }
    // Right wall (X=W)
    for (int i = 0; i < N-1; ++i) {
        float y0 = i*(H/(N-1)), y1 = (i+1)*(H/(N-1));
        float z0 = z[i][N-1],   z1 = z[i+1][N-1];
        addTri(tris, {W,y0,z0}, {W,y1,z1}, {W,y0,0});
        addTri(tris, {W,y0,0},  {W,y1,z1}, {W,y1,0});
    }
    // Front wall (Y=0)
    for (int j = 0; j < N-1; ++j) {
        float x0 = j*(W/(N-1)), x1 = (j+1)*(W/(N-1));
        float z0 = z[0][j],     z1 = z[0][j+1];
        addTri(tris, {x0,0,z0}, {x1,0,z1}, {x0,0,0});
        addTri(tris, {x0,0,0},  {x1,0,z1}, {x1,0,0});
    }
    // Back wall  (Y=H)
    for (int j = 0; j < N-1; ++j) {
        float x0 = j*(W/(N-1)), x1 = (j+1)*(W/(N-1));
        float z0 = z[N-1][j],   z1 = z[N-1][j+1];
        addTri(tris, {x0,H,z0}, {x0,H,0},  {x1,H,z1});
        addTri(tris, {x0,H,0},  {x1,H,0},  {x1,H,z1});
    }

    // ---- Base (Z=0) ------------------------------------------------------
    addTri(tris, {0,0,0}, {W,0,0}, {W,H,0});
    addTri(tris, {0,0,0}, {W,H,0}, {0,H,0});

    writeSTL(tris, "cone.stl");
    std::cout << "Wrote cone.stl (" << tris.size() << " triangles, "
              << N << "x" << N << " grid, edge Z=" << baseZ
              << "\", peak Z~" << peakZ << "\")" << std::endl;
    return 0;
}
