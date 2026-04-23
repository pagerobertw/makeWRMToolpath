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

    return 0;
}
