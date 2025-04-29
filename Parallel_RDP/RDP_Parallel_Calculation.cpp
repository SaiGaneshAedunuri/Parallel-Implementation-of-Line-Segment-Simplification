
#include <chrono>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <omp.h>

struct Point {
    double x, y;
};

// Calculate the perpendicular distance from a point to a line segment
double pointLineDistance(const Point& p, const Point& start, const Point& end) {
    double A = p.x - start.x;
    double B = p.y - start.y;
    double C = end.x - start.x;
    double D = end.y - start.y;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = (len_sq != 0) ? dot / len_sq : -1;

    double xx, yy;
    if (param < 0) {
        xx = start.x;
        yy = start.y;
    } else if (param > 1) {
        xx = end.x;
        yy = end.y;
    } else {
        xx = start.x + param * C;
        yy = start.y + param * D;
    }

    double dx = p.x - xx;
    double dy = p.y - yy;
    return std::sqrt(dx * dx + dy * dy);
}

// Parallel Recursive RDP Algorithm
std::vector<Point> rdpRecursive(const std::vector<Point>& points, double epsilon) {
    if (points.size() < 2) {
        return points;
    }

    double dmax = 0;
    size_t index = 0;
    size_t end = points.size();

    // Find the point with the maximum distance
    for (size_t i = 1; i < end - 1; ++i) {
        double d = pointLineDistance(points[i], points[0], points.back());
        if (d > dmax) {
            dmax = d;
            index = i;
        }
    }

    // If max distance is greater than epsilon, recursively simplify
    if (dmax > epsilon) {
        std::vector<Point> result1, result2;

        #pragma omp task shared(result1)
        result1 = rdpRecursive(std::vector<Point>(points.begin(), points.begin() + index + 1), epsilon);

        #pragma omp task shared(result2)
        result2 = rdpRecursive(std::vector<Point>(points.begin() + index, points.end()), epsilon);

        #pragma omp taskwait
        result1.pop_back(); // Remove duplicate point
        result1.insert(result1.end(), result2.begin(), result2.end());
        return result1;
    } else {
        return {points.front(), points.back()};
    }
}

// Wrapper function to simplify input/output and start the parallel region
std::vector<Point> rdp(const std::vector<Point>& points, double epsilon) {
    std::vector<Point> result;
    #pragma omp parallel
    {
        #pragma omp single
        result = rdpRecursive(points, epsilon);
    }
    return result;
}
// Function to read points from a file
std::vector<Point> readPointsFromFile(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream infile(filename);
    std::string line;

    if (!infile.is_open()) {
        std::cerr << "Error opening input file: " << filename << std::endl;
        return points;
    }

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        Point p;
        char comma;
        if (iss >> p.x >> comma >> p.y) {
            points.push_back(p);
        } else {
            std::cerr << "Error reading line: " << line << std::endl;
        }
    }

    infile.close();
    return points;
}

// Function to write points to a file
void writePointsToFile(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << filename << std::endl;
        return;
    }

    for (const auto& p : points) {
        file << p.x << "," << p.y << "\n";
    }

    file.close();
}

int main() {
    std::string inputFile, outputFile;
    std::cout << "Enter the input file: ";
    std::cin >> inputFile;

    std::cout << "Enter the output file: ";
    std::cin >> outputFile;

    double epsilon = 0.001;

    std::vector<Point> points = readPointsFromFile(inputFile);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    for (int numCores = 1; numCores <= 10; ++numCores) {
        omp_set_num_threads(numCores);

        auto start = std::chrono::high_resolution_clock::now();
        std::vector<Point> simplified_points = rdp(points, epsilon);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration = end - start;
        std::cout << "Execution time with " << numCores << " core(s): "
                  << std::fixed << std::setprecision(9) << duration.count() << " seconds" << std::endl;

        if (numCores == 10) {
            writePointsToFile(simplified_points, outputFile);
        }
    }

    std::cout << "RDP simplification complete. Simplified points saved to " << outputFile << std::endl;
    return 0;
}
