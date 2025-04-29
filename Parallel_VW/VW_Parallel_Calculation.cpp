#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip> // For std::fixed and std::setprecision
#include <filesystem> // For filesystem operations
#include <algorithm> // For std::min_element

namespace fs = std::filesystem;

struct Point {
    double x, y;
};

// Function to read points from a file (same as in the original code)
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

// Function to write points to a file (same as in the original code)
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

// Function to calculate the area of the triangle formed by three points
double triangleArea(const Point& p1, const Point& p2, const Point& p3) {
    return std::abs(0.5 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)));
}

// Parallel Visvalingam-Whyatt Algorithm
std::vector<Point> visvalingam_whyatt(std::vector<Point>& points, double epsilon, int num_cores) {
    if (points.size() < 3) {
        std::cerr << "Not enough points for VW algorithm." << std::endl;
        return points;
    }

    std::vector<double> areas(points.size(), 0);
    // Initial area calculation in parallel
    #pragma omp parallel for num_threads(num_cores)
    for (size_t i = 1; i < points.size() - 1; ++i) {
        areas[i] = triangleArea(points[i - 1], points[i], points[i + 1]);
    }

    while (true) {
        if (points.size() < 3) break;

        double min_area = *std::min_element(areas.begin() + 1, areas.end() - 1);
        if (min_area >= epsilon) break;

        size_t min_index = std::distance(areas.begin(), std::find(areas.begin() + 1, areas.end() - 1, min_area));

        points.erase(points.begin() + min_index);
        areas.erase(areas.begin() + min_index);

        // Update areas for affected triangles in parallel
        #pragma omp parallel sections num_threads(num_cores)
        {
            if (min_index > 1 && min_index - 1 < areas.size()) {
                #pragma omp task
                areas[min_index - 1] = triangleArea(points[min_index - 2], points[min_index - 1], points[min_index]);
            }
            if (min_index < points.size()) {
                #pragma omp task
                areas[min_index] = triangleArea(points[min_index - 1], points[min_index], points[min_index + 1]);
            }
        }
    }

    return points;
}


int main() {
    std::string inputFile, output_vw_File;
    std::string path = "/home/vygunda/ganesh/project_M";

    std::cout << "Enter the input file name: ";
    std::cin >> inputFile;
    std::string inputFilePath = path + "/" + inputFile;

    std::cout << "Enter the VW output file name: ";
    std::cin >> output_vw_File;
    if (output_vw_File.find(".txt") == std::string::npos) {
        output_vw_File += ".txt";
    }
    std::string outputFilePath = path + "/" + output_vw_File;

    double epsilon = 0.0001;

    std::vector<Point> points = readPointsFromFile(inputFilePath);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    std::vector<double> executionTimes(4);
    for (int cores = 1; cores <= 10; ++cores) {
	auto start_vw = std::chrono::high_resolution_clock::now();
        
	std::vector<Point> simplified_points = visvalingam_whyatt(points, epsilon, cores);
	
	auto end_vw = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> vw_duration = end_vw - start_vw;
        executionTimes[cores - 1] = vw_duration.count();

        if (cores == 10) { // Only store the points from the last core
            writePointsToFile(simplified_points, outputFilePath);
        }

        std::cout << "VW execution time with " << cores << " core(s): "
                  << std::fixed << std::setprecision(9) << vw_duration.count() << " seconds" << std::endl;
    }


    std::cout << "VW simplification complete. Simplified points saved to " << outputFilePath << std::endl;

    return 0;
}

