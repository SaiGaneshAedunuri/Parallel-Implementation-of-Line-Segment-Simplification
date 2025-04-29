#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <iomanip>
#include <algorithm>

struct Point {
    double x, y;
};

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

// Calculate the area of the triangle formed by three points
double triangleArea(const Point& p1, const Point& p2, const Point& p3) {
    return std::abs(0.5 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)));
}

// VW Algorithm (iteratively removes points based on calculated areas)
std::vector<Point> visvalingam_whyatt(std::vector<Point>& points, double epsilon) {
    if (points.size() < 3) {
        std::cerr << "Not enough points for VW algorithm." << std::endl;
        return points;
    }

    std::vector<double> areas(points.size(), 0);
    for (size_t i = 1; i < points.size() - 1; ++i) {
        areas[i] = triangleArea(points[i - 1], points[i], points[i + 1]);
    }

    while (true) {
        if (points.size() < 3) break;

        auto min_it = std::min_element(areas.begin() + 1, areas.end() - 1);
        double min_area = *min_it;

        if (min_area >= epsilon) break;

        size_t min_index = std::distance(areas.begin(), min_it);

        points.erase(points.begin() + min_index);
        areas.erase(areas.begin() + min_index);

        if (min_index > 1 && min_index - 1 < areas.size()) {
            areas[min_index - 1] = triangleArea(points[min_index - 2], points[min_index - 1], points[min_index]);
        }
        if (min_index < points.size() - 1) {
            areas[min_index] = triangleArea(points[min_index - 1], points[min_index], points[min_index + 1]);
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

    double epsilon = 0.001;

    std::vector<Point> points = readPointsFromFile(inputFilePath);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }


    // Perform VW simplification in parallel
    std::vector<Point> simplifiedPoints;
    std::vector<double> executionTimes(4);

    for (int numCores = 1; numCores <= 10; ++numCores) {
	auto start_vw = std::chrono::high_resolution_clock::now();
        omp_set_num_threads(numCores);
        int numPoints = points.size();
        int chunkSize = numPoints / numCores;
        int remainder = numPoints % numCores;

#pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            int startIdx = threadId * chunkSize + std::min(threadId, remainder);
            int endIdx = (threadId + 1) * chunkSize + std::min(threadId + 1, remainder);

            if (threadId == numCores - 1) {
                endIdx = numPoints;
            }

            std::vector<Point> localPoints(points.begin() + startIdx, points.begin() + endIdx);
            std::vector<Point> localSimplified = visvalingam_whyatt(localPoints, epsilon);

#pragma omp critical
            {
                simplifiedPoints.insert(simplifiedPoints.end(), localSimplified.begin(), localSimplified.end());
            }
        }
	auto end_vw = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> vw_duration = end_vw - start_vw;
        executionTimes[numCores - 1] = vw_duration.count();

	std::cout << "VW execution time with " << numCores << " core(s): "
                  << std::fixed << std::setprecision(9) << vw_duration.count() << " seconds" << std::endl;
    }

    // Write only the points from the third core's output to the file
    if (simplifiedPoints.size() > 10) {
        writePointsToFile(simplifiedPoints, outputFilePath);
    }
    std::cout << "VW simplification complete. Simplified points saved to " << outputFilePath << std::endl;

    return 0;
}

