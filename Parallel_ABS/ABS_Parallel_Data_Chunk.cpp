#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

struct Point {
    double x, y;
};

double calculateAngle(const Point& v1, const Point& v2) {
    double dotProduct = v1.x * v2.x + v1.y * v2.y;
    double magnitudeV1 = std::sqrt(v1.x * v1.x + v1.y * v1.y);
    double magnitudeV2 = std::sqrt(v2.x * v2.x + v2.y * v2.y);
    double cosTheta = dotProduct / (magnitudeV1 * magnitudeV2);
    return std::acos(cosTheta);
}

std::vector<Point> angle_based_simplification(const std::vector<Point>& points, double theta, int start, int end) {
    if (points.size() < 3) {
        return points;
    }
    std::vector<Point> result;
    result.push_back(points[start]);

    for (int i = start + 1; i < end - 1; ++i) {
        Point v1 = { points[i].x - points[i - 1].x, points[i].y - points[i - 1].y };
        Point v2 = { points[i + 1].x - points[i].x, points[i + 1].y - points[i].y };
        double angle = calculateAngle(v1, v2);

        if (angle >= theta) {
            result.push_back(points[i]);
        }
    }

    result.push_back(points[end]);
    return result;
}

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
    std::string inputFile, output_angle_File;
    std::string path = "/home/vygunda/ganesh/project_M";

    std::cout << "Enter the input file name: ";
    std::cin >> inputFile;
    std::string inputFilePath = path + "/" + inputFile;

    std::cout << "Enter the Angle Simplification output file name: ";
    std::cin >> output_angle_File;
    if (output_angle_File.find(".txt") == std::string::npos) {
        output_angle_File += ".txt";
    }
    std::string outputFilePath = path + "/" + output_angle_File;

    double theta = 0.349;  // Threshold angle in radians

    std::vector<Point> points = readPointsFromFile(inputFilePath);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    std::vector<double> executionTimes(4);

    for (int numCores = 2; numCores <= 10; ++numCores) {

        auto start_angle = std::chrono::high_resolution_clock::now();

        std::vector<Point> simplifiedPoints;

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
                endIdx = numPoints - 1;
            }

            std::vector<Point> localSimplified = angle_based_simplification(points, theta, startIdx, endIdx);

#pragma omp critical
            {
                if (threadId != 0) localSimplified.erase(localSimplified.begin()); // Avoid duplicate boundary points
                simplifiedPoints.insert(simplifiedPoints.end(), localSimplified.begin(), localSimplified.end());
            }
        }

        auto end_angle = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> angle_duration = end_angle - start_angle;
        executionTimes[numCores - 1] = angle_duration.count();

        if (numCores == 10) { // Only store the points from the last core
            writePointsToFile(simplifiedPoints, outputFilePath);
        }

        std::cout << "Angle-based simplification execution time with " << numCores << " core(s): " 
                  << std::fixed << std::setprecision(9) << angle_duration.count() << " seconds" << std::endl;
    }

    std::cout << "Angle-based simplification complete. Simplified points saved to " << outputFilePath << std::endl;

    return 0;
}

