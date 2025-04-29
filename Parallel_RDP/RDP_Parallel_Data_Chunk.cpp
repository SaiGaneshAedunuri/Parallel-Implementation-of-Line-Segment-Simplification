#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip> // For std::fixed and std::setprecision
#include <filesystem> // For filesystem operations

namespace fs = std::filesystem;

struct Point {
    double x, y;
};

double perpendicularDistance(const Point &pt, const Point &lineStart, const Point &lineEnd) {
    double dx = lineEnd.x - lineStart.x;
    double dy = lineEnd.y - lineStart.y;

    double mag = std::sqrt(dx * dx + dy * dy);
    if (mag > 0.0) {
        dx /= mag;
        dy /= mag;
    }

    double pvx = pt.x - lineStart.x;
    double pvy = pt.y - lineStart.y;

    double pvDot = dx * pvx + dy * pvy;

    double ax = pvDot * dx;
    double ay = pvDot * dy;

    double pdx = pvx - ax;
    double pdy = pvy - ay;

    return std::sqrt(pdx * pdx + pdy * pdy);
}

void rdp(const std::vector<Point> &points, double epsilon, int start, int end, std::vector<Point> &out) {
    if (end <= start) return;

    double dmax = 0.0;
    int index = start;

    for (int i = start + 1; i < end; ++i) {
        double d = perpendicularDistance(points[i], points[start], points[end]);
        if (d > dmax) {
            index = i;
            dmax = d;
        }
    }

    if (dmax > epsilon) {
        std::vector<Point> recResults1;
        std::vector<Point> recResults2;

        rdp(points, epsilon, start, index, recResults1);
        rdp(points, epsilon, index, end, recResults2);

        out.insert(out.end(), recResults1.begin(), recResults1.end());
        out.pop_back(); // Remove the duplicate point that is added twice
        out.insert(out.end(), recResults2.begin(), recResults2.end());
    } else {
        out.push_back(points[start]);
        out.push_back(points[end]);
    }
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
    std::string inputFile, output_rdp_File;
    std::string path = "/home/vygunda/ganesh/project_M";

    std::cout << "Enter the input file name: ";
    std::cin >> inputFile;
    std::string inputFilePath = path + "/" + inputFile;

    std::cout << "Enter the RDP output file name: ";
    std::cin >> output_rdp_File;
    if (output_rdp_File.find(".txt") == std::string::npos) {
        output_rdp_File += ".txt";
    }
    std::string outputFilePath = path + "/" + output_rdp_File;

    double epsilon = 0.001;

    std::vector<Point> points = readPointsFromFile(inputFilePath);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    // Perform RDP simplification in parallel
    std::vector<Point> simplifiedPoints;
    std::vector<double> executionTimes(4);
    for (int numCores = 1; numCores <= 10; ++numCores) {
	auto start_rdp = std::chrono::high_resolution_clock::now();
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

            std::vector<Point> localSimplified;

            rdp(points, epsilon, startIdx, endIdx - 1, localSimplified);

#pragma omp critical
            {
                if (threadId != 0) localSimplified.erase(localSimplified.begin());
                simplifiedPoints.insert(simplifiedPoints.end(), localSimplified.begin(), localSimplified.end());
            }	    
        }
	auto end_rdp = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> rdp_duration = end_rdp - start_rdp;
            executionTimes[numCores - 1] = rdp_duration.count();
            std::cout << "RDP execution time with " << numCores << " core(s): "
                  << std::fixed << std::setprecision(9) << rdp_duration.count() << " seconds" << std::endl;

    }
    writePointsToFile(simplifiedPoints, outputFilePath);
    std::cout << "RDP simplification complete. Simplified points saved to " << outputFilePath << std::endl;

    return 0;
}
