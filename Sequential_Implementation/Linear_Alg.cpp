#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <algorithm>

struct Point {
    double x, y;
};

// Calculate distance from a point to a line segment
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

// RDP Algorithm
std::vector<Point> rdp(const std::vector<Point>& points, double epsilon) {
    if (points.size() < 2) {
        return points;
    }

    double dmax = 0;
    size_t index = 0;
    size_t end = points.size();

    for (size_t i = 1; i < end - 1; ++i) {
        double d = pointLineDistance(points[i], points[0], points.back());
        if (d > dmax) {
            index = i;
            dmax = d;
        }
    }

    if (dmax > epsilon) {
        std::vector<Point> recResults1(points.begin(), points.begin() + index + 1);
        recResults1 = rdp(recResults1, epsilon);
        std::vector<Point> recResults2(points.begin() + index, points.end());
        recResults2 = rdp(recResults2, epsilon);
        recResults1.pop_back();
        recResults1.insert(recResults1.end(), recResults2.begin(), recResults2.end());
        return recResults1;
    } else {
        std::vector<Point> result;
        result.push_back(points.front());
        result.push_back(points.back());
        return result;
    }
}

// Calculate the area of the triangle formed by three points
double triangleArea(const Point& p1, const Point& p2, const Point& p3) {
    return std::abs(0.5 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)));
}

// VW Algorithm
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

// Function to calculate the angle between two vectors
double calculateAngle(const Point& v1, const Point& v2) {
    double dotProduct = v1.x * v2.x + v1.y * v2.y;
    double magnitudeV1 = std::sqrt(v1.x * v1.x + v1.y * v1.y);
    double magnitudeV2 = std::sqrt(v2.x * v2.x + v2.y * v2.y);
    double cosTheta = dotProduct / (magnitudeV1 * magnitudeV2);
    return std::acos(cosTheta);
}

// Angle-based simplification algorithm
std::vector<Point> angle_based_simplification(const std::vector<Point>& points, double theta) {
    if (points.size() < 3) {
        return points;
    }

    std::vector<Point> result;
    result.push_back(points.front());

    for (size_t i = 1; i < points.size() - 1; ++i) {
        Point v1 = { points[i].x - points[i - 1].x, points[i].y - points[i - 1].y };
        Point v2 = { points[i + 1].x - points[i].x, points[i + 1].y - points[i].y };
        double angle = calculateAngle(v1, v2);

        if (angle >= theta) {
            result.push_back(points[i]);
        }
    }

    result.push_back(points.back());
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
    std::string inputFile, output_rdp_File, output_vw_File, output_angle_File;
    std::string path = "/home/vygunda/ganesh/project_M";

    std::cout << "Enter the input file: ";
    std::cin >> inputFile;
    std::string ipFilePath = path + "/" + inputFile;

    std::cout << "Enter the RDP output file: ";
    std::cin >> output_rdp_File;
    if (output_rdp_File.find(".txt") == std::string::npos) {
        output_rdp_File += ".txt";
    }
    std::string outrdpFilePath = path + "/" + output_rdp_File;

    double epsilon = 0.001;

    std::vector<Point> points = readPointsFromFile(ipFilePath);
    if (points.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    auto start_rdp = std::chrono::high_resolution_clock::now();
    std::vector<Point> simplified_rdp_points = rdp(points, epsilon);
    auto end_rdp = std::chrono::high_resolution_clock::now();
    writePointsToFile(simplified_rdp_points, outrdpFilePath);

    std::chrono::duration<double> rdp_duration = end_rdp - start_rdp;
    //std::cout << "RDP execution time: " << rdp_duration.count() << " seconds" << std::endl;
    std::cout << "RDP execution time: " << std::fixed << std::setprecision(9) << rdp_duration.count() << " seconds" << std::endl;

    std::cout << "RDP simplification complete. Simplified points saved to " << outrdpFilePath << std::endl;


    double vw_epsilon = 0.0001;
    
    std::vector<Point> points2 = readPointsFromFile(ipFilePath);
    if (points2.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    std::cout << "Enter the VW output file: ";
    std::cin >> output_vw_File;
    if (output_vw_File.find(".txt") == std::string::npos) {
        output_vw_File += ".txt";
    }
    std::string outvwFilePath = path + "/" + output_vw_File;

    auto start_vw = std::chrono::high_resolution_clock::now();
    std::vector<Point> simplified_vw_points = visvalingam_whyatt(points2, vw_epsilon);
    auto end_vw = std::chrono::high_resolution_clock::now();
    writePointsToFile(simplified_vw_points, outvwFilePath);

    std::chrono::duration<double> vw_duration = end_vw - start_vw;
    //std::cout << "VW execution time: " << vw_duration.count() << " seconds" << std::endl;
    std::cout << "RDP execution time: " << std::fixed << std::setprecision(9) << vw_duration.count() << " seconds" << std::endl;

    std::cout << "VW simplification complete. Simplified points saved to " << outvwFilePath << std::endl;

    double angle_threshold = 0.349;
    //std::cout << "Enter the angle threshold (in radians): ";
    //std::cin >> angle_threshold;

    std::vector<Point> points3 = readPointsFromFile(ipFilePath);
    if (points3.empty()) {
        std::cerr << "No points read from file. Exiting." << std::endl;
        return 1;
    }

    std::cout << "Enter the angle-based simplification output file: ";
    std::cin >> output_angle_File;
    if (output_angle_File.find(".txt") == std::string::npos) {
        output_angle_File += ".txt";
    }
    std::string outangleFilePath = path + "/" + output_angle_File;

    auto start_angle = std::chrono::high_resolution_clock::now();
    std::vector<Point> simplified_angle_points = angle_based_simplification(points3, angle_threshold);
    auto end_angle = std::chrono::high_resolution_clock::now();
    writePointsToFile(simplified_angle_points, outangleFilePath);

    std::chrono::duration<double> angle_duration = end_angle - start_angle;
    //std::cout << "Angle-based execution time: " << angle_duration.count() << " seconds" << std::endl;
    std::cout << "RDP execution time: " << std::fixed << std::setprecision(9) << angle_duration.count() << " seconds" << std::endl;

    std::cout << "Angle-based simplification complete. Simplified points saved to " << outangleFilePath << std::endl;

    return 0;
}



