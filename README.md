# Parallel-Implementation-of-Line-Segment-Simplification
This repository presents optimized parallel implementations of three popular line simplification algorithms using C++ and OpenMP, focusing on improving execution time and scalability for large datasets. The project compares Visvalingam-Whyatt (VW), Angle-Based Simplification (ABS), and Ramer-Douglas-Peucker (RDP) algorithms across sequential and Parallel (parallel calculation, and parallel data chunking) strategies.

# 🔍 Overview
Line segment simplification reduces the number of points in a polyline while preserving its essential shape and features. This process is crucial for applications such as:
Digital cartography
GPS trajectory simplification
Spatial data visualization
Network optimization

While sequential implementations are often slow for large datasets, this project explores parallel processing techniques to address performance bottlenecks using:
Parallel Calculation
Parallel Data Chunking

# 🧠 Key Features
✅ Parallelized algorithms using OpenMP for efficient multithreaded execution

📊 Benchmarking and performance analysis for datasets with 2,000 and 20,000 points

📈 Comparison of execution times across 1 to 10 cores for each approach

🔍 Accuracy and scalability evaluation of each algorithm

🧩 Modular code structure for easy testing and extension

# 💡 Insights from the Study
VW: Best choice when accuracy matters, with excellent performance in multithreaded environments.

ABS: Ideal for speed and scalability; handles large datasets efficiently via chunking.

RDP: Performs well only on small datasets; recursion makes it less efficient in parallel setups.

# 🛠 Technologies Used
C++ for core implementations

OpenMP for multithreading

Custom timing utilities for performance analysis

.txt files for input/output data handling

Results(Execution times) were shown in Graphs.

# 🚀 Getting Started
Clone the repository and compile the programs using a C++ compiler with OpenMP support (e.g., g++ -fopenmp). Provide input files with coordinate data and follow prompts to generate output and benchmark files.

g++ -fopenmp file_name.cpp -o output_file_name
./output_file_name
