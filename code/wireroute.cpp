/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <omp.h>
#include <unistd.h>

#define MAX_COST 999999999

typedef std::vector<Wire> wire_set_t;
typedef std::vector<std::vector<int>> matrix_t;

inline bool on_same_line(Point start, Point end)  {
  return start.x == end.x || start.y == end.y;
}

wire_set_t get_all_wires(Point &start, Point &end, int num_threads)
{
  std::vector<wire_set_t> ws_per_thread(num_threads);
  wire_set_t ws;
  int dx = std::abs(end.x - start.x);
  int dy = std::abs(end.y - start.y);
  int explore = dx + dx + 2 * (dx * dy);
  if (num_threads > 1)
    ws.reserve(explore * 3);

  // first check if the endpoints lie on the same line.
  if (start.x == end.x || start.y == end.y) {
    Wire w{};
    w.num_pts = 2;
    w.pts[0] = start;
    w.pts[1] = end;
    ws.push_back(w);
    return ws;
  }

  // step one: horizontal and vertical bend-one turns
  Wire h, v;
  h.num_pts = 3;
  v.num_pts = 3;
  h.pts[0] = start;
  h.pts[1] = { start.x, end.y };
  h.pts[2] = end;

  v.pts[0] = start;
  v.pts[1] = { end.x, start.y };
  v.pts[2] = end;

  ws.push_back(v);
  ws.push_back(h);

  // step two: all bend-two turns
  for (int a = start.x + 1; a < end.x; a++) {
    Wire w;
    w.num_pts = 4;
    w.pts[0] = start;
    w.pts[1] = { a, start.y };
    w.pts[2] = { a, end.y };
    w.pts[3] = end;
    ws.push_back(w);
  }

  for (int b = start.y + 1; b < end.y; b++) {
    Wire w;
    w.num_pts = 4;
    w.pts[0] = start;
    w.pts[1] = { start.x, b };
    w.pts[2] = { end.x, b };
    w.pts[3] = end;
    ws.push_back(w);
  }

  #pragma omp parallel num_threads(num_threads)
  {
    int threadid = omp_get_thread_num();
    wire_set_t &wst = ws_per_thread[threadid];
    // step three: all bend-three turns
    #pragma omp for schedule(static, 64)
    for (int j = start.x + 1; j < end.x; j++) {
      for (int k = start.y + 1; k < end.y; k++) {
      // intermediary point is (j, k)
        Wire h, v;
        h.num_pts = 5;
        v.num_pts = 5;

      // double-horizontal
        h.pts[0] = start;
        h.pts[1] = { j, start.y };
        h.pts[2] = { j, k };
        h.pts[3] = { end.x, k };
        h.pts[4] = end;

      // double-vertical
        v.pts[0] = start;
        v.pts[1] = { start.x, k };
        v.pts[2] = { j, k };
        v.pts[3] = { j, end.y };
        v.pts[4] = end;
         // horizontal-vertical or vertical-horizontal collapse into a bend-2 turn!
        wst.push_back(h);
        wst.push_back(v);
      }
    }
  }

  for (auto &wst: ws_per_thread) {
    ws.insert(ws.end(), std::make_move_iterator(wst.begin()),
                        std::make_move_iterator(wst.end()));
  }

  return ws;
}

// calculate the cost for a new wire n, ignoring a past wire o,
// given the occupancy matrix
int cost_for_path(const Wire &o, const Wire &n, const matrix_t &occupancy) {
  int cost = 0;
  for (const Point &p: n) {
    int occ = occupancy[p.y][p.x];
    cost += (occ + 1) * (occ + 1);
  }
  return cost;
}

void reroute(Wire old, Wire n, matrix_t &occupancy) {
  for (Point p: old)
  {
    occupancy[p.y][p.x]--;
  }

  for (Point p: n)
  {
    occupancy[p.y][p.x]++;
  }

  return;
}

// WITHIN WIRES SOLUTION
void solve_within_wires(
    matrix_t &occupancy,
    wire_set_t &wires,
    int dim_x, int dim_y, int num_wires,
    int num_threads, float prob,
    int iters) {

    Wire empty{};
    std::cout << "solving within wires\n";
    for (int t = 0; t < iters; t++) {
      // TIME STEP LOOP
      for (Wire &wire: wires) { // holy shit auto is a thing
        Point &start = wire.pts[0];
        Point &end = wire.pts[wire.num_pts - 1];
        if (on_same_line(start, end)) continue;
        wire_set_t all_wires;
        int min_cost;
        Wire best_path;
        reroute(wire, empty, occupancy); // unroute the normal wire
        min_cost = MAX_COST;
        best_path = wire;
        all_wires = get_all_wires(start, end, num_threads);
        static std::random_device rd;
        static std::mt19937 gen(rd());

        float p = (double(std::rand()) / (RAND_MAX));
        if (0.f <= p && p <= prob)
          best_path = all_wires[std::uniform_int_distribution<>(
                  0, all_wires.size()-1)(gen)];
        else {
          #pragma omp parallel for schedule(static) num_threads(num_threads)
          for (size_t i = 0; i < all_wires.size(); i++) {
            const Wire new_path = all_wires[i];
            int new_cost;
            new_cost = cost_for_path(wire, new_path, occupancy);
            if (new_cost < min_cost  || (new_cost == min_cost && new_path.num_pts >= wire.num_pts))
            #pragma omp critical
            {
              min_cost = new_cost;
              best_path = new_path;
            }
          }
         }

        if (best_path == wire) {
          continue;
        } else {
          reroute(empty, best_path, occupancy);
          wire = best_path;
        }
      }
    }
}


// ACROSS WIRES SOLUTION
void solve_across_wires(
    matrix_t &occupancy,
    wire_set_t &wires,
    int dim_x, int dim_y, int num_wires,
    int num_threads, float prob,
    int iters, int batch_size) {

    Wire empty{};
    std::cout << "solving across wires\n";
    static std::random_device rd;
    static std::mt19937 gen(rd());

    for (int t = 0; t < iters; t++) {
      // TIME STEP LOOP
      #pragma omp parallel for schedule(dynamic, batch_size) num_threads(num_threads)
      for (int i = 0; i < wires.size(); i++) { // holy shit auto is a thing
        Wire &wire = wires[i];
        Point &start = wire.pts[0];
        Point &end = wire.pts[wire.num_pts - 1];
        if (on_same_line(start, end)) continue;
        wire_set_t all_wires;
        int min_cost;
        Wire best_path;
        #pragma omp critical
        {
          reroute(wire, empty, occupancy); // unroute the normal wire
        }
        min_cost = MAX_COST;
        best_path = wire;
        all_wires = get_all_wires(start, end, num_threads);
        float p = (double(std::rand()) / (RAND_MAX));
        if (0.f <= p && p <= prob)
          best_path = all_wires[std::uniform_int_distribution<>(
                  0, all_wires.size()-1)(gen)];
        else {
          for (Wire &new_path: all_wires) {
            int new_cost;
            new_cost = cost_for_path(wire, new_path, occupancy);
            if (new_cost < min_cost)
            {
              min_cost = new_cost;
              best_path = new_path;
            }
          }
        }

        if (best_path == wire) {
          continue;
        } else
        #pragma omp critical
        {
          reroute(empty, best_path, occupancy);
          wire = best_path;
        }
      }
    }
}

void print_stats(const matrix_t &occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto &row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

/* This function write the output into 2 files
(1) It write occupancy grids into a file
(2) It convert wires from Wire to validate_wire_t by to_validate_format
(2) It write wires into another file
*/
void write_output(
    const std::vector<Wire> &wires, const int num_wires,
    const std::vector<std::vector<int>> &occupancy, const int dim_x,
    const int dim_y,
    std::string wires_output_file_path = "outputs/wire_output.txt",
    std::string occupancy_output_file_path = "outputs/occ_output.txt") {

  std::ofstream out_occupancy(occupancy_output_file_path, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }
  out_occupancy << dim_x << ' ' << dim_y << '\n';

  for (const auto &row : occupancy) {
    for (size_t i = 0; i < row.size(); ++i)
      out_occupancy << row[i] << (i == row.size() - 1 ? "" : " ");
    out_occupancy << '\n';
  }
  out_occupancy.close();

  std::ofstream out_wires(wires_output_file_path, std::fstream::out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n';
  out_wires << num_wires << '\n';

  for (const auto &wire : wires) {
    // NOTICE: we convert to keypoint representation here, using
    // to_validate_format which need to be defined in the bottom of this file
    validate_wire_t keypoints = wire.to_validate_format();
    for (int i = 0; i < keypoints.num_pts; ++i) {
      out_wires << keypoints.p[i].x << ' ' << keypoints.p[i].y;
      if (i < keypoints.num_pts - 1)
        out_wires << ' ';
    }
    out_wires << '\n';
  }

  out_wires.close();
}


int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1) {
    switch (opt) {
    case 'f':
      input_filename = optarg;
      break;
    case 'n':
      num_threads = atoi(optarg);
      break;
    case 'p':
      SA_prob = atof(optarg);
      break;
    case 'i':
      SA_iters = atoi(optarg);
      break;
    case 'm':
      parallel_mode = *optarg;
      break;
    case 'b':
      batch_size = atoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -f input_filename -n num_threads [-p SA_prob] [-i "
                   "SA_iters] -m parallel_mode -b batch_size\n";
      exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 ||
      (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    std::cerr << "Usage: " << argv[0]
              << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] "
                 "-m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);
//  omp_set_num_threads(num_threads);
  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x, 0)); // give each value 418 which is high
  std::cout << "Question Spec: dim_x=" << dim_x << ", dim_y=" << dim_y
            << ", number of wires=" << num_wires << '\n';

  // TODO (student code start): Read the wire information from file,
  // you may need to change this if you define the wire structure differently.
  Wire empty{};
  empty.num_pts = 0;
  for (auto &wire: wires) {
    fin >> wire.pts[0].x >> wire.pts[0].y >> wire.pts[2].x >> wire.pts[2].y;
    wire.num_pts = 3;
    if (wire.pts[0].x == wire.pts[2].x || wire.pts[0].y == wire.pts[2].y) {
      wire.pts[1] = wire.pts[2];
      wire.num_pts = 2;
    } else {
      wire.pts[1].x = wire.pts[0].x;
      wire.pts[1].y = wire.pts[2].y;
    }
    reroute(empty, wire, occupancy);
  }

  /* Initialize any additional data structures needed in the algorithm */

  // Student code end
  const double init_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - init_start)
          .count();
  std::cout << "Initialization time (sec): " << std::fixed
            << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  /* TODO (student code start): Implement the wire routing algorithm here and
    feel free to structure the algorithm into different functions.
    Don't use global variables.
    Use OpenMP to parallelize the algorithm.
  */

//  solve_sequential(occupancy, wires, dim_x, dim_y, num_wires);
  // initialize wires
  // Within wires
  if (parallel_mode == 'W') {
    solve_within_wires(occupancy, wires, dim_x, dim_y, num_wires, num_threads, SA_prob, SA_iters);
    // within wires
  } else {
    // across wires
    solve_across_wires(occupancy, wires, dim_x, dim_y, num_wires, num_threads, SA_prob, SA_iters, batch_size);
  }

  // Student code end
  // DON'T CHANGE THE FOLLOWING CODE
  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* wire to run check on wires and occupancy */
  wr_checker checker(wires, occupancy);
  checker.validate();

  /* Write wires and occupancy matrix to files */
  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y);
}

/* TODO (student): implement to_validate_format to convert Wire to
  validate_wire_t keypoint representation in order to run checker and
  write output
*/
