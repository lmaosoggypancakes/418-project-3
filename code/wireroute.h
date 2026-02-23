/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <cstdint>
#include <omp.h>
#include <vector>

#define MAX_PTS_PER_WIRE 5
#define COST_REPORT_DEPTH 10

/** README(student):
 We provide a way to validate consistency between wire layout
 and occupancy. Within the program, a wire checker is provided
 by:
 wr_checker Checker(wires, occupancy);
 and its validate() method can be called to validate the consistency.

 The struct below is the standard format for wires used by the wire checker.
 It contains a buffer that holds up to MAX_PTS_PER_WIRE points, and a num_pts
 field that specifies the number of points. Regardless of what representation you use for your wires, you should
 implement the Wire::to_validate_format method to convert your Wire
 to a validate_wire_t if you wish to use the checker.
*/

/* validate_wire_t is a format to represent wires by key points

For example a wire with 3 bend has 5 key poinst:
start, bend1, bend2, bend3, end.
Notice this implies two consecutive keypoints share at least same x or same y

num_pts is number of keypoints
p is an array store the xy coordinates of each keypoints

Requires (see validate.cpp validate_wire_t::cleanup() for how is this
checked）： (1) two consecutive keypoints share at least same x or same y, (2) 2
<= num_pts <= MAX_PTS_PER_WIRE.
*/
struct validate_wire_t {
  uint8_t num_pts;
  struct {
    uint16_t x;
    uint16_t y;
  } p[MAX_PTS_PER_WIRE];
  validate_wire_t &cleanup(void);
  void print_wire(void) const;
};

/* TODO (student): Define the data structure for wire here.
Feel free to redefine this, notice you also need to change how bend is read in
inside main of wireroute.cpp.

A better readable way of define this is:
int start_x, start_y, end_x, end_y, bend1_x, bend1_y, bend2_x, bend2_y, bend3_x,
bend3_y; but this might not be the most efficient way to define the solution
space for a wire with <= 3 bends.
*/
struct Point {
  int x;
  int y;
  bool operator==(const Point& o) const {
    return x == o.x && y == o.y;
  }
};

static inline int sgn (int v) {
  if (v > 0) return 1;
  if (v == 0) return 0;
  return -1;
}
struct Wire {
    int num_pts;                     // number of keypoints actually used
    Point pts[MAX_PTS_PER_WIRE];     // pts[0]=start, pts[num_pts-1]=end

    validate_wire_t to_validate_format() const {
        validate_wire_t v{};
        v.num_pts = num_pts;
        for (int i = 0; i < num_pts; i++) {
            v.p[i].x = pts[i].x;
            v.p[i].y = pts[i].y;
        }
        return v;
    }

    struct Iterator {
        const Wire* w;
        int seg;        // current segment index
        Point cur;      // current lattice point
        int dx, dy;     // step direction
        bool done;

        Iterator(const Wire* wire, bool is_end)
            : w(wire), seg(0), dx(0), dy(0), done(is_end)
        {
            if (!done && w->num_pts > 0) {
                cur = w->pts[0];
                load_dir();
            } else {
                done = true;
            }
        }

        void load_dir() {
            if (seg >= w->num_pts - 1) {
                done = true;
                return;
            }
            Point nxt = w->pts[seg + 1];
            dx = sgn(nxt.x - cur.x);
            dy = sgn(nxt.y - cur.y);
        }

        Point operator*() const { return cur; }

        Iterator& operator++() {
            if (done) return *this;

            Point target = w->pts[seg + 1];

            // If reached end of this segment, move to next
            if (cur == target) {
                seg++;
                if (seg >= w->num_pts - 1) {
                    done = true;
                    return *this;
                }
                load_dir();
                target = w->pts[seg + 1];
            }

            cur.x += dx;
            cur.y += dy;
            return *this;
        }

        bool operator!=(const Iterator& other) const {
            return done != other.done;
        }
    };

    Iterator begin() const { return Iterator(this, false); }
    Iterator end()   const { return Iterator(this, true); }
};


// Definition of the wire checker
struct wr_checker {
  std::vector<Wire> wires;
  std::vector<std::vector<int>> occupancies;
  const int nwires;
  const int dim_x;
  const int dim_y;
  wr_checker(std::vector<Wire> &wires,
             std::vector<std::vector<int>> &occupancies)
      : wires(wires), occupancies(occupancies), nwires(wires.size()),
        dim_x(occupancies[0].size()), dim_y(occupancies.size()) {}
  void validate() const;
};

const char *get_option_string(const char *option_name,
                              const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
