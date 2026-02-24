#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <utility>
#include <vector>

// Platform-specific virtual memory support
#if defined(_WIN32)
#if !defined(WIN32_LEAN_AND_MEAN)
#define WIN32_LEAN_AND_MEAN
#define MAPBOX_EARCUT_UNDEF_WIN32_LEAN_AND_MEAN
#endif
#if !defined(NOMINMAX)
#define NOMINMAX
#define MAPBOX_EARCUT_UNDEF_NOMINMAX
#endif
#include <windows.h>
#ifdef MAPBOX_EARCUT_UNDEF_WIN32_LEAN_AND_MEAN
#undef WIN32_LEAN_AND_MEAN
#undef MAPBOX_EARCUT_UNDEF_WIN32_LEAN_AND_MEAN
#endif
#ifdef MAPBOX_EARCUT_UNDEF_NOMINMAX
#undef NOMINMAX
#undef MAPBOX_EARCUT_UNDEF_NOMINMAX
#endif
#elif defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#endif

namespace mapbox {

namespace util {

template <std::size_t I, typename T>
struct nth {
    inline static typename std::tuple_element<I, T>::type get(const T& t) { return std::get<I>(t); };
};

} // namespace util

namespace detail {

template <typename N = uint32_t>
class Earcut {
public:
    // Index type used for all node references. 0 = null sentinel.
    using NodeIndex = uint32_t;

    std::vector<N> indices;
    std::size_t vertices = 0;

    template <typename Polygon>
    void operator()(const Polygon& points);

private:
    struct Node {
        Node() : x(0), y(0), i(0), steiner(0) {}
        Node(N index, double x_, double y_) : x(x_), y(y_), i(index), steiner(0) {}

        double x;
        double y;

        // previous and next vertice nodes in a polygon ring (indices, 0 = null)
        NodeIndex prev = 0;
        NodeIndex next = 0;

        // z-order curve value
        int32_t z = 0;

        // original index in polygon
        N i : (sizeof(N) * 8 - 1);

        // indicates whether this is a steiner point
        N steiner : 1;

        // previous and next nodes in z-order (indices, 0 = null)
        NodeIndex prevZ = 0;
        NodeIndex nextZ = 0;
    };

    // Cache-optimized Triangle structure for repeated geometric tests
    struct Triangle {
        const double ax, ay;
        const double bx, by;
        const double cx, cy;

        Triangle(const Node& a, const Node& b, const Node& c)
            : ax(a.x), ay(a.y), bx(b.x), by(b.y), cx(c.x), cy(c.y) {}

        inline double area() const { return (by - ay) * (cx - bx) - (bx - ax) * (cy - by); }

        inline bool containsPoint(double px, double py) const {
            return (cx - px) * (ay - py) >= (ax - px) * (cy - py) && (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
                   (bx - px) * (cy - py) >= (cx - px) * (by - py);
        }
    };

    template <typename Ring>
    NodeIndex linkedList(const Ring& points, const bool clockwise);
    NodeIndex filterPoints(NodeIndex start, NodeIndex end = 0);
    void earcutLinked(NodeIndex ear, int pass = 0);
    bool isEar(NodeIndex ear);
    bool isEarHashed(NodeIndex ear);
    NodeIndex cureLocalIntersections(NodeIndex start);
    void splitEarcut(NodeIndex start);
    template <typename Polygon>
    NodeIndex eliminateHoles(const Polygon& points, NodeIndex outerNode);
    NodeIndex eliminateHole(NodeIndex hole, NodeIndex outerNode);
    NodeIndex findHoleBridge(NodeIndex hole, NodeIndex outerNode);
    bool sectorContainsSector(NodeIndex m, NodeIndex p);
    void indexCurve(NodeIndex start);
    NodeIndex sortLinked(NodeIndex list);
    int32_t zOrder(const double x_, const double y_);
    NodeIndex getLeftmost(NodeIndex start);
    bool pointInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) const;
    bool isValidDiagonal(NodeIndex a, NodeIndex b);
    double area(NodeIndex p, NodeIndex q, NodeIndex r) const;
    bool equals(NodeIndex p1, NodeIndex p2) const;
    bool intersects(NodeIndex p1, NodeIndex q1, NodeIndex p2, NodeIndex q2) const;
    bool onSegment(NodeIndex p, NodeIndex q, NodeIndex r) const;
    int sign(double val) const;
    bool intersectsPolygon(NodeIndex a, NodeIndex b) const;
    bool locallyInside(NodeIndex a, NodeIndex b) const;
    bool middleInside(NodeIndex a, NodeIndex b) const;
    NodeIndex splitPolygon(NodeIndex a, NodeIndex b);
    template <typename Point>
    NodeIndex insertNode(std::size_t i, const Point& p, NodeIndex last);
    void removeNode(NodeIndex p);

    bool hashing;
    double minX, maxX;
    double minY, maxY;
    double inv_size = 0;

    // Arena allocator for nodes using virtual memory with demand paging.
    // Reserves a large virtual address range upfront. On POSIX, pages are
    // automatically backed by physical memory on first write. On Windows,
    // pages are explicitly committed as the arena grows.
    // Index 0 is a null sentinel; real nodes start at index 1.
    class NodeArena {
        // Maximum number of nodes the arena can hold.
        // 4M nodes * 40 bytes = ~160MB of virtual address space.
        static constexpr std::size_t kMaxNodes = 1u << 22;
        static constexpr std::size_t kPageSize = 4096;

    public:
        NodeArena() { initReserve(); }

        ~NodeArena() { releaseAll(); }

        NodeArena(const NodeArena&) = delete;
        NodeArena& operator=(const NodeArena&) = delete;
        NodeArena(NodeArena&&) = delete;
        NodeArena& operator=(NodeArena&&) = delete;

        template <typename... Args>
        NodeIndex construct(Args&&... args) {
            NodeIndex idx;
            if (freeHead_ != 0) {
                idx = freeHead_;
                freeHead_ = data_[freeHead_].next;
            } else {
                idx = static_cast<NodeIndex>(size_);
                ensureCommitted(size_ + 1);
                size_++;
            }
            new (&data_[idx]) Node(std::forward<Args>(args)...);
            return idx;
        }

        void freeNode(NodeIndex idx) {
            assert(idx != 0);
            data_[idx].next = freeHead_;
            freeHead_ = idx;
        }

        void reset() {
            size_ = 1; // keep sentinel at index 0
            freeHead_ = 0;
        }

        Node& operator[](NodeIndex idx) { return data_[idx]; }
        const Node& operator[](NodeIndex idx) const { return data_[idx]; }

    private:
        Node* data_ = nullptr;
        std::size_t size_ = 0;
        NodeIndex freeHead_ = 0;
#if defined(_WIN32)
        std::size_t committed_ = 0; // bytes committed so far
#endif

        void initReserve() {
            const std::size_t totalBytes = kMaxNodes * sizeof(Node);
#if defined(_WIN32)
            // Reserve virtual address space without committing physical memory
            data_ = static_cast<Node*>(::VirtualAlloc(nullptr, totalBytes, MEM_RESERVE, PAGE_NOACCESS));
            if (!data_) throw std::bad_alloc();
            committed_ = 0;
            // Commit first page for the null sentinel
            ensureCommitted(1);
#elif defined(__unix__) || defined(__APPLE__)
            // mmap reserves address space; pages are demand-faulted on first write
            data_ = static_cast<Node*>(
                ::mmap(nullptr, totalBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
            if (data_ == MAP_FAILED) {
                data_ = nullptr;
                throw std::bad_alloc();
            }
#else
            // Fallback: malloc
            data_ = static_cast<Node*>(std::malloc(totalBytes));
            if (!data_) throw std::bad_alloc();
#endif
            new (&data_[0]) Node(); // index 0 = null sentinel
            size_ = 1;
        }

        void ensureCommitted(std::size_t nodeCount) {
#if defined(_WIN32)
            std::size_t neededBytes = nodeCount * sizeof(Node);
            if (neededBytes > committed_) {
                // Round up to page boundary
                std::size_t newCommitted = (neededBytes + kPageSize - 1) & ~(kPageSize - 1);
                if (!::VirtualAlloc(data_, newCommitted, MEM_COMMIT, PAGE_READWRITE)) throw std::bad_alloc();
                committed_ = newCommitted;
            }
#else
            (void)nodeCount;
#endif
        }

        void releaseAll() {
            if (!data_) return;
            const std::size_t totalBytes = kMaxNodes * sizeof(Node);
#if defined(_WIN32)
            (void)totalBytes;
            ::VirtualFree(data_, 0, MEM_RELEASE);
#elif defined(__unix__) || defined(__APPLE__)
            ::munmap(data_, totalBytes);
#else
            std::free(data_);
#endif
            data_ = nullptr;
        }
    };

    NodeArena arena;
    std::vector<NodeIndex> holeQueue;
};

template <typename N>
template <typename Polygon>
void Earcut<N>::operator()(const Polygon& points) {
    // reset
    indices.clear();
    vertices = 0;

    if (points.empty()) return;

    double x;
    double y;
    int threshold = 80;
    std::size_t len = 0;

    for (size_t i = 0; threshold >= 0 && i < points.size(); i++) {
        threshold -= static_cast<int>(points[i].size());
        len += points[i].size();
    }

    indices.reserve(len + points[0].size());

    NodeIndex outerNode = linkedList(points[0], true);
    if (outerNode == 0 || arena[outerNode].prev == arena[outerNode].next) return;

    if (points.size() > 1) outerNode = eliminateHoles(points, outerNode);

    // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
    hashing = threshold < 0;
    if (hashing) {
        NodeIndex p = arena[outerNode].next;
        minX = maxX = arena[outerNode].x;
        minY = maxY = arena[outerNode].y;
        do {
            x = arena[p].x;
            y = arena[p].y;
            minX = std::min<double>(minX, x);
            minY = std::min<double>(minY, y);
            maxX = std::max<double>(maxX, x);
            maxY = std::max<double>(maxY, y);
            p = arena[p].next;
        } while (p != outerNode);

        // minX, minY and inv_size are later used to transform coords into integers for z-order calculation
        inv_size = std::max<double>(maxX - minX, maxY - minY);
        inv_size = inv_size != .0 ? (32767. / inv_size) : .0;
    }

    earcutLinked(outerNode);

    arena.reset();
    holeQueue.clear();
}

// create a circular doubly linked list from polygon points in the specified winding order
template <typename N>
template <typename Ring>
typename Earcut<N>::NodeIndex Earcut<N>::linkedList(const Ring& points, const bool clockwise) {
    using Point = typename Ring::value_type;
    double sum = 0;
    const std::size_t len = points.size();
    std::size_t i, j;
    NodeIndex last = 0;

    // calculate original winding order of a polygon ring
    for (i = 0, j = len > 0 ? len - 1 : 0; i < len; j = i++) {
        const auto& p1 = points[i];
        const auto& p2 = points[j];
        const double p20 = util::nth<0, Point>::get(p2);
        const double p10 = util::nth<0, Point>::get(p1);
        const double p11 = util::nth<1, Point>::get(p1);
        const double p21 = util::nth<1, Point>::get(p2);
        sum += (p20 - p10) * (p11 + p21);
    }

    // link points into circular doubly-linked list in the specified winding order
    if (clockwise == (sum > 0)) {
        for (i = 0; i < len; i++) last = insertNode(vertices + i, points[i], last);
    } else {
        for (i = len; i-- > 0;) last = insertNode(vertices + i, points[i], last);
    }

    if (last != 0 && equals(last, arena[last].next)) {
        removeNode(last);
        last = arena[last].next;
    }

    vertices += len;

    return last;
}

// eliminate colinear or duplicate points
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::filterPoints(NodeIndex start, NodeIndex end) {
    if (end == 0) end = start;

    NodeIndex p = start;
    bool again;
    do {
        again = false;

        if (!arena[p].steiner && (equals(p, arena[p].next) || area(arena[p].prev, p, arena[p].next) == 0)) {
            removeNode(p);
            p = end = arena[p].prev;

            if (p == arena[p].next) break;
            again = true;

        } else {
            p = arena[p].next;
        }
    } while (again || p != end);

    return end;
}

// main ear slicing loop which triangulates a polygon (given as a linked list)
template <typename N>
void Earcut<N>::earcutLinked(NodeIndex ear, int pass) {
    if (ear == 0) return;

    // interlink polygon nodes in z-order
    if (!pass && hashing) indexCurve(ear);

    NodeIndex stop = ear;
    NodeIndex prev;
    NodeIndex next;

    // iterate through ears, slicing them one by one
    while (arena[ear].prev != arena[ear].next) {
        prev = arena[ear].prev;
        next = arena[ear].next;

        if (hashing ? isEarHashed(ear) : isEar(ear)) {
            // cut off the triangle
            indices.emplace_back(static_cast<N>(arena[prev].i));
            indices.emplace_back(static_cast<N>(arena[ear].i));
            indices.emplace_back(static_cast<N>(arena[next].i));

            removeNode(ear);

            // skipping the next vertice leads to less sliver triangles
            ear = arena[next].next;
            stop = arena[next].next;

            continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if (ear == stop) {
            // try filtering points and slicing again
            if (!pass) earcutLinked(filterPoints(ear), 1);

            // if this didn't work, try curing all small self-intersections locally
            else if (pass == 1) {
                ear = cureLocalIntersections(filterPoints(ear));
                earcutLinked(ear, 2);

                // as a last resort, try splitting the remaining polygon into two
            } else if (pass == 2)
                splitEarcut(ear);

            break;
        }
    }
}

// check whether a polygon node forms a valid ear with adjacent nodes
template <typename N>
bool Earcut<N>::isEar(NodeIndex ear) {
    const NodeIndex a = arena[ear].prev;
    const NodeIndex c = arena[ear].next;

    // Create triangle with cached coordinates and bounding box
    const Triangle tri(arena[a], arena[ear], arena[c]);
    if (tri.area() >= 0) return false; // reflex, can't be an ear

    // now make sure we don't have other points inside the potential ear
    NodeIndex p = arena[c].next;

    while (p != a) {
        if (tri.containsPoint(arena[p].x, arena[p].y) && area(arena[p].prev, p, arena[p].next) >= 0) return false;
        p = arena[p].next;
    }

    return true;
}

template <typename N>
bool Earcut<N>::isEarHashed(NodeIndex ear) {
    const NodeIndex a = arena[ear].prev;
    const NodeIndex c = arena[ear].next;

    // Create triangle with cached coordinates and bounding box
    const Triangle tri(arena[a], arena[ear], arena[c]);
    if (tri.area() >= 0) return false; // reflex, can't be an ear

    // triangle bbox; min & max are calculated like this for speed
    const double minTX = std::min<double>(tri.ax, std::min<double>(tri.bx, tri.cx));
    const double minTY = std::min<double>(tri.ay, std::min<double>(tri.by, tri.cy));
    const double maxTX = std::max<double>(tri.ax, std::max<double>(tri.bx, tri.cx));
    const double maxTY = std::max<double>(tri.ay, std::max<double>(tri.by, tri.cy));

    // z-order range for the current triangle bbox;
    const int32_t minZ = zOrder(minTX, minTY);
    const int32_t maxZ = zOrder(maxTX, maxTY);

    // first look for points inside the triangle in increasing z-order
    NodeIndex p = arena[ear].nextZ;

    while (p != 0 && arena[p].z <= maxZ) {
        if (p != a && p != c && tri.containsPoint(arena[p].x, arena[p].y) &&
            area(arena[p].prev, p, arena[p].next) >= 0)
            return false;
        p = arena[p].nextZ;
    }

    // then look for points in decreasing z-order
    p = arena[ear].prevZ;

    while (p != 0 && arena[p].z >= minZ) {
        if (p != a && p != c && tri.containsPoint(arena[p].x, arena[p].y) &&
            area(arena[p].prev, p, arena[p].next) >= 0)
            return false;
        p = arena[p].prevZ;
    }

    return true;
}

// go through all polygon nodes and cure small local self-intersections
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::cureLocalIntersections(NodeIndex start) {
    NodeIndex p = start;
    do {
        NodeIndex a = arena[p].prev;
        NodeIndex pn = arena[p].next;
        NodeIndex b = arena[pn].next;

        // a self-intersection where edge (v[i-1],v[i]) intersects (v[i+1],v[i+2])
        if (!equals(a, b) && intersects(a, p, pn, b) && locallyInside(a, b) && locallyInside(b, a)) {
            indices.emplace_back(static_cast<N>(arena[a].i));
            indices.emplace_back(static_cast<N>(arena[p].i));
            indices.emplace_back(static_cast<N>(arena[b].i));

            // remove two nodes involved
            removeNode(p);
            removeNode(pn);

            p = start = b;
        }
        p = arena[p].next;
    } while (p != start);

    return filterPoints(p);
}

// try splitting polygon into two and triangulate them independently
template <typename N>
void Earcut<N>::splitEarcut(NodeIndex start) {
    // look for a valid diagonal that divides the polygon into two
    NodeIndex a = start;
    do {
        NodeIndex b = arena[arena[a].next].next;
        while (b != arena[a].prev) {
            if (arena[a].i != arena[b].i && isValidDiagonal(a, b)) {
                // split the polygon in two by the diagonal
                NodeIndex c = splitPolygon(a, b);

                // filter colinear points around the cuts
                a = filterPoints(a, arena[a].next);
                c = filterPoints(c, arena[c].next);

                // run earcut on each half
                earcutLinked(a);
                earcutLinked(c);
                return;
            }
            b = arena[b].next;
        }
        a = arena[a].next;
    } while (a != start);
}

// link every hole into the outer loop, producing a single-ring polygon without holes
template <typename N>
template <typename Polygon>
typename Earcut<N>::NodeIndex Earcut<N>::eliminateHoles(const Polygon& points, NodeIndex outerNode) {
    const size_t len = points.size();

    holeQueue.clear();
    for (size_t i = 1; i < len; i++) {
        NodeIndex list = linkedList(points[i], false);
        if (list != 0) {
            if (list == arena[list].next) arena[list].steiner = true;
            holeQueue.push_back(getLeftmost(list));
        }
    }
    std::sort(holeQueue.begin(), holeQueue.end(),
              [this](NodeIndex a, NodeIndex b) { return arena[a].x < arena[b].x; });

    // process holes from left to right
    for (size_t i = 0; i < holeQueue.size(); i++) {
        outerNode = eliminateHole(holeQueue[i], outerNode);
    }

    return outerNode;
}

// find a bridge between vertices that connects hole with an outer ring and and link it
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::eliminateHole(NodeIndex hole, NodeIndex outerNode) {
    NodeIndex bridge = findHoleBridge(hole, outerNode);
    if (bridge == 0) {
        return outerNode;
    }

    NodeIndex bridgeReverse = splitPolygon(bridge, hole);

    // filter collinear points around the cuts
    filterPoints(bridgeReverse, arena[bridgeReverse].next);

    // Check if input node was removed by the filtering
    return filterPoints(bridge, arena[bridge].next);
}

// David Eberly's algorithm for finding a bridge between hole and outer polygon
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::findHoleBridge(NodeIndex hole, NodeIndex outerNode) {
    NodeIndex p = outerNode;
    double hx = arena[hole].x;
    double hy = arena[hole].y;
    double qx = -std::numeric_limits<double>::infinity();
    NodeIndex m = 0;

    // find a segment intersected by a ray from the hole's leftmost Vertex to the left;
    // segment's endpoint with lesser x will be potential connection Vertex
    do {
        NodeIndex pNext = arena[p].next;
        if (hy <= arena[p].y && hy >= arena[pNext].y && arena[pNext].y != arena[p].y) {
            double x = arena[p].x + (hy - arena[p].y) * (arena[pNext].x - arena[p].x) / (arena[pNext].y - arena[p].y);
            if (x <= hx && x > qx) {
                qx = x;
                m = arena[p].x < arena[pNext].x ? p : pNext;
                if (x == hx) return m; // hole touches outer segment; pick leftmost endpoint
            }
        }
        p = pNext;
    } while (p != outerNode);

    if (m == 0) return 0;

    // look for points inside the triangle of hole Vertex, segment intersection and endpoint;
    // if there are no points found, we have a valid connection;
    // otherwise choose the Vertex of the minimum angle with the ray as connection Vertex

    const NodeIndex stop = m;
    double tanMin = std::numeric_limits<double>::infinity();
    double tanCur = 0;

    p = m;
    double mx = arena[m].x;
    double my = arena[m].y;

    do {
        if (hx >= arena[p].x && arena[p].x >= mx && hx != arena[p].x &&
            pointInTriangle(
                hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, arena[p].x, arena[p].y)) {
            tanCur = std::abs(hy - arena[p].y) / (hx - arena[p].x); // tangential

            if (locallyInside(p, hole) &&
                (tanCur < tanMin ||
                 (tanCur == tanMin && (arena[p].x > arena[m].x || sectorContainsSector(m, p))))) {
                m = p;
                tanMin = tanCur;
            }
        }

        p = arena[p].next;
    } while (p != stop);

    return m;
}

// whether sector in vertex m contains sector in vertex p in the same coordinates
template <typename N>
bool Earcut<N>::sectorContainsSector(NodeIndex m, NodeIndex p) {
    return area(arena[m].prev, m, arena[p].prev) < 0 && area(arena[p].next, m, arena[m].next) < 0;
}

// interlink polygon nodes in z-order
template <typename N>
void Earcut<N>::indexCurve(NodeIndex start) {
    assert(start != 0);
    NodeIndex p = start;

    do {
        if (arena[p].z == 0) arena[p].z = zOrder(arena[p].x, arena[p].y);
        arena[p].prevZ = arena[p].prev;
        arena[p].nextZ = arena[p].next;
        p = arena[p].next;
    } while (p != start);

    arena[arena[p].prevZ].nextZ = 0;
    arena[p].prevZ = 0;

    sortLinked(p);
}

// Simon Tatham's linked list merge sort algorithm
// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::sortLinked(NodeIndex list) {
    assert(list != 0);
    NodeIndex p;
    NodeIndex q;
    NodeIndex e;
    NodeIndex tail;
    int i, numMerges, pSize, qSize;
    int inSize = 1;

    for (;;) {
        p = list;
        list = 0;
        tail = 0;
        numMerges = 0;

        while (p != 0) {
            numMerges++;
            q = p;
            pSize = 0;
            for (i = 0; i < inSize; i++) {
                pSize++;
                q = arena[q].nextZ;
                if (q == 0) break;
            }

            qSize = inSize;

            while (pSize > 0 || (qSize > 0 && q != 0)) {
                if (pSize == 0) {
                    e = q;
                    q = arena[q].nextZ;
                    qSize--;
                } else if (qSize == 0 || q == 0) {
                    e = p;
                    p = arena[p].nextZ;
                    pSize--;
                } else if (arena[p].z <= arena[q].z) {
                    e = p;
                    p = arena[p].nextZ;
                    pSize--;
                } else {
                    e = q;
                    q = arena[q].nextZ;
                    qSize--;
                }

                if (tail != 0)
                    arena[tail].nextZ = e;
                else
                    list = e;

                arena[e].prevZ = tail;
                tail = e;
            }

            p = q;
        }

        arena[tail].nextZ = 0;

        if (numMerges <= 1) return list;

        inSize *= 2;
    }
}

// z-order of a Vertex given coords and size of the data bounding box
template <typename N>
int32_t Earcut<N>::zOrder(const double x_, const double y_) {
    // coords are transformed into non-negative 15-bit integer range
    int32_t x = static_cast<int32_t>((x_ - minX) * inv_size);
    int32_t y = static_cast<int32_t>((y_ - minY) * inv_size);

    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    return x | (y << 1);
}

// find the leftmost node of a polygon ring
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::getLeftmost(NodeIndex start) {
    NodeIndex p = start;
    NodeIndex leftmost = start;
    do {
        if (arena[p].x < arena[leftmost].x || (arena[p].x == arena[leftmost].x && arena[p].y < arena[leftmost].y))
            leftmost = p;
        p = arena[p].next;
    } while (p != start);

    return leftmost;
}

// check if a point lies within a convex triangle
template <typename N>
bool Earcut<N>::pointInTriangle(
    double ax, double ay, double bx, double by, double cx, double cy, double px, double py) const {
    return (cx - px) * (ay - py) >= (ax - px) * (cy - py) && (ax - px) * (by - py) >= (bx - px) * (ay - py) &&
           (bx - px) * (cy - py) >= (cx - px) * (by - py);
}

// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
template <typename N>
bool Earcut<N>::isValidDiagonal(NodeIndex a, NodeIndex b) {
    return arena[arena[a].next].i != arena[b].i && arena[arena[a].prev].i != arena[b].i &&
           !intersectsPolygon(a, b) && // dones't intersect other edges
           ((locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b) && // locally visible
             (area(arena[a].prev, a, arena[b].prev) != 0.0 ||
              area(a, arena[b].prev, b) != 0.0)) || // does not create opposite-facing sectors
            (equals(a, b) && area(arena[a].prev, a, arena[a].next) > 0 &&
             area(arena[b].prev, b, arena[b].next) > 0)); // special zero-length case
}

// signed area of a triangle
template <typename N>
double Earcut<N>::area(NodeIndex p, NodeIndex q, NodeIndex r) const {
    return (arena[q].y - arena[p].y) * (arena[r].x - arena[q].x) -
           (arena[q].x - arena[p].x) * (arena[r].y - arena[q].y);
}

// check if two points are equal
template <typename N>
bool Earcut<N>::equals(NodeIndex p1, NodeIndex p2) const {
    return arena[p1].x == arena[p2].x && arena[p1].y == arena[p2].y;
}

// check if two segments intersect
template <typename N>
bool Earcut<N>::intersects(NodeIndex p1, NodeIndex q1, NodeIndex p2, NodeIndex q2) const {
    int o1 = sign(area(p1, q1, p2));
    int o2 = sign(area(p1, q1, q2));
    int o3 = sign(area(p2, q2, p1));
    int o4 = sign(area(p2, q2, q1));

    if (o1 != o2 && o3 != o4) return true; // general case

    if (o1 == 0 && onSegment(p1, p2, q1)) return true; // p1, q1 and p2 are collinear and p2 lies on p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true; // p1, q1 and q2 are collinear and q2 lies on p1q1
    if (o3 == 0 && onSegment(p2, p1, q2)) return true; // p2, q2 and p1 are collinear and p1 lies on p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true; // p2, q2 and q1 are collinear and q1 lies on p2q2

    return false;
}

// for collinear points p, q, r, check if point q lies on segment pr
template <typename N>
bool Earcut<N>::onSegment(NodeIndex p, NodeIndex q, NodeIndex r) const {
    return arena[q].x <= std::max<double>(arena[p].x, arena[r].x) &&
           arena[q].x >= std::min<double>(arena[p].x, arena[r].x) &&
           arena[q].y <= std::max<double>(arena[p].y, arena[r].y) &&
           arena[q].y >= std::min<double>(arena[p].y, arena[r].y);
}

template <typename N>
int Earcut<N>::sign(double val) const {
    return (0.0 < val) - (val < 0.0);
}

// check if a polygon diagonal intersects any polygon segments
template <typename N>
bool Earcut<N>::intersectsPolygon(NodeIndex a, NodeIndex b) const {
    NodeIndex p = a;
    do {
        if (arena[p].i != arena[a].i && arena[arena[p].next].i != arena[a].i && arena[p].i != arena[b].i &&
            arena[arena[p].next].i != arena[b].i && intersects(p, arena[p].next, a, b))
            return true;
        p = arena[p].next;
    } while (p != a);

    return false;
}

// check if a polygon diagonal is locally inside the polygon
template <typename N>
bool Earcut<N>::locallyInside(NodeIndex a, NodeIndex b) const {
    return area(arena[a].prev, a, arena[a].next) < 0
               ? area(a, b, arena[a].next) >= 0 && area(a, arena[a].prev, b) >= 0
               : area(a, b, arena[a].prev) < 0 || area(a, arena[a].next, b) < 0;
}

// check if the middle Vertex of a polygon diagonal is inside the polygon
template <typename N>
bool Earcut<N>::middleInside(NodeIndex a, NodeIndex b) const {
    NodeIndex p = a;
    bool inside = false;
    double px = (arena[a].x + arena[b].x) / 2;
    double py = (arena[a].y + arena[b].y) / 2;
    do {
        NodeIndex pNext = arena[p].next;
        if (((arena[p].y > py) != (arena[pNext].y > py)) && arena[pNext].y != arena[p].y &&
            (px < (arena[pNext].x - arena[p].x) * (py - arena[p].y) / (arena[pNext].y - arena[p].y) + arena[p].x))
            inside = !inside;
        p = pNext;
    } while (p != a);

    return inside;
}

// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits
// polygon into two; if one belongs to the outer ring and another to a hole, it merges it into a
// single ring
template <typename N>
typename Earcut<N>::NodeIndex Earcut<N>::splitPolygon(NodeIndex a, NodeIndex b) {
    NodeIndex a2 = arena.construct(static_cast<N>(arena[a].i), arena[a].x, arena[a].y);
    NodeIndex b2 = arena.construct(static_cast<N>(arena[b].i), arena[b].x, arena[b].y);
    NodeIndex an = arena[a].next;
    NodeIndex bp = arena[b].prev;

    arena[a].next = b;
    arena[b].prev = a;

    arena[a2].next = an;
    arena[an].prev = a2;

    arena[b2].next = a2;
    arena[a2].prev = b2;

    arena[bp].next = b2;
    arena[b2].prev = bp;

    return b2;
}

// create a node and optionally link it with previous one (in a circular doubly linked list)
template <typename N>
template <typename Point>
typename Earcut<N>::NodeIndex Earcut<N>::insertNode(std::size_t i, const Point& pt, NodeIndex last) {
    NodeIndex p = arena.construct(static_cast<N>(i), util::nth<0, Point>::get(pt), util::nth<1, Point>::get(pt));

    if (last == 0) {
        arena[p].prev = p;
        arena[p].next = p;

    } else {
        arena[p].next = arena[last].next;
        arena[p].prev = last;
        arena[arena[last].next].prev = p;
        arena[last].next = p;
    }
    return p;
}

template <typename N>
void Earcut<N>::removeNode(NodeIndex p) {
    arena[arena[p].next].prev = arena[p].prev;
    arena[arena[p].prev].next = arena[p].next;

    if (arena[p].prevZ != 0) arena[arena[p].prevZ].nextZ = arena[p].nextZ;
    if (arena[p].nextZ != 0) arena[arena[p].nextZ].prevZ = arena[p].prevZ;
}
} // namespace detail

template <typename N = uint32_t, typename Polygon>
std::vector<N> earcut(const Polygon& poly) {
    mapbox::detail::Earcut<N> earcut;
    earcut(poly);
    return std::move(earcut.indices);
}
} // namespace mapbox
