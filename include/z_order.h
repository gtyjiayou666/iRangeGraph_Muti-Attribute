#pragma once
#include <iostream>
#include <vector>
#include <utility>
#include <cstdint>
#include <algorithm>

using Range = std::pair<uint64_t, uint64_t>;
using ZRangeList = std::vector<Range>;



uint64_t z_order_encode(uint32_t x, uint32_t y) {
    auto spread = [](uint32_t a) -> uint64_t {
        uint64_t out = 0;
        for (int i = 0; i < 16; ++i) {
            if (a & (1U << i)) {
                out |= (1ULL << (2 * i));
            }
        }
        return out;
    };
    return spread(x) | (spread(y) << 1);
}

struct ZInterval {
    uint64_t start;
    uint64_t end;
};

void decompose_optimized(
    uint32_t x1, uint32_t x2, 
    uint32_t y1, uint32_t y2,
    uint32_t bx, uint32_t by, 
    uint32_t width,
    uint32_t min_block_size,
    uint32_t max_z_value,
    std::vector<ZInterval>& result
) {
    if (x1 > x2 || y1 > y2) {
        return;
    }
    uint64_t current_min_z = z_order_encode(bx, by);
    
    if (current_min_z > max_z_value) {
        return;
    }
    uint32_t end_x = bx + width - 1;
    uint32_t end_y = by + width - 1;
    
    if (end_x < x1 || end_y < y1 || bx > x2 || by > y2) return;

    if (bx >= x1 && end_x <= x2 && by >= y1 && end_y <= y2) {
        result.push_back({z_order_encode(bx, by), z_order_encode(end_x, end_y)});
        return;
    }

    if (width <= min_block_size) {
        result.push_back({z_order_encode(bx, by), z_order_encode(end_x, end_y)});
        return;
    }

    uint32_t half = width / 2;
    decompose_optimized(x1, x2, y1, y2, bx, by, half, min_block_size, max_z_value, result);
    decompose_optimized(x1, x2, y1, y2, bx + half, by, half, min_block_size, max_z_value, result);
    decompose_optimized(x1, x2, y1, y2, bx, by + half, half, min_block_size, max_z_value, result);
    decompose_optimized(x1, x2, y1, y2, bx + half, by + half, half, min_block_size, max_z_value, result);
}


struct IndexRange {
    size_t start_idx;
    size_t end_idx;
};

std::vector<IndexRange> get_index_ranges_merged_optimized(
    const std::vector<int>& sorted_data, 
    uint32_t x1, uint32_t x2, 
    uint32_t y1, uint32_t y2
) {
    std::vector<IndexRange> index_ranges;
    
    if (sorted_data.empty()) {
        return index_ranges;
    }

    std::vector<ZInterval> z_intervals;
    // uint32_t max_width = 1 << 14; 
    uint32_t max_width = 1 << 16; 

    decompose_optimized(x1, x2, y1, y2, 0, 0, max_width, 64, sorted_data.back(), z_intervals);

    if (z_intervals.empty()) {
        return index_ranges;
    }
    index_ranges.reserve(sorted_data.size());

    auto search_start_it = sorted_data.begin();

    for (const auto& z_interval : z_intervals) {
        auto it_start = std::lower_bound(
            search_start_it, sorted_data.end(), z_interval.start,
            [](const int& p, uint64_t val) { return p < val; }
        );

        if (it_start == sorted_data.end()) {
            break;
        }

        auto it_end = std::upper_bound(
            it_start, sorted_data.end(), z_interval.end,
            [](uint64_t val, const int& p) { return val < p; }
        );

        if (it_start < it_end) {
            size_t start_idx = std::distance(sorted_data.begin(), it_start);
            size_t end_idx = std::distance(sorted_data.begin(), it_end) - 1;
            
            if (!index_ranges.empty()) {
                IndexRange& last = index_ranges.back();
                
                if (start_idx == last.end_idx + 1) {
                    last.end_idx = end_idx;
                    search_start_it = it_end; 
                    continue; 
                }
            }
            index_ranges.push_back({start_idx, end_idx});
        }

        search_start_it = it_end;
    }

    return index_ranges;
}

std::vector<ZInterval> query_z_order_or(uint32_t x1, uint32_t x2, uint32_t y1, uint32_t y2, uint32_t x1_global, uint32_t x2_global,  uint32_t y1_global, uint32_t y2_global, uint32_t max_z_value) {
    std::vector<ZInterval> result;
    
    uint32_t max_width = 1 << 14; 
    
    decompose_optimized(x1, x2, y1_global, y2_global, 0, 0, max_width, 64, max_z_value, result);
    decompose_optimized(x1_global, x2_global, y1, y2, 0, 0, max_width, 64, max_z_value, result);
    
    return result;
}

void merge_z_intervals_inplace(std::vector<ZInterval>& intervals) {
    if (intervals.size() <= 1) return;

    std::sort(intervals.begin(), intervals.end(),
              [](const ZInterval& a, const ZInterval& b) { return a.start < b.start; });

    size_t write_idx = 0;
    for (size_t read_idx = 1; read_idx < intervals.size(); ++read_idx) {
        if (intervals[read_idx].start <= intervals[write_idx].end + 1) {
            intervals[write_idx].end = std::max(intervals[write_idx].end, intervals[read_idx].end);
        } else {
            intervals[++write_idx] = intervals[read_idx];
        }
    }
    intervals.resize(write_idx + 1); 
}


std::vector<IndexRange> get_index_ranges_merged_or(
    const std::vector<int>& sorted_data, 
    uint32_t x1, uint32_t x2, 
    uint32_t y1, uint32_t y2,
    uint32_t x1_global, uint32_t x2_global, 
    uint32_t y1_global, uint32_t y2_global
) {
    int xxx = sorted_data.size();
    std::vector<ZInterval> raw_intervals = query_z_order_or(x1, x2, y1, y2, x1_global, x2_global, y1_global, y2_global, sorted_data[xxx - 1]);


    merge_z_intervals_inplace(raw_intervals);

    auto search_start_it = sorted_data.begin();

    std::vector<IndexRange> index_ranges;
    index_ranges.reserve(raw_intervals.size());
    for (const auto& z_interval : raw_intervals) {
        auto it_start = std::lower_bound(
            search_start_it, sorted_data.end(), z_interval.start,
            [](const int& p, uint64_t val) { return p < val; }
        );

        if (it_start == sorted_data.end()) {
            break;
        }

        auto it_end = std::upper_bound(
            it_start, sorted_data.end(), z_interval.end,
            [](uint64_t val, const int& p) { return val < p; }
        );

        if (it_start < it_end) {
            size_t start_idx = std::distance(sorted_data.begin(), it_start);
            size_t end_idx = std::distance(sorted_data.begin(), it_end) - 1;
            
            if (!index_ranges.empty()) {
                IndexRange& last = index_ranges.back();
                
                if (start_idx == last.end_idx + 1) {
                    last.end_idx = end_idx; 
                    search_start_it = it_end; 
                    continue; 
                }
            }

            index_ranges.push_back({start_idx, end_idx});
        }

        search_start_it = it_end;
    }
    return index_ranges;
}

