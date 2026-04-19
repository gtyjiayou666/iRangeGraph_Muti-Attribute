#pragma once

#include <vector>
#include "utils.h"
#include "searcher.hpp"
#include "memory.hpp"
#include <bitset>
#include <stack>
#include <cmath>
#include <climits>
#include "z_order.h"

namespace iRangeGraph
{
    template <typename dist_t>
    class iRangeGraph_Search_Zorder
    {
    public:
        iRangeGraph_multi::DataLoader *storage;
        SegmentTree *tree;
        size_t max_elements_{0};
        size_t dim_{0};
        size_t M_out{0};
        size_t ef_construction{0};

        size_t size_data_per_element_{0};
        size_t size_links_per_element_{0};
        size_t data_size_{0};

        size_t size_links_per_layer_{0};
        size_t offsetData_{0};
        int range_num_avg{0};

        char *data_memory_{nullptr};

        hnswlib::L2Space *space;
        hnswlib::DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_{nullptr};

        size_t metric_distance_computations{0};
        size_t metric_hops{0};

        int prefetch_lines{0};

        iRangeGraph_Search_Zorder(std::string edgefilename, iRangeGraph_multi::DataLoader *store, int M) : storage(store)
        {
            std::ifstream edgefile(edgefilename, std::ios::in | std::ios::binary);
            if (!edgefile.is_open())
                throw Exception("cannot open " + edgefilename);

            max_elements_ = storage->data_nb;
            dim_ = storage->Dim;
            tree = new iRangeGraph::SegmentTree(max_elements_);
            tree->BuildTree(tree->root);

            space = new hnswlib::L2Space(dim_);
            fstdistfunc_ = space->get_dist_func();
            dist_func_param_ = space->get_dist_func_param();
            M_out = M;

            data_size_ = dim_ * sizeof(float);
            size_links_per_layer_ = M_out * sizeof(tableint) + sizeof(linklistsizeint);
            size_links_per_element_ = size_links_per_layer_ * (tree->max_depth + 1);
            size_data_per_element_ = size_links_per_element_ + data_size_;
            offsetData_ = size_links_per_element_;

            data_memory_ = (char *)malloc(max_elements_ * size_data_per_element_);
            if (data_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            for (int pid = 0; pid < max_elements_; pid++)
            {
                for (int layer = 0; layer <= tree->max_depth; layer++)
                {
                    linklistsizeint *data = get_linklist(pid, layer);
                    edgefile.read((char *)data, sizeof(tableint));
                    int size = getListCount(data);
                    if (size > M_out)
                        throw Exception("real linklist size is bigger than defined M_out");
                    for (int i = 0; i < size; i++)
                    {
                        char *current_neighbor_ = (char *)(data + 1 + i);
                        edgefile.read(current_neighbor_, sizeof(tableint));
                    }
                }

                size_t size_in_bytes = dim_ * sizeof(float);
                char *data = getDataByInternalId(pid);
                std::memcpy(data, reinterpret_cast<char *>(storage->data_points[pid].data()), size_in_bytes);
            }
            edgefile.close();
        }

        ~iRangeGraph_Search_Zorder()
        {
            free(data_memory_);
            data_memory_ = nullptr;
        }

        inline char *getDataByInternalId(tableint internal_id) const
        {
            return (data_memory_ + internal_id * size_data_per_element_ + offsetData_);
        }

        linklistsizeint *get_linklist(tableint internal_id, int layer) const
        {
            return (linklistsizeint *)(data_memory_ + internal_id * size_data_per_element_ + layer * size_links_per_layer_);
        }

        int getListCount(linklistsizeint *ptr) const
        {
            return *((int *)ptr);
        }

        int GetOverLap(int l, int r, int ql, int qr)
        {
            int L = std::max(l, ql);
            int R = std::min(r, qr);
            return R - L + 1;
        }

        struct GapInfo
        {
            int index;
            int len;
        };
        struct CompareGap
        {
            bool operator()(const GapInfo &a, const GapInfo &b)
            {
                return a.len < b.len;
            }
        };
        void merge_small_gaps(const std::vector<IndexRange> &Q,
                              std::vector<IndexRange> &Q1,
                              std::vector<int> &mapping, double p)
        {
            int n = Q.size();
            mapping.assign(n, 0);

            double total_data = Q[n - 1].end_idx - Q[0].start_idx + 1;
            double total_gaps = 0;

            std::priority_queue<GapInfo, std::vector<GapInfo>, CompareGap> max_heap_gaps;

            for (int i = 0; i < n - 1; ++i)
            {
                int g_l = Q[i].end_idx + 1;
                int g_r = Q[i + 1].start_idx - 1;
                max_heap_gaps.push({i, g_r - g_l + 1});
                total_gaps = total_gaps + g_r - g_l + 1;
            }
            int total_nogaps = total_data - total_gaps;

            std::vector<bool> is_keeper(n - 1, false);

            int s = total_data - (p * total_nogaps);
            if (s <= 0)
            {
                Q1.push_back({Q[0].start_idx, Q[n - 1].end_idx});
                return;
            }
            int current_gaps = 0;
            while (!max_heap_gaps.empty())
            {
                GapInfo gap = max_heap_gaps.top();
                current_gaps = current_gaps + gap.len;
                if (current_gaps >= s)
                {
                    is_keeper[gap.index] = true;
                    break;
                }
                max_heap_gaps.pop();
                is_keeper[gap.index] = true;
            }
            int current_l = Q[0].start_idx;
            int current_r = Q[0].end_idx;

            int new_idx = 0;
            mapping[0] = new_idx;
            for (int i = 0; i < n - 1; ++i)
            {
                if (is_keeper[i])
                {
                    Q1.push_back({current_l, current_r});
                    new_idx++;
                    current_l = Q[i + 1].start_idx;
                    current_r = Q[i + 1].end_idx;
                    mapping[i + 1] = new_idx;
                }
                else
                {
                    current_r = Q[i + 1].end_idx;
                    mapping[i + 1] = new_idx;
                }
            }
            Q1.push_back({current_l, current_r});
        }
        inline bool CheckInQueryRange(int pid, std::pair<int, int> &queryrange)
        {
            int val = storage->attributes[pid][1];
            if (val < queryrange.first || val > queryrange.second)
                return false;
            return true;
        }
        inline bool CheckInQueryRangeOr(int pid, std::vector<std::pair<int, int>> &queryrange)
        {
            int val = storage->attributes[pid][0];
            if (val >= queryrange[0].first && val <= queryrange[0].second)
                return true;
            val = storage->attributes[pid][1];
            if (val >= queryrange[1].first && val <= queryrange[1].second)
                return true;
            return false;
        }
        std::vector<tableint> SelectEdge(int pid, int ql, int qr, int edge_limit, searcher::Bitset<uint64_t> &visited_set)
        {
            TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<tableint> selected_edges;
            selected_edges.reserve(edge_limit);
            do
            {
                cur_node = nxt_node;
                bool contain = false;
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, ql, qr) == GetOverLap(nxt_node->lbound, nxt_node->rbound, ql, qr))
                        {
                            cur_node = nxt_node;
                            contain = true;
                        }
                    }
                } while (contain);

                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);

                for (size_t j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    // if (visitedpool[neighborId] == visited_tag)
                    //     continue;
                    if (visited_set.get(neighborId))
                        continue;
                    selected_edges.emplace_back(neighborId);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }

            } while (cur_node->lbound < ql || cur_node->rbound > qr);

            return selected_edges;
        }

        std::vector<tableint> SelectEdgeAnd(int pid, int z_key, int ql, int qr, int edge_limit, searcher::Bitset<uint64_t> &visited_set)
        {
            TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<tableint> selected_edges;
            selected_edges.reserve(edge_limit);
            do
            {
                cur_node = nxt_node;
                bool contain = false;
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(storage->zorder[cur_node->lbound], storage->zorder[cur_node->rbound], ql, qr) == GetOverLap(storage->zorder[nxt_node->lbound], storage->zorder[nxt_node->rbound], ql, qr))
                        {
                            if (cur_node->lbound < nxt_node->lbound)
                            {
                                if (storage->zorder[nxt_node->lbound] != ql || storage->zorder[nxt_node->lbound - 1] != ql)
                                {
                                    cur_node = nxt_node;
                                    contain = true;
                                }
                            }
                            else
                            {
                                if (storage->zorder[nxt_node->rbound] != qr || storage->zorder[nxt_node->rbound + 1] != qr)
                                {
                                    cur_node = nxt_node;
                                    contain = true;
                                }
                            }
                        }
                    }
                } while (contain);

                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);

                for (size_t j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    int z_key1 = storage->zorder[neighborId];
                    if (z_key1 < ql || z_key1 > qr)
                        continue;
                    // if (visitedpool[neighborId] == visited_tag)
                    //     continue;
                    if (visited_set.get(neighborId))
                        continue;
                    selected_edges.emplace_back(neighborId);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }

            } while (storage->zorder[cur_node->lbound] < ql || storage->zorder[cur_node->rbound] > qr);

            return selected_edges;
        }

        std::priority_queue<PFI> TopDown_nodeentries_search_zorder_and(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<IndexRange> &Q, int edge_limit, int id, std::pair<int, int> queryrange)
        {
            std::vector<IndexRange> Q1;
            std::vector<int> mapping;
            merge_small_gaps(Q, Q1, mapping, 1.5);
            int range_num = Q1.size();
            range_num_avg = range_num_avg + range_num;

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFIMI, std::vector<PFIMI>, Compare> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                {
                    continue;
                }
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    if (CheckInQueryRange(pid, queryrange))
                    {
                        top_candidates.emplace(dis, pid);
                        lowerBound = top_candidates.top().first;
                    }
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    if (CheckInQueryRange(pid, queryrange))
                    {
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            bool erfen = range_num < 16;
            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.dis > lowerBound)
                    break;
                candidate_set.pop();
                int current_pid = current_point_pair.id;
                int range_id = current_point_pair.range_id;
                if (range_id == -1)
                {
                    if (erfen)
                    {
                        for (int i = 0; i < range_num; i++)
                        {
                            if (current_pid >= Q1[i].start_idx)
                            {
                                if (current_pid <= Q1[i].end_idx)
                                {
                                    range_id = i;
                                    break;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        auto it = std::lower_bound(Q1.begin(), Q1.end(), current_pid + 1,
                                                   [](const IndexRange &a, size_t key)
                                                   {
                                                       return a.start_idx < key;
                                                   });

                        if (it != Q1.begin())
                        {
                            IndexRange &candidate = *(it - 1);
                            if (current_pid <= candidate.end_idx)
                            {
                                range_id = static_cast<int>(std::distance(Q1.begin(), it - 1));
                            }
                            else
                            {
                                range_id = -1; 
                            }
                        }
                        else
                        {
                            range_id = -1;
                        }
                    }
                }
                auto selected_edges = SelectEdge(current_pid, Q1[range_id].start_idx, Q1[range_id].end_idx, edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (CheckInQueryRange(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (CheckInQueryRange(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }

                if (range_num == 1)
                {
                    continue;
                }
                std::vector<tableint> edges;
                edges.reserve(edge_limit / 2);
                JINS_and(Q, current_pid, edge_limit / 2, 2, edges, queryrange);
                int num_edges1 = edges.size();
                for (int i = 0; i < std::min(num_edges1, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges1; ++i)
                {
                    int neighbor_id = edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;
                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (CheckInQueryRange(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (CheckInQueryRange(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }
        std::priority_queue<PFI> TopDown_nodeentries_search_zorder_or(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, std::vector<IndexRange> &Q, int edge_limit, int id, std::vector<std::pair<int, int>> queryrange)
        {
            std::vector<IndexRange> Q1;
            std::vector<int> mapping;
            merge_small_gaps(Q, Q1, mapping, 1.5);
            int range_num = Q1.size();
            range_num_avg = range_num_avg + range_num;

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFIMI, std::vector<PFIMI>, Compare> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            float lowerBound = std::numeric_limits<float>::max();
            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                {
                    continue;
                }
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                if (top_candidates.size() < ef)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    if (CheckInQueryRangeOr(pid, queryrange))
                    {
                        top_candidates.emplace(dis, pid);
                        lowerBound = top_candidates.top().first;
                    }
                }
                else if (dis < lowerBound)
                {
                    candidate_set.emplace(dis, pid, mapping[id]);
                    if (CheckInQueryRangeOr(pid, queryrange))
                    {
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            bool erfen = range_num < 16;
            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.dis > lowerBound)
                    break;
                candidate_set.pop();
                int current_pid = current_point_pair.id;
                int range_id = current_point_pair.range_id;
                if (range_id == -1)
                {
                    if (erfen)
                    {
                        for (int i = 0; i < range_num; i++)
                        {
                            if (current_pid >= Q1[i].start_idx)
                            {
                                if (current_pid <= Q1[i].end_idx)
                                {
                                    range_id = i;
                                    break;
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        auto it = std::lower_bound(Q1.begin(), Q1.end(), current_pid + 1,
                                                   [](const IndexRange &a, size_t key)
                                                   {
                                                       return a.start_idx < key;
                                                   });

                        if (it != Q1.begin())
                        {
                            IndexRange &candidate = *(it - 1);
                            if (current_pid <= candidate.end_idx)
                            {
                                range_id = static_cast<int>(std::distance(Q1.begin(), it - 1));
                            }
                            else
                            {
                                range_id = -1;
                            }
                        }
                        else
                        {
                            range_id = -1;
                        }
                    }
                }
                auto selected_edges = SelectEdge(current_pid, Q1[range_id].start_idx, Q1[range_id].end_idx, edge_limit, visited_set);
                int num_edges = selected_edges.size();
                for (int i = 0; i < std::min(num_edges, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(selected_edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges; ++i)
                {
                    int neighbor_id = selected_edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;

                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (CheckInQueryRangeOr(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, range_id);
                        if (CheckInQueryRangeOr(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }

                if (range_num == 1)
                {
                    continue;
                }
                std::vector<tableint> edges;
                edges.reserve(edge_limit / 2);
                JINS_or(Q, current_pid, edge_limit / 2, 2, edges, queryrange);
                int num_edges1 = edges.size();
                for (int i = 0; i < std::min(num_edges1, 3); ++i)
                {
                    memory::mem_prefetch_L1(getDataByInternalId(edges[i]), this->prefetch_lines);
                }
                for (int i = 0; i < num_edges1; ++i)
                {
                    int neighbor_id = edges[i];
                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    ++metric_distance_computations;
                    if (top_candidates.size() < ef)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (CheckInQueryRangeOr(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            lowerBound = top_candidates.top().first;
                        }
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id, -1);
                        if (CheckInQueryRangeOr(neighbor_id, queryrange))
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            top_candidates.pop();
                            lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }
        void JINS_or(std::vector<IndexRange> &Q, int pid, int edge_limit, int deep, std::vector<tableint> &edges, std::vector<std::pair<int, int>> &queryrange)
        {
            int Qsize = Q.size();
            if (deep == 0 || Qsize <= 1)
            {
                return;
            }
            int start_idx = 0;
            int end_idx = Qsize - 1;
            int left = Q[0].start_idx;
            int right = Q[end_idx].end_idx;
            TreeNode *cur_node = nullptr;
            TreeNode *nxt_node = tree->root;
            cur_node = nxt_node;
            bool contain = false;
            while (cur_node)
            {
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, left, right) == GetOverLap(nxt_node->lbound, nxt_node->rbound, left, right))
                        {
                            if (cur_node->lbound < nxt_node->lbound)
                            {
                                cur_node = nxt_node;
                                contain = true;
                            }
                        }
                    }
                } while (contain);
                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);
                for (int j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < left || neighborId > right)
                        continue;
                    if (CheckInQueryRangeOr(neighborId, queryrange))
                    {
                        edges.emplace_back(neighborId);
                    }
                    else
                    {
                        JINS_or(Q, neighborId, edge_limit, deep - 1, edges, queryrange);
                    }
                    if (edges.size() >= edge_limit)
                    {
                        return;
                    }
                }
                if (nxt_node == nullptr)
                {
                    return;
                }
                cur_node = nxt_node;
                int low = start_idx;
                int high = end_idx;
                start_idx = end_idx + 1;

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1); 
                    if (Q[mid].end_idx >= cur_node->lbound)
                    {
                        start_idx = mid;
                        high = mid - 1;
                    }
                    else
                    {
                        low = mid + 1;
                    }
                }

                if (start_idx > end_idx)
                {
                    return;
                }
                low = start_idx;
                high = end_idx;
                end_idx = start_idx;

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1);
                    if (Q[mid].start_idx > cur_node->rbound)
                    {
                        high = mid - 1;
                    }
                    else
                    {
                        end_idx = mid;
                        low = mid + 1;
                    }
                }

                left = Q[start_idx].start_idx;
                right = Q[end_idx].end_idx;
            }
        }

        void JINS_and(std::vector<IndexRange> &Q, int pid, int edge_limit, int deep, std::vector<tableint> &edges, std::pair<int, int> &queryrange)
        {
            int Qsize = Q.size();
            if (deep == 0 || Qsize <= 1)
            {
                return;
            }
            int start_idx = 0;
            int end_idx = Qsize - 1;
            int left = Q[0].start_idx;
            int right = Q[end_idx].end_idx;
            TreeNode *cur_node = nullptr;
            TreeNode *nxt_node = tree->root;
            cur_node = nxt_node;
            bool contain = false;
            while (cur_node)
            {
                do
                {
                    contain = false;
                    if (cur_node->childs.size() == 0)
                        nxt_node = nullptr;
                    else
                    {
                        for (int i = 0; i < cur_node->childs.size(); ++i)
                        {
                            if (cur_node->childs[i]->lbound <= pid && cur_node->childs[i]->rbound >= pid)
                            {
                                nxt_node = cur_node->childs[i];
                                break;
                            }
                        }
                        if (GetOverLap(cur_node->lbound, cur_node->rbound, left, right) == GetOverLap(nxt_node->lbound, nxt_node->rbound, left, right))
                        {
                            if (cur_node->lbound < nxt_node->lbound)
                            {
                                cur_node = nxt_node;
                                contain = true;
                            }
                        }
                    }
                } while (contain);
                int *data = (int *)get_linklist(pid, cur_node->depth);
                size_t size = getListCount((linklistsizeint *)data);
                for (int j = 1; j <= size; ++j)
                {
                    int neighborId = *(data + j);
                    if (neighborId < left || neighborId > right)
                        continue;
                    if (CheckInQueryRange(neighborId, queryrange))
                    {
                        edges.emplace_back(neighborId);
                    }
                    else
                    {
                        JINS_and(Q, neighborId, edge_limit, deep - 1, edges, queryrange);
                    }
                    if (edges.size() >= edge_limit)
                    {
                        return;
                    }
                }
                if (nxt_node == nullptr)
                {
                    return;
                }
                cur_node = nxt_node;
                int low = start_idx;
                int high = end_idx;
                start_idx = end_idx + 1; 

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1); 
                    if (Q[mid].end_idx >= cur_node->lbound)
                    {
                        start_idx = mid;
                        high = mid - 1;
                    }
                    else
                    {
                        low = mid + 1;
                    }
                }

                if (start_idx > end_idx)
                {
                    return;
                }
                low = start_idx;
                high = end_idx;
                end_idx = start_idx;

                while (low <= high)
                {
                    int mid = low + ((high - low) >> 1);
                    if (Q[mid].start_idx > cur_node->rbound)
                    {
                        high = mid - 1;
                    }
                    else
                    {
                        end_idx = mid;
                        low = mid + 1;
                    }
                }

                left = Q[start_idx].start_idx;
                right = Q[end_idx].end_idx;
            }
        }
        void searchor(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int x_min_global = 1000000;
            int x_max_global = -1;
            int y_min_global = 1000000;
            int y_max_global = -1;
            for (int i = 0; i < storage->data_nb; i++)
            {
                int val1 = storage->attributes[i][0];
                if (val1 > x_max_global)
                    x_max_global = val1;
                if (val1 < x_min_global)
                    x_min_global = val1;
                int val2 = storage->attributes[i][1];
                if (val2 > y_max_global)
                    y_max_global = val2;
                if (val2 < y_min_global)
                    y_min_global = val2;
            }
            for (auto range : storage->query_range_or)
            {
                std::string domain = range.first;
                std::vector<std::vector<int>> &gt = storage->ground_truth[domain];
                std::string savepath = saveprefix + domain + "-zorder-or.csv";
                CheckPath(savepath);

                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }

                std::vector<int> HOP;
                std::vector<int> DCO;
                std::vector<float> QPS;
                std::vector<float> RECALL;

                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;
                    range_num_avg = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        auto cons = range.second[i].attr_constraints;
                        std::vector<int> ql, qr;
                        timeval t1, t2;
                        int x_min = cons[0].first;
                        int x_max = cons[0].second;
                        int y_min = cons[1].first;
                        int y_max = cons[1].second;
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::default_random_engine e(seed);
                        gettimeofday(&t1, NULL);
                        std::vector<IndexRange> intervals = get_index_ranges_merged_or(storage->zorder, x_min, x_max, y_min, y_max, x_min_global, x_max_global, y_min_global, y_max_global);
                        int range_num = intervals.size();
                        // gettimeofday(&t1, NULL); //cache

                        std::uniform_int_distribution<int> dist(0, range_num - 1);
                        int id = dist(e);
                        std::vector<TreeNode *> node = tree->range_filter(tree->root, intervals[id].start_idx, intervals[id].end_idx);
                        std::priority_queue<PFI> res = TopDown_nodeentries_search_zorder_or(node, storage->query_points[i].data(), ef, storage->query_K, intervals, edge_limit, id, cons);
                        gettimeofday(&t2, NULL);
                        searchtime += GetTime(t1, t2);

                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }
                    range_num_avg = range_num_avg / storage->query_nb;

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << "," << range_num_avg << std::endl;
                    if (recall < 0.9)
                    {
                        break;
                    }
                }

                outfile.close();
            }
        }
        void searchand(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            int x_min = 1000000;
            int x_max = -1;
            for (int i = 0; i < storage->data_nb; i++)
            {
                int val = storage->attributes[i][0];
                if (val > x_max)
                    x_max = val;
                if (val < x_min)
                    x_min = val;
            }

            for (auto range : storage->query_range)
            {
                std::string domain = range.first;
                std::vector<std::vector<int>> &gt = storage->ground_truth[domain];
                std::string savepath = saveprefix + domain + "-zorder_and.csv";
                CheckPath(savepath);

                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }

                std::vector<int> HOP;
                std::vector<int> DCO;
                std::vector<float> QPS;
                std::vector<float> RECALL;

                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;
                    range_num_avg = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {
                        auto cons = range.second[i];
                        std::vector<int> ql, qr;
                        timeval t1, t2;
                        int y_min = cons.first;
                        int y_max = cons.second;
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::default_random_engine e(seed);

                        gettimeofday(&t1, NULL);
                        std::vector<IndexRange> intervals = get_index_ranges_merged_optimized(storage->zorder, x_min, x_max, y_min, y_max);
                        // gettimeofday(&t1, NULL); //cache
                        int range_num = intervals.size();
                        std::uniform_int_distribution<int> dist(0, range_num - 1);
                        int id = dist(e);
                        std::vector<TreeNode *> node = tree->range_filter(tree->root, intervals[id].start_idx, intervals[id].end_idx);
                        std::priority_queue<PFI> res = TopDown_nodeentries_search_zorder_and(node, storage->query_points[i].data(), ef, storage->query_K, intervals, edge_limit, id, cons);
                        gettimeofday(&t2, NULL);
                        searchtime += GetTime(t1, t2);

                        std::map<int, int> record;
                        while (res.size())
                        {
                            auto x = res.top().second;
                            res.pop();
                            if (record.count(x))
                                throw Exception("repetitive search results");
                            record[x] = 1;
                            if (std::find(gt[i].begin(), gt[i].end(), x) != gt[i].end())
                                tp++;
                        }
                    }
                    range_num_avg = range_num_avg / storage->query_nb;

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << "," << range_num_avg << std::endl;
                    if (recall < 0.9)
                    {
                        break;
                    }
                }

                outfile.close();
            }
        }
    };
}