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
    typedef std::pair<float, std::pair<int, int>> PFII;
    template <typename dist_t>
    class iRangeGraph_Search_Lex
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

        char *data_memory_{nullptr};

        hnswlib::L2Space *space;
        hnswlib::DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_{nullptr};

        size_t metric_distance_computations{0};
        size_t metric_hops{0};

        int prefetch_lines{0};

        iRangeGraph_Search_Lex(std::string edgefilename, iRangeGraph_multi::DataLoader *store, int M) : storage(store)
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

        ~iRangeGraph_Search_Lex()
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

        inline bool CheckInQueryRange(int pid, std::pair<int, int> &queryrange)
        {
            int val = storage->attributes[pid][1];
            if (val < queryrange.first || val > queryrange.second)
                return false;
            return true;
        }

        std::vector<std::pair<tableint, bool>> SelectEdge(int pid, int ql, int qr, int edge_limit, std::pair<int, int> &queryrange)
        {
            iRangeGraph::TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<std::pair<tableint, bool>> selected_edges;
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
                        for (int i = 0; i < cur_node->childs.size(); i++)
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

                for (size_t j = 1; j <= size; j++)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    // int prob = 1;
                    bool inrange = CheckInQueryRange(neighborId, queryrange);
                    // if (!inrange)
                    //     prob = ProbFunc(next_step);
                    // if (!prob)
                    //     continue;

                    selected_edges.emplace_back(neighborId, inrange);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }
            } while (cur_node->lbound < ql || cur_node->rbound > qr);
            return selected_edges;
        }

        std::vector<tableint> SelectEdgeOr(int pid, int ql, int qr, int edge_limit, searcher::Bitset<uint64_t> &visited_set)
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

        std::priority_queue<PFI> TopDown_nodeentries_search(std::vector<TreeNode *> &filterednodes, const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit)
        {
            // To fix the starting points for different 'ef' parameter, set seed to a fixed number, e.g., seed =0
            // unsigned seed = 0;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFI, std::vector<PFI>, std::greater<PFI>> candidate_set;
            std::priority_queue<PFI> top_candidates;
            searcher::Bitset<uint64_t> visited_set(max_elements_);

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(dis, pid);
                top_candidates.emplace(dis, pid);
            }

            float lowerBound = top_candidates.top().first;

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                ++metric_hops;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second;
                auto selected_edges = SelectEdgeOr(current_pid, QL, QR, edge_limit, visited_set);
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
                        candidate_set.emplace(dis, neighbor_id);
                        top_candidates.emplace(dis, neighbor_id);
                        lowerBound = top_candidates.top().first;
                    }
                    else if (dis < lowerBound)
                    {
                        candidate_set.emplace(dis, neighbor_id);
                        top_candidates.emplace(dis, neighbor_id);
                        top_candidates.pop();
                        lowerBound = top_candidates.top().first;
                    }
                }
            }
            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }
        std::priority_queue<PFI> Postfilter(const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit, std::pair<int, int> queryrange, std::vector<iRangeGraph::TreeNode *> &filterednodes)
        {
            searcher::Bitset<uint64_t> visited_set(max_elements_);
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFI, std::vector<PFI>, std::greater<PFI>> candidate_set;
            std::priority_queue<PFI> top_candidates;

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                    continue;
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(dis, pid);
                if (CheckInQueryRange(pid, queryrange))
                    top_candidates.emplace(dis, pid);
            }

            float lowerBound = std::numeric_limits<float>::max();

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                metric_hops++;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second;
                auto selected_edges = SelectEdge(current_pid, QL, QR, edge_limit, queryrange);

                while (selected_edges.size())
                {
                    auto neighbor_pair = selected_edges.back();
                    selected_edges.pop_back();
                    int neighbor_id = neighbor_pair.first;
                    bool inrange = neighbor_pair.second;

                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    metric_distance_computations++;

                    if (top_candidates.size() < ef || dis < lowerBound)
                    {
                        if (inrange)
                        {
                            top_candidates.emplace(dis, neighbor_id);
                        }
                        candidate_set.emplace(dis, neighbor_id);

                        if (top_candidates.size() > ef)
                            top_candidates.pop();
                        if (top_candidates.size())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }

            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }

        void searchfilter(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            for (auto range : storage->query_range)
            {
                std::string domain = range.first;
                std::vector<std::vector<int>> &gt = storage->ground_truth[domain];
                std::string savepath = saveprefix + domain + "-lex-postfilter.csv";
                CheckPath(savepath);

                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {

                        auto cons = range.second[i];
                        int ql = 0;
                        int qr = 999999;

                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        auto filterednodes = tree->range_filter(tree->root, ql, qr);
                        auto res = Postfilter(storage->query_points[i].data(), ef, storage->query_K, ql, qr, edge_limit, cons, filterednodes);
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

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
                }

                outfile.close();
            }
        }


        int MaxStep{20};
        std::vector<double> probability;
        void setprob()
        {
            probability.resize(MaxStep);

            for (int x = 0; x < MaxStep; x++)
            {
                probability[x] = 1 / (1 + exp(x));
            }
        }

        // purepost = True -> p=1   purepost =  False -> 0<=p<=1
        bool purepost{false};

        int ProbFunc(int x)
        {
            if (purepost)
                return 1;
            if (x >= MaxStep)
                return 0;
            double randNum = (double)rand() / RAND_MAX;
            return randNum < probability[x] ? 1 : 0;
        }
        std::vector<bool> flags = std::vector<bool>(10000000, false);
        std::vector<std::pair<tableint, bool>> SelectEdgeF(int pid, int ql, int qr, int edge_limit, std::pair<int, int> &queryrange, int current_step)
        {
            iRangeGraph::TreeNode *cur_node = nullptr, *nxt_node = tree->root;
            std::vector<std::pair<tableint, bool>> selected_edges;
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
                        for (int i = 0; i < cur_node->childs.size(); i++)
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

                for (size_t j = 1; j <= size; j++)
                {
                    int neighborId = *(data + j);
                    if (neighborId < ql || neighborId > qr)
                        continue;
                    int prob = 1;
                    int next_step = current_step + 1;
                    bool inrange = CheckInQueryRange(neighborId, queryrange);
                    if (!inrange)
                        prob = ProbFunc(next_step);
                    if (!prob)
                        continue;

                    selected_edges.emplace_back(neighborId, inrange);
                    if (selected_edges.size() == edge_limit)
                        return selected_edges;
                }
            } while (cur_node->lbound < ql || cur_node->rbound > qr);
            return selected_edges;
        }

        std::priority_queue<PFI> Postfilter_P(const void *query_data, int ef, int query_k, int QL, int QR, int edge_limit, std::pair<int, int> queryrange, std::vector<iRangeGraph::TreeNode *> &filterednodes)
        {
            searcher::Bitset<uint64_t> visited_set(max_elements_);
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine e(seed);

            std::priority_queue<PFII, std::vector<PFII>, std::greater<PFII>> candidate_set;
            std::priority_queue<PFI> top_candidates;

            for (auto u : filterednodes)
            {
                std::uniform_int_distribution<int> u_start(u->lbound, u->rbound);
                int pid = u_start(e);
                if (visited_set.get(pid))
                    continue;
                visited_set.set(pid);
                char *ep_data = getDataByInternalId(pid);
                float dis = fstdistfunc_(query_data, ep_data, dist_func_param_);
                candidate_set.emplace(std::make_pair(dis, std::make_pair(pid, -1)));
                if (CheckInQueryRange(pid, queryrange))
                    top_candidates.emplace(dis, pid);
            }

            float lowerBound = std::numeric_limits<float>::max();

            while (!candidate_set.empty())
            {
                auto current_point_pair = candidate_set.top();
                metric_hops++;
                if (current_point_pair.first > lowerBound)
                {
                    break;
                }
                candidate_set.pop();
                int current_pid = current_point_pair.second.first;
                int current_step = current_point_pair.second.second;
                auto selected_edges = SelectEdgeF(current_pid, QL, QR, edge_limit, queryrange, current_step);

                while (selected_edges.size())
                {
                    auto neighbor_pair = selected_edges.back();
                    selected_edges.pop_back();
                    int neighbor_id = neighbor_pair.first;
                    bool inrange = neighbor_pair.second;

                    if (visited_set.get(neighbor_id))
                        continue;
                    visited_set.set(neighbor_id);
                    char *neighbor_data = getDataByInternalId(neighbor_id);
                    float dis = fstdistfunc_(query_data, neighbor_data, dist_func_param_);
                    metric_distance_computations++;

                    if (top_candidates.size() < ef || dis < lowerBound)
                    {
                        int next_step = current_step + 1;
                        if (inrange)
                        {
                            top_candidates.emplace(dis, neighbor_id);
                            next_step = -1;
                        }
                        candidate_set.emplace(std::make_pair(dis, std::make_pair(neighbor_id, next_step)));

                        if (top_candidates.size() > ef)
                            top_candidates.pop();
                        if (top_candidates.size())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }

            while (top_candidates.size() > query_k)
                top_candidates.pop();
            return top_candidates;
        }
        void searchfilterP(std::vector<int> &SearchEF, std::string saveprefix, int edge_limit)
        {
            for (auto range : storage->query_range)
            {
                std::string domain = range.first;
                std::vector<std::vector<int>> &gt = storage->ground_truth[domain];
                std::string savepath = saveprefix + domain + "-lex-postfilter-p.csv";
                CheckPath(savepath);

                std::ofstream outfile(savepath);
                if (!outfile.is_open())
                {
                    throw Exception("cannot open " + savepath);
                }
                for (auto ef : SearchEF)
                {
                    int tp = 0;
                    float searchtime = 0;

                    metric_hops = 0;
                    metric_distance_computations = 0;

                    for (int i = 0; i < storage->query_nb; i++)
                    {

                        auto cons = range.second[i];
                        int ql = 0;
                        int qr = 999999;

                        timeval t1, t2;
                        gettimeofday(&t1, NULL);
                        auto filterednodes = tree->range_filter(tree->root, ql, qr);
                        std::priority_queue<PFI> res = Postfilter_P(storage->query_points[i].data(), ef, storage->query_K, ql, qr, edge_limit, cons, filterednodes);
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

                    float recall = 1.0 * tp / storage->query_nb / storage->query_K;
                    float qps = storage->query_nb / searchtime;
                    float dco = metric_distance_computations * 1.0 / storage->query_nb;
                    float hop = metric_hops * 1.0 / storage->query_nb;

                    outfile << ef << "," << recall << "," << qps << "," << dco << "," << hop << std::endl;
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