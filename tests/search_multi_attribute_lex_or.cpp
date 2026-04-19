
#include "iRG_search_multi.h"
#include "iRG_search_lex.h"

std::unordered_map<std::string, std::string> paths;

int query_K = 10;
int M;

void Generate(iRangeGraph_multi::DataLoader &storage)
{
    storage.LoadRangesOr(paths["range_prefix1"], paths["range_prefix2"]);
    storage.Generate_GroundtruthOr(paths["groundtruth_prefix"]);
}

int main(int argc, char **argv)
{
    for (int i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--data_path1")
            paths["data_vector1"] = argv[i + 1];
        if (arg == "--data_path2")
            paths["data_vector2"] = argv[i + 1];
        if (arg == "--query_path")
            paths["query_vector"] = argv[i + 1];
        if (arg == "--range_saveprefix1")
            paths["range_prefix1"] = argv[i + 1];
        if (arg == "--range_saveprefix2")
            paths["range_prefix2"] = argv[i + 1];
        if (arg == "--groundtruth_saveprefix")
            paths["groundtruth_prefix"] = argv[i + 1];
        if (arg == "--index_file1")
            paths["index1"] = argv[i + 1];
        if (arg == "--index_file2")
            paths["index2"] = argv[i + 1];
        if (arg == "--result_saveprefix")
            paths["result_saveprefix"] = argv[i + 1];
        if (arg == "--attribute11")
            paths["attribute11"] = argv[i + 1];
        if (arg == "--attribute12")
            paths["attribute12"] = argv[i + 1];
        if (arg == "--attribute21")
            paths["attribute21"] = argv[i + 1];
        if (arg == "--attribute22")
            paths["attribute22"] = argv[i + 1];
        if (arg == "--vid1")
            paths["vid1"] = argv[i + 1];
        if (arg == "--vid2")
            paths["vid2"] = argv[i + 1];
        if (arg == "--M")
            M = std::stoi(argv[i + 1]);
        if (arg == "--query_K")
            query_K = std::stoi(argv[i + 1]);
    }

    iRangeGraph_multi::DataLoader storage1;
    storage1.query_K = query_K;
    storage1.LoadQuery(paths["query_vector"]);
    storage1.LoadData(paths["data_vector1"]);
    storage1.LoadAttribute(paths["attribute11"]);
    storage1.LoadAttribute(paths["attribute12"]);
    storage1.LoadKey(paths["attribute11"]);
    storage1.LoadVid(paths["vid1"]);

    iRangeGraph_multi::DataLoader storage2;
    storage2.query_K = query_K;
    storage2.LoadQuery(paths["query_vector"]);
    storage2.LoadData(paths["data_vector2"]);
    storage2.LoadAttribute(paths["attribute21"]);
    storage2.LoadAttribute(paths["attribute22"]);
    storage2.LoadKey(paths["attribute22"]);
    storage2.LoadVid(paths["vid2"]);

    // Generate should be called when running for the first time; otherwise, it can be skipped.
    Generate(storage1);
    storage1.LoadRangesOr(paths["range_prefix1"], paths["range_prefix2"]);
    storage1.LoadGroundtruthOr(paths["groundtruth_prefix"]);
    storage2.LoadRangesOr(paths["range_prefix1"], paths["range_prefix2"]);
    storage2.LoadGroundtruthOr(paths["groundtruth_prefix"]);
    iRangeGraph::iRangeGraph_Search_Lex<float> index1(paths["index1"], &storage1, M);
    iRangeGraph::iRangeGraph_Search_Lex<float> index2(paths["index2"], &storage2, M);
    index1.setprob();
    index2.setprob();
    std::vector<int> SearchEF = {300, 250, 200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10};

    for (auto range : storage1.query_range_or)
    {
        std::string domain = range.first;
        std::vector<std::vector<int>> &gt = storage1.ground_truth[domain];
        std::string savepath = paths["result_saveprefix"] + domain + "-lex-or.csv";
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
            auto &col_data1 = storage1.keys;

            auto &col_data2 = storage2.keys;

            for (int i = 0; i < storage1.query_nb; i++)
            {
                Attr_Constraint cons = range.second[i];
                int x_1 = cons.attr_constraints[0].first;
                int x_2 = cons.attr_constraints[0].second;
                int y_1 = cons.attr_constraints[1].first;
                int y_2 = cons.attr_constraints[1].second;

                timeval t1, t2;
                gettimeofday(&t1, NULL);
                int l1, r1, l2, r2;
                auto it_start1 = std::lower_bound(col_data1.begin(), col_data1.end(), x_1);
                auto it_end1 = std::upper_bound(col_data1.begin(), col_data1.end(), x_2);
                l1 = it_start1 - col_data1.begin();
                r1 = it_end1 - col_data1.begin() - 1;
                auto filterednodes1 = index1.tree->range_filter(index1.tree->root, l1, r1);
                auto res1 = index1.TopDown_nodeentries_search(filterednodes1, storage1.query_points[i].data(), ef, storage1.query_K, l1, r1, 16);
                auto it_start2 = std::lower_bound(col_data2.begin(), col_data2.end(), y_1);
                auto it_end2 = std::upper_bound(col_data2.begin(), col_data2.end(), y_2);
                l2 = it_start2 - col_data2.begin();
                r2 = it_end2 - col_data2.begin() - 1;
                auto filterednodes2 = index2.tree->range_filter(index2.tree->root, l2, r2);
                auto res2 = index2.TopDown_nodeentries_search(filterednodes2, storage2.query_points[i].data(), ef, storage2.query_K, l2, r2, 16);
                std::vector<PFI> final_results;
                final_results.reserve(res1.size() + res2.size());
                while (!res1.empty())
                {
                    auto x = res1.top();
                    x.second = storage1.vids[x.second];
                    final_results.push_back(x);
                    res1.pop();
                }
                while (!res2.empty())
                {
                    auto x = res2.top();
                    x.second = storage2.vids[x.second];
                    final_results.push_back(x);
                    res2.pop();
                }
                std::sort(final_results.begin(), final_results.end(), [](const PFI &a, const PFI &b)
                          { return a.second < b.second; });
                auto last = std::unique(final_results.begin(), final_results.end(), [](const PFI &a, const PFI &b)
                                        {
                                            return a.second == b.second;
                                        });
                final_results.erase(last, final_results.end());
                std::sort(final_results.begin(), final_results.end(), [](const PFI &a, const PFI &b)
                          { return a.first < b.first; });
                gettimeofday(&t2, NULL);
                searchtime += GetTime(t1, t2);

                std::map<int, int> record;

                for (int x = 0; x < storage1.query_K; x++)
                {
                    auto y = final_results[x].second;
                    if (record.count(y))
                        throw Exception("repetitive search results");
                    record[y] = 1;
                    if (std::find(gt[i].begin(), gt[i].end(), y) != gt[i].end())
                        tp++;
                }
            }

            float recall = 1.0 * tp / storage1.query_nb / storage1.query_K;
            float qps = storage1.query_nb / searchtime;

            outfile << ef << "," << recall << "," << qps << std::endl;
        }

        outfile.close();
    }
}