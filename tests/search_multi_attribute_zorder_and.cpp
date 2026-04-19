
#include "iRG_search_multi.h"
#include "iRG_search_zorder.h"

std::unordered_map<std::string, std::string> paths;


void Generate(iRangeGraph_multi::DataLoader &storage)
{
    storage.LoadRanges(paths["range_prefix"]);
    storage.Generate_Groundtruth(paths["groundtruth_prefix"]);
}
int query_K = 10;
int M;

int main(int argc, char **argv)
{
    for (int i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--data_path")
            paths["data_vector"] = argv[i + 1];
        if (arg == "--query_path")
            paths["query_vector"] = argv[i + 1];
        if (arg == "--range_saveprefix")
            paths["range_prefix"] = argv[i + 1];
        if (arg == "--groundtruth_saveprefix")
            paths["groundtruth_prefix"] = argv[i + 1];
        if (arg == "--index_file")
            paths["index"] = argv[i + 1];
        if (arg == "--result_saveprefix")
            paths["result_saveprefix"] = argv[i + 1];
        if (arg == "--attribute1")
            paths["attribute1"] = argv[i + 1];
        if (arg == "--attribute2")
            paths["attribute2"] = argv[i + 1];
        if (arg == "--M")
            M = std::stoi(argv[i + 1]);
        if (arg == "--query_K")
            query_K = std::stoi(argv[i + 1]);
    }

    iRangeGraph_multi::DataLoader storage;
    storage.query_K = query_K;
    storage.LoadQuery(paths["query_vector"]);
    storage.LoadData(paths["data_vector"]);
    // the order of load in attribute values should not be switched
    // the first attribute should be first loaded
    storage.LoadAttribute(paths["attribute1"]);
    storage.LoadAttribute(paths["attribute2"]);
    storage.CalZorder();


    Generate(storage);
    storage.LoadRanges(paths["range_prefix"]);
    // storage.Generate_Groundtruth(paths["groundtruth_prefix"]);

    storage.LoadGroundtruth(paths["groundtruth_prefix"]);

    iRangeGraph::iRangeGraph_Search_Zorder<float> index(paths["index"], &storage, M);
    std::vector<int> SearchEF = {300, 250, 200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10};
    index.searchand(SearchEF, paths["result_saveprefix"], M);
}