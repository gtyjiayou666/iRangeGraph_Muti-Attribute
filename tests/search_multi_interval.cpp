#include "iRG_search_multi_interval.h"

std::unordered_map<std::string, std::string> paths;

int query_K;
int interval_num;
int M;

void Generate(iRangeGraph::DataLoaderMultiInterval &storage)
{
    storage.LoadData(paths["data_vector"]);
    iRangeGraph::QueryGeneratorMultiInterval generator(storage.data_nb, storage.query_nb, interval_num);
    generator.GenerateMultiInterval(paths["range_saveprefix"], 1000000);
    storage.LoadQueryRange(paths["range_saveprefix"]);
    generator.GenerateGroundtruthThread(paths["groundtruth_saveprefix"], storage);
    // generator.GenerateMultiIntervalLargeGaps(paths["range_saveprefix"]);
    // storage.LoadQueryRangeLargeGaps(paths["range_saveprefix"]);
    // generator.GenerateGroundtruthThreadLargeGaps(paths["groundtruth_saveprefix"], storage);
}

void init()
{
    // data vectors should be sorted by the attribute values in ascending order
    paths["data_vector"] = "";

    paths["query_vector"] = "";
    // the path of document where range files are saved
    paths["range_saveprefix"] = "";
    // the path of document where groundtruth files are saved
    paths["groundtruth_saveprefix"] = "";
    // the path where index file is saved
    paths["index"] = "";
    // the path of document where search result files are saved
    paths["result_saveprefix"] = "";
    paths["interval_num"] = "";
    paths["query_K"] = "";
    // M is the maximum out-degree same as index build
}

int main(int argc, char **argv)
{
    // init();

    for (int i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "--data_path")
            paths["data_vector"] = argv[i + 1];
        if (arg == "--query_path")
            paths["query_vector"] = argv[i + 1];
        if (arg == "--range_saveprefix")
            paths["range_saveprefix"] = argv[i + 1];
        if (arg == "--groundtruth_saveprefix")
            paths["groundtruth_saveprefix"] = argv[i + 1];
        if (arg == "--index_file")
            paths["index"] = argv[i + 1];
        if (arg == "--result_saveprefix")
            paths["result_saveprefix"] = argv[i + 1];
        if (arg == "--M")
            M = std::stoi(argv[i + 1]);
        if (arg == "--interval_num")
            interval_num = std::stoi(argv[i + 1]);
        if (arg == "--query_K")
            query_K = std::stoi(argv[i + 1]);
    }

    if (argc != 19)
        throw Exception("please check input parameters");

    iRangeGraph::DataLoaderMultiInterval storage(interval_num);
    storage.query_K = query_K;
    storage.LoadQuery(paths["query_vector"]);
    // If it is the first run, Generate shall be called; otherwise, Generate can be skipped
    Generate(storage);
    storage.LoadQueryRange(paths["range_saveprefix"]);
    storage.LoadGroundtruth(paths["groundtruth_saveprefix"]);
    // storage.LoadQueryRangepro(paths["range_saveprefix"]);
    // storage.LoadGroundtruthpro(paths["groundtruth_saveprefix"]);

    iRangeGraph::iRangeGraph_Search_Muti_Interval<float> index(paths["data_vector"], paths["index"], &storage, M);
    std::vector<int> SearchEF = { 400, 350, 300, 250, 200, 180, 160, 140, 120, 100, 90, 80, 70, 60, 55, 45, 40, 35, 30, 25, 20, 15, 10};
    index.setprob();
    index.search_Prefilter(SearchEF, paths["result_saveprefix"], M);
    index.search_Postfilter(SearchEF, paths["result_saveprefix"], M);
    index.search_MIDG(SearchEF, paths["result_saveprefix"], M);
    index.search_MIDG_P(SearchEF, paths["result_saveprefix"], M);
    index.searchor_MIDG_G(SearchEF, paths["result_saveprefix"], M);
    index.search_Postfilter_P(SearchEF, paths["result_saveprefix"], M);
}
