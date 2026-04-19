# iRangeGraph_Muti-Attribute


## Quick Start

### Build

```bash
mkdir build && cd build && cmake .. && make
```

##MIANNS

### Construct Index

#### parameters:

**`--data_path`**: The input data over which to build an index, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time. 
The data points should be already sorted in ascending order by the attribute.

**`--index_file`**: The constructed index will be saved to this file, in .bin format.

**`--M`**: The degree of the graph index.

**`--ef_construction`**: The size of result set during index building.

**`--threads`**: The number of threads for index building.


#### command:
```bash
./tests/buildindex --data_path [path to data points] --index_file [file path to save index] --M [integer] --ef_construction [integer] --threads [integer]
```


### Search For Muti-Interval


#### parameters:

**`--data_path`**: The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
The data points should be already sorted in ascending order by the attribute.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix`**: The path of folder where query range files will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved, in .bin format. 

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index.

**`--interval_num`**: the number of intervals.

**`--query_K`**: TOP-k.

#### command:
```bash
./tests/search_multi_interval --data_path [path to data points] --query_path [path to query points] --range_saveprefix [folder path to save query ranges] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --M [integer] --interval_num [integer] --query_K [integer]
```



##MAANNS

Please note that the current code implementation of MAANNS is restricted to two attributes

### dataset Sort

The data needs to be pre-sorted according to your dataset.

#### command:
```bash
./tests/lex
./tests/zorder
```

### Construct Index

#### parameters:

**`--data_path`**: The input data over which to build an index, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time. 
The data points should be already sorted in ascending order by the attribute.

**`--index_file`**: The constructed index will be saved to this file, in .bin format.

**`--M`**: The degree of the graph index.

**`--ef_construction`**: The size of result set during index building.

**`--threads`**: The number of threads for index building.


#### command:
```bash
./tests/buildindex --data_path [path to data points] --index_file [file path to save index] --M [integer] --ef_construction [integer] --threads [integer]
```



### Search For Multi-Attribute


Theoretically, the Z-order implementation and the predicate logic (supporting arbitrary combinations of AND/OR) could be unified. However, for the sake of simplicity and modularity, we have implemented them separately.

#### z-order-and:

#### parameter:
**`--data_path`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
There is no need to pre-sort the data points by any attribute. Just make sure data points and attribute1 and attribute2 match one by one in order.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix`**: The path of folder where query range files(.csv) will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved.

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--attribute1_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--attribute2_file`**: The path of the second attribute file, in .bin format. `n*sizeof(int)` bytes contain the second attributes of the data for one data point in a time.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index by the first attribute.

**`--query_K`**: Top-K.


#### command:
```bash
./tests/search_multi_attribute_zorder_and --data_path [path to data points] --query_path [path to query points] --range_saveprefix [folder path to save query ranges] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --attribute1 [path to first attributes] --attribute2 [path to second attributes] --M [integer] --query_K [integer]
```


#### z-order-or:

#### parameter:
**`--data_path`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
There is no need to pre-sort the data points by any attribute. Just make sure data points and attribute1 and attribute2 match one by one in order.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix1`**: The path of folder where the first attributes of query range files(.csv) will be saved.

**`--range_saveprefix2`**: The path of folder where the second attributes of query range files(.csv) will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved.

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--attribute1_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--attribute2_file`**: The path of the second attribute file, in .bin format. `n*sizeof(int)` bytes contain the second attributes of the data for one data point in a time.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index by the first attribute.

**`--query_K`**: Top-K.


#### command:
```bash
./tests/search_multi_attribute_zorder_or --data_path [path to data points] --query_path [path to query points] --range_saveprefix1 [folder path to save query ranges of first attribute] --range_saveprefix2 [folder path to save query ranges of second attribute] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --attribute1 [path to first attributes] --attribute2 [path to second attributes] --M [integer] --query_K [integer]
```


#### lex-and:

#### parameter:
**`--data_path`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
There is no need to pre-sort the data points by any attribute. Just make sure data points and attribute1 and attribute2 match one by one in order.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix`**: The path of folder where query range files(.csv) will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file`**: The file path where the constructed index is saved.

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--attribute1_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--attribute2_file`**: The path of the second attribute file, in .bin format. `n*sizeof(int)` bytes contain the second attributes of the data for one data point in a time.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index by the first attribute.

**`--query_K`**: Top-K.


#### command:
```bash
./tests/search_multi_attribute_lex_and --data_path [path to data points] --query_path [path to query points] --range_saveprefix [folder path to save query ranges] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --attribute1 [path to first attributes] --attribute2 [path to second attributes] --M [integer] --query_K [integer]
```


#### lex-or:

#### parameter:

**`--data_path1`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
There is no need to pre-sort the data points by any attribute. Just make sure data points and attribute11 and attribute12 match one by one in order.

**`--data_path2`**:  The data points over which the index is built, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the data one data point in a time.
There is no need to pre-sort the data points by any attribute. Just make sure data points and attribute21 and attribute22 match one by one in order.

**`--query_path`**: The query vectors, in .bin format. The first 4 bytes represent number of points as integer. The next 4 bytes represent the dimension of data as integer. The following `n*d*sizeof(float)` bytes contain the contents of the query one query point in a time.

**`--range_saveprefix1`**: The path of folder where the first attributes of query range files(.csv) will be saved.

**`--range_saveprefix2`**: The path of folder where the second attributes of query range files(.csv) will be saved.

**`--groundtruth_saveprefix`**: The path of folder where groundtruth files will be saved.

**`--index_file1`**: The file path where the constructed index1 is saved.

**`--index_file2`**: The file path where the constructed index2 is saved.

**`--result_saveprefix`**: The path of folder where result files will be saved.

**`--attribute11_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--attribute12_file`**: The path of the second attribute file, in .bin format. `n*sizeof(int)` bytes contain the second attributes of the data for one data point in a time.

**`--attribute21_file`**: The path of the first attribute file, in .bin format. `n*sizeof(int)` bytes contain the first attributes of the data for one data point in a time.

**`--vid1`**: vid in index1.

**`--vid2`**: vid in index2.

**`--M`**: The degree of the graph index. It should equal the 'M' used for constructing index by the first attribute.

**`--query_K`**: Top-K.


#### command:
```bash
./tests/search_multi_attribute_zorder_or --data_path1 [path to data points in index1] --data_path2 [path to data points in index2] --query_path [path to query points] --range_saveprefix1 [folder path to save query ranges of first attribute] --range_saveprefix2 [folder path to save query ranges of second attribute] --groundtruth_saveprefix [folder path to save groundtruth] --index_file [path of the index file] --result_saveprefix [folder path to save results] --attribute11 [path to first attributes in index1] --attribute12 [path to second attributes in index2] --attribute21 [path to first attributes in index2] --attribute22 [path to second attributes in index2] --vid1 [path to vid in index1] --vid2 [path to vid in index2] --M [integer] --query_K [integer]
```

## MA Datasets
| Dataset |Vector Type| Dimension | Attribute Type |
|---------|-----------|-----------|----------------|
|   [WIT](https://github.com/google-research-datasets/wit)   |   image   |   2048    |  height, width   |



