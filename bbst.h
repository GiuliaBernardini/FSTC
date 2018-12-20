/* Taken from https://github.com/kowallus/BbST under GNU GPL*/

#ifndef SBRMA2_NOC_H
#define SBRMA2_NOC_H

#define MAX_T_VALUE INT64_MAX
#define MAX_T_ARRAYSIZE INT64_MAX
typedef long int t_value;
typedef long int t_array_size;

#include <vector>

using namespace std;

class BbST {
public:
#ifdef MINI_BLOCKS
    BbST(const t_value* valuesArray, const t_array_size n, int kExp, int miniKExp);
    BbST(int kExp, int miniKExp);
#else 
    BbST(const t_value* valuesArray, const t_array_size n, int kExp);
    BbST(int kExp);
#endif
    void rmqBatch(const t_value* valuesArray, const t_array_size n, const vector<t_array_size> &queries, t_array_size *resultLoc);
    void rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc);
    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);

    virtual ~BbST();

    size_t memUsageInBytes();

private:
    t_array_size blocksCount;
    int k, kExp, D, miniK, miniKExp;

    const t_value *valuesArray;
    t_array_size n;

    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    t_array_size miniBlocksCount;
    int miniBlocksInBlock;
    uint8_t* miniBlocksLoc = 0;

    void getBlocksMinsBase();
    void getBlocksSparseTable();

    t_array_size scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& smallerThanVal, bool orEqual);
    inline t_array_size rawScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& smallerThanVal, bool orEqual);
    inline t_array_size miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& smallerThanVal, bool orEqual);

    bool batchMode = false;
    void cleanup();

};

#endif //SBRMA2_NOC_H
