using namespace std;

inline float pairRatio (int n1, int n2, float pair_ratio) {
    if (n2 > n1) {
        return float (n2 * pair_ratio);
    } else {
        return float (n1 * pair_ratio);
    }
}

