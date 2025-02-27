using namespace std;

struct reportStruct {
    int hPlus;
    int hMinus; 
    int hFull;
    int oPlus;
    int oMinus; 
    int oFull;
    int lPlus;
    int lMinus; 
    int lFull;
};

struct reportSEQ {
    string h5;
    string h3;
    string loop;
};

struct shape_inf {
    int swL;
    int swU;
    int swStep;
    int spStep;
    int stemLower;
    int stemUpper;
    int loopLower;
    int loopUpper;
    int bulge;
    int pushNo;
    string report_descr;
    string seq_info_descr;
    int wts_min;
    int wts_max;
    int wte_min;
    int wte_max;
};

struct resultOfPair {
    bool conti;
    vector < vector <int> > arrayPaired;
    vector <int> sec_struct;
    string h5_pairing_seq;
    string h3_pairing_seq;
    string loop_interval;
    string outer_h5;
    string outer_h3;
};
