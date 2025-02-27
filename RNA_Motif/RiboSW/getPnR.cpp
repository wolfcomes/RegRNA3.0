using namespace std;

const float pair_ratio   = 2/3;

resultOfPair getRnP (const string operate, const shape_inf mySI, const int startP, const string seq_before_op, const string seq_after_op){
    resultOfPair myRop;
    myRop.conti = false;

    vector < int > sec_struct(operate.size(), '.');  // Recording 2nd structure

    if (seq_after_op.size() > mySI.wte_max) { myRop.conti = true; return myRop; }
    if (operate.size() < (mySI.stemLower*2+1)) { myRop.conti = true; return myRop; }

/*for (int pos=mySI.stemLower; pos<=mySI.stemUpper+1; pos++) {
    string xxx = operate.substr(0, operate.size()-pos);
    if (! pattern_match("au$|ac$|gc$|uc$", xxx) ) { myRop.conti = true; return myRop; }
}*/

    string stem_h5 = operate.substr(0, mySI.stemUpper);
    string stem_h3 = operate.substr(operate.size()-mySI.stemUpper, mySI.stemUpper);
    reverse(stem_h3.begin(), stem_h3.end());    // reverse helix 3 prime of stem for pairing
    vector < vector <int> > arrayPaired;
    arrayPaired = searchStemPair(stem_h5, stem_h3, mySI.stemUpper, mySI.stemLower, operate.length(), mySI.bulge);

    // check pair number - threshold
    int stemPairNo = (int)arrayPaired.size();
    if (stemPairNo > mySI.stemUpper || stemPairNo < mySI.stemLower) { myRop.conti = true; return myRop; }

    // check pairing is not over head
    int numOfLoopNuc = arrayPaired[arrayPaired.size()-1][2] - arrayPaired[arrayPaired.size()-1][0] -1;
    if ( (numOfLoopNuc > mySI.loopUpper) || (numOfLoopNuc < mySI.loopLower) ) { myRop.conti = true; return myRop; }
    string loop_interval = "";
    loop_interval = operate.substr(arrayPaired[arrayPaired.size()-1][0]+1, numOfLoopNuc);
    
    // check for pairing number, if pairing number is less than half of stem length, refuse it
    int pairNumOfStem = (int)arrayPaired.size();
    if ( pairNumOfStem < (mySI.stemUpper/2) ) { myRop.conti = true; return myRop; }
    float min_pair = pairRatio ( ((int)arrayPaired[arrayPaired.size()-1][0])-((int)arrayPaired[0][0]), ((int)arrayPaired[0][2])-                     ((int)arrayPaired[arrayPaired.size()-1][2]), pair_ratio );
    if ( pairNumOfStem < min_pair ) { myRop.conti = true; return myRop; }

/*/ check pair number - threshold
int stemPairNo = (int)arrayPaired.size();
if (stemPairNo > mySI.stemUpper || stemPairNo < mySI.stemLower) { myRop.conti = true; return myRop; }*/

    // saving pairing information
    string h5_pairing_seq = "", h3_pairing_seq = "";
    for (int i=0; i<arrayPaired.size(); i++) {
        h5_pairing_seq.push_back(arrayPaired[i][1]);
        h3_pairing_seq.push_back(arrayPaired[i][3]);    // h3 seq is from 3' to 5'
    }
    
//string loop_interval = "";
//loop_interval = operate.substr(arrayPaired[arrayPaired.size()-1][0]+1, numOfLoopNuc);

    // 2nd structure creating
    //sec_struct.push_back((int)(lastPairedNoIn3prime-lengthOfStem1Interval+4));
    // ABOVE IS FOR SAVING NEXT PAIRING SEARCH START POINT
    for (int c=0; c<(int)arrayPaired.size(); c++) {     // Creating 2nd structure
        sec_struct[ (int)arrayPaired[c][0] ] = '(';
        sec_struct[ (int)arrayPaired[c][2] ] = ')';
    }

    string outh5 = operate.substr(0, arrayPaired[0][0]);
    string outh3 = operate.substr(arrayPaired[0][2]+1);
    if (seq_after_op.size()+outh3.size() > mySI.wte_max) { myRop.conti = true; return myRop; }
    if ((mySI.wte_min > 0) && (seq_after_op.size()+outh3.size() < mySI.wte_min)) { myRop.conti = true; return myRop; }

    myRop.arrayPaired       = arrayPaired;
    myRop.sec_struct        = sec_struct;
    myRop.h5_pairing_seq    = h5_pairing_seq;
    myRop.h3_pairing_seq    = h3_pairing_seq;
    myRop.loop_interval     = loop_interval;
    myRop.outer_h5          = seq_before_op + outh5;
    myRop.outer_h3          = outh3 + seq_after_op;

    // Checking if sequence pattern is corresponding provided restrict
    if ( (mySI.seq_info_descr != "") && !(report_seq_info(myRop, mySI.seq_info_descr)) ) { myRop.conti = true; return myRop; }

    return myRop;
}
