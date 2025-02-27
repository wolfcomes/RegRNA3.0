#include "include.h"
#include <stdexcept>

const int DEBUG_struct = false;
const int DEBUG = false;
const int STRUCT_LIMIT = 100;

using namespace std;
using namespace pcrecpp;

// subfunction commence
#include "c_def_struct.h"
#include "prototype.h"
#include "seq_pattern.cpp"
#include "stempair.cpp"
#include "getPairRatio.cpp"
#include "tokenSplit.cpp"
#include "getPnR.cpp"
#include "seq_info.cpp"
// subfunction fini
/*
inline void searchStem_start(const string seqStr, const int);
inline void searchStem_gt2(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
inline void searchStem_gt3(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
inline void searchStem_gt4(const string seqStr, const int);
*/
void searchStem_start(const string seqStr, const int);
//void searchStem_gt2(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
void searchStem_gt2(vector <string> const &, vector <string> const &, string const &, const string &, vector <int> const &, const int &);
//void searchStem_gt3(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
void searchStem_gt3(vector <string> const &, vector <string> const &, string const &, const string &, vector <int> const &, const int &);
void searchStem_gt4(const string seqStr, const int);

//const float pair_ratio   = 2/3;

int globalStartPoint    = 0;
int baseOfGSP           = 0;
int struct_counter      = 0;
int seq_size = 0;


string complete_seq;

string fasta;

/*
 * Nouveau variables commence ici
 */

int global_depth = 0;   // Global depth, for result reporting
int startPstruct =0;    // control structures of each startP
vector <string> v_descr;

/*
 * Nouveau variables fini ici
 */

template <class T>
T from_string(  const std::string& s,
                std::ios_base& (*f) (std::ios_base&))
{
    std::istringstream iss(s);
    //return !(iss >> f >> t).fail();
    T t;
    (iss >> f >> t).fail();
    return t;
}
void seq_checking (reportSEQ myReportSeq) {
    cout << "h5: '" << myReportSeq.h5 << "'\th3: '" << myReportSeq.h3 << "'\tloop: '" << myReportSeq.loop << "'\t now seq. saved" << endl;
}
void report_checking (reportStruct myRS) {

    cout << "IN REPORTER\n";

        if (myRS.lPlus!=-1) {
            cout << "GOT loop plus, value= " << myRS.lPlus << endl;
        } if (myRS.lMinus!=-1) {
            cout << "GOT loop minus, value= " << myRS.lMinus << endl;
        } if (myRS.lFull!=-1) {
            cout << "GOT loop full, value= " << myRS.lFull << endl;
        } if (myRS.hPlus!=-1) {
            cout << "GOT helix plus, value= " << myRS.hPlus << endl;
        } if (myRS.hMinus!=-1) {
            cout << "GOT helix minus, value= " << myRS.hMinus << endl;
        } if (myRS.hFull!=-1) {
            cout << "GOT helix full, value= " << myRS.hFull << endl;
        } if (myRS.oPlus!=-1) {
            cout << "GOT outer plus, value= " << myRS.oPlus << endl;
        } if (myRS.oMinus!=-1) {
            cout << "GOT outer minus, value= " << myRS.oMinus << endl;
        } if (myRS.oFull!=-1) {
            cout << "GOT outer full, value= " << myRS.oFull << endl;
        }
}

void report_pairing_result4debug (string operate, vector < vector <int> > arrayPaired, vector <int> sec_struct) {
    //*     ^^^FOR DEBUG^^^
    if (true) {
        cout << "operate: " << operate << endl;
        for ( int ii=0; ii<arrayPaired.size(); ii++ ) {
            for ( int jj=0; jj<arrayPaired[ii].size(); jj++ ) {
                if (jj%2==0) cout << (int)arrayPaired[ii][jj] << " "; 
                if (jj%2==1) cout << (char)arrayPaired[ii][jj] << " "; 
            }
            cout << endl;
        }
    }
    //*/    ^^^FOR DEBUG^^^
    if (true) {
        for (int c=0; c<=sec_struct.size()-1; c++) { cout << (char)sec_struct[c]; };
        cout << endl << operate << endl;
    }
    //*** DEBUG fini ***//
}

template <class T>
void dump_vector (vector <T> myV) {
    for (int i=0; i<myV.size(); i++) {
        cout << (char)myV[i];
    }
    cout << endl;
}


#include "reportStruct.cpp"

int main(int argc, char *argv[]) {

    char* fastaName     = argv[1];
    char* FileNameIn    = argv[2];
    char* struct_descr  = argv[3];
    char* SPofWholeFas  = argv[4];
    char* seq_length    = argv[5];

    if (seq_length != NULL) {
        seq_size = from_string<int>(seq_length, std::dec);
//        cout << "GOT LENGTH = " << seq_size << endl;
    }

    baseOfGSP = from_string<int>(SPofWholeFas, std::dec);
    
    string sequence = "";
    fasta = fastaName;

    ifstream FileInput;

    FileInput.open (FileNameIn);
    if (!FileInput) {
        cout << "File: "    << FileNameIn
             << " failed to be open!!"  << endl;
        exit(1);
    }

    char * pch;
    pch = strtok(struct_descr, "%");
    if (pch == NULL) {
        cout << "NO DESCRIPTION IN FILE" << endl;
        exit(1);
    }
    string descr;
    while (pch != NULL ) {
    //    printf ("%s ", pch);
        descr = pch;
        v_descr.push_back(descr);
        pch = strtok(NULL, "%");
        global_depth++;
    }

    while (FileInput.good()) {
        FileInput >> sequence;
        if (FileInput.good() && DEBUG) cout << sequence << " for test" << endl;
    }

    FileInput.close();

    int start_type = from_string<int>((v_descr[0]).substr(0, 1), std::dec);
    if (start_type == 4 ) {
        searchStem_gt4(sequence, 1);
        return 0;
    } else if (start_type == 1 ) {
        searchStem_start(sequence, 1);
    }

    return 0;
}

void searchStem_start (const string seqStr, const int depth) {

//cout << v_descr[depth-1] << endl;       // debug

    // Remove the structure type identifier
    // throw description sequence to 'tokenSplit' process
    shape_inf mySI = tokenSplit (string(v_descr[depth-1]));

    string operate;
    vector <string> v_seq_left;
vector <string> v_seq_wait_combine;

// move forward
    reportStruct myRS;
    int outer_stem1 = 0;
    string after_stem = "";
    resultOfPair myRop;
    reportSEQ myRSeq;
// move forward

    struct_counter = 0;

    for (int slideW=mySI.swL; slideW<=mySI.swU; slideW+=mySI.swStep) {
// Every start point, get only 1000 structures
//startPstruct = 0;
        for (int startP=0; startP<=(int)seqStr.length()-slideW; startP+=mySI.spStep) {

startPstruct = 0;
// Every start point, get only 1000 structures
//startPstruct = 0;

            // SET THE START POINT FOR BACKTRACK
            globalStartPoint = startP + baseOfGSP;  // Global variable
            //v_seq_left.clear();
            vector<string>().swap(v_seq_left);

            operate = seqStr.substr(startP, slideW);
            complete_seq = operate;

            // Reporting section!!!     report_descr
            //reportStruct myRS = reportParse (mySI.report_descr);
            myRS = reportParse (mySI.report_descr);
            //int outer_stem1 = 0;
            outer_stem1 = 0;
            //string after_stem = "";
            after_stem = "";
            if ((myRS.oMinus != -1) && (startP+slideW < seqStr.size())) {
                // for riboswitches like PreQ1 and SAM_alpha,
                //  those have outer sequence after stem 1
                outer_stem1 = myRS.oMinus;
                after_stem += seqStr.substr(startP+slideW, myRS.oMinus);
//                after_stem = myRop.outer_h3.substr(0, myRS.oMinus);
                complete_seq += seqStr.substr(startP+slideW, outer_stem1);
//                for (int i=0; i<complete_seq.size()-myRop.sec_struct.size(); i++) {
//                    myRop.sec_struct.push_back('.');
//                }
            }

            //resultOfPair myRop = getRnP (operate, mySI, 0, "", after_stem);
            myRop = getRnP (operate, mySI, 0, "", after_stem);
//cout << "HI! my conti is " << myRop.conti << endl;
//report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);
            if (myRop.conti) { continue; }
//cout << "\033[32mxxxxxxxxxxxxxxxxxxx LoopLower: " << mySI.loopLower << "\t LoopUpper: " << mySI.loopUpper << "\033[m" << endl;
            if ( (myRop.loop_interval.size() > mySI.loopUpper) || (myRop.loop_interval.size() < mySI.loopLower) ) { continue; }


//*
            // Reporting section!!!     report_descr
            if ((myRS.oMinus != -1) && (startP+slideW < seqStr.size())) {
                // for riboswitches like PreQ1 and SAM_alpha,
                //  those have outer sequence after stem 1
                for (int i=0, full_size=myRop.sec_struct.size(); i<complete_seq.size()-full_size; i++) {
                    myRop.sec_struct.push_back('.');
                }
            }
//*/

            //reportSEQ myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);
            myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);

// FOR DEBUG
/*
report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);
report_checking (myRS);    // test report value and print the values
seq_checking (myRSeq);
dump_vector (myRop.sec_struct);
cout << complete_seq << "\ndepth: " << depth << "\tglobal: " << global_depth << endl;
*/

            // setting push No.
            if (mySI.pushNo == 1) {
                if (myRS.lPlus == -1 && myRS.lMinus == -1) {    // no need loop
                    myRop.sec_struct.push_back(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1);
                    v_seq_left.push_back(operate.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1, myRop.loop_interval.size()));
                    //myRop.arrayPaired[myRop.arrayPaired.size()-1][2] - myRop.arrayPaired[myRop.arrayPaired.size()-1][0] -1));
                } else if (myRS.lMinus == -1){                  // need loop plus
                    myRop.sec_struct.push_back(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] + myRS.lPlus +1);
                    v_seq_left.push_back(operate.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] + myRS.lPlus +1, myRop.arrayPaired[myRop.arrayPaired.size()-1][2] - myRop.arrayPaired[myRop.arrayPaired.size()-1][0] -1 -myRS.lPlus));
                } else if (myRS.lPlus == -1){                   // need loop minus
                    myRop.sec_struct.push_back(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1);
                    v_seq_left.push_back(operate.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1, myRop.arrayPaired[myRop.arrayPaired.size()-1][2] - myRop.arrayPaired[myRop.arrayPaired.size()-1][0] -1 -myRS.lMinus));
                } else {                                        // need loop plus and minus
                    myRop.sec_struct.push_back(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] + myRS.lPlus +1);
                    v_seq_left.push_back(operate.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] + myRS.lPlus +1, myRop.loop_interval.size() - myRS.lPlus - myRS.lMinus));
                    //myRop.arrayPaired[myRop.arrayPaired.size()-1][2] - myRop.arrayPaired[myRop.arrayPaired.size()-1][0] -1 - myRS.lPlus - myRS.lMinus));
                }
            }
//cout << "\033[33m>>>>>>>>>>>> startP: " << startP << "\tpush back in sec_struct: " << myRop.sec_struct.back() << "\033[m" <<endl;

//            vector <string> v_seq_wait_combine;

            // check next calling method
            if (depth < global_depth) {
                //int gtype;// = atoi(((v_descr[depth]).substr(0, 1)).c_str());
                int gtype = from_string<int>((v_descr[depth]).substr(0, 1), std::dec);
                //cout << "=================NEXT TYPE: " << gtype << endl;

                switch (gtype)
                {
                    case 2:
//cout << "\033[34mSEND:\n" << v_seq_left[0]  << "\ncomplete:\n" << seqStr <<  "\033[m" << endl;
                        searchStem_gt2(v_seq_left, v_seq_wait_combine, myRSeq.h5, myRSeq.h3, myRop.sec_struct, depth+1);
                        break;
                    case 3:
                        searchStem_gt3(v_seq_left, v_seq_wait_combine, myRSeq.h5, myRSeq.h3, myRop.sec_struct, depth+1);
                        break;
                }
            } else {    // report the final result
                struct_counter++;   // Found 1 more structrue
                startPstruct++;

                //cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                if (seq_size == 0) {
                    cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                } nclude "include.h"
#include <stdexcept>

const int DEBUG_struct = false;
const int DEBUG = false;
const int STRUCT_LIMIT = 100;

using namespace std;
using namespace pcrecpp;

// subfunction commence
// #include "c_def_struct.h"
// #include "prototype.h"
// #include "seq_pattern.cpp"
// #include "stempair.cpp"
// #include "getPairRatio.cpp"
// #include "tokenSplit.cpp"
// #include "getPnR.cpp"
// #include "seq_info.cpp"
// // subfunction fini
// /*
// inline void searchStem_start(const string seqStr, const int);
// inline void searchStem_gt2(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
// inline void searchStem_gt3(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
// inline void searchStem_gt4(const string seqStr, const int);
// */
// void searchStem_start(const string seqStr, const int);
// //void searchStem_gt2(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
// void searchStem_gt2(vector <string> const &, vector <string> const &, string const &, const string &, vector <int> const &, const int &);
// //void searchStem_gt3(const vector <string>, const vector <string>, const string, const string, vector <int>, const int);
// void searchStem_gt3(vector <string> const &, vector <string> const &, string const &, const string &, vector <int> const &, const int &);
// void searchStem_gt4(const string seqStr, const int);
//
// //const float pair_ratio   = 2/3;
//
// int globalStartPoint    = 0;
// int baseOfGSP           = 0;
// int struct_counter      = 0;
// int seq_size = 0;
//
//
// string complete_seq;
//
// string fasta;
//
// /*
//  * Nouveau variables commence ici
//   */
//
//   int global_depth = 0;   // Global depth, for result reporting
//   int startPstruct =0;    // control structures of each startP
//   vector <string> v_descr;
//
//   /*
//    * Nouveau variables fini ici
//    else {
                    cout << ">" << fasta << "=anti=" << seq_size - globalStartPoint << "s" << struct_counter << endl;
                }
                cout << myRSeq.h5 + myRSeq.loop + myRSeq.h3 << endl;
                cout << "# ";
                dump_vector (myRop.sec_struct);
                cout << "# " << complete_seq << endl;
            }
            if (mySI.pushNo == 1) {
                myRop.sec_struct.pop_back();
            }

            operate.clear();
        }
    }
}

//void searchStem_gt2(vector <string> const &v_seq_left_ORI, const vector <string> v_seq_wait_combine_ORI, const string h5loop_ORI, const string h3loop_ORI, const vector <int> sec_struct_ORI, const int depth) {
void searchStem_gt2(vector <string> const &v_seq_left_ORI, vector <string> const &v_seq_wait_combine_ORI, string const &h5loop_ORI, string const &h3loop_ORI, vector <int> const &sec_struct_ORI, int const &depth) {

//if (struct_counter >= 1000 ) { return; }
if (startPstruct >= STRUCT_LIMIT ) { return; }

    shape_inf mySI = tokenSplit (string(v_descr[depth-1]));
    
    string operate;

vector <string> v_seq_left = v_seq_left_ORI;
vector <string> v_seq_wait_combine;
vector <int> sec_struct;

// move forward
    int backOfSS = 0;
    string seq_left = "";
    string seq_before_op = "";
    string seq_after_op = "";
    resultOfPair myRop;
    reportStruct myRS;
    reportSEQ myRSeq;
// move forward

    for (int slideW=mySI.swU; slideW>=mySI.swL; slideW--) {
        for (int startP=mySI.wts_min; startP<=mySI.wts_max; startP++) {

//            vector <string> v_seq_left = v_seq_left_ORI;
//            vector <string> v_seq_wait_combine = v_seq_wait_combine_ORI;

            // clear vector contain
            vector<string>().swap(v_seq_left);
            vector<string>().swap(v_seq_wait_combine);
            vector<int>().swap(sec_struct);

            v_seq_left = v_seq_left_ORI;
            v_seq_wait_combine = v_seq_wait_combine_ORI;

            //string seq_left = v_seq_left.back();
            seq_left = v_seq_left.back();
//cout << v_seq_left.back() << " " ;
            v_seq_left.pop_back();
//            vector <int> sec_struct = sec_struct_ORI;
            sec_struct = sec_struct_ORI;

            //int backOfSS = sec_struct.back();
            backOfSS = sec_struct.back();
//cout << backOfSS << endl;
            if (seq_left.size() < startP) { continue; }
            operate = seq_left.substr(startP, slideW);
            if (slideW > seq_left.size()) { continue; }
//cout << v_seq_left.back() << " " << endl;
            sec_struct.pop_back();

            //string seq_before_op = seq_left.substr(0, startP);
            //string seq_after_op = seq_left.substr(startP + operate.size());
            seq_before_op = seq_left.substr(0, startP);
            seq_after_op = seq_left.substr(startP + operate.size());

            //resultOfPair myRop = getRnP (operate, mySI, startP, seq_before_op, seq_after_op);
            myRop = getRnP (operate, mySI, startP, seq_before_op, seq_after_op);
            if (myRop.conti) { continue; }
            if (myRop.arrayPaired[0][0]!=0) { continue; }
            if ( (myRop.loop_interval.size() > mySI.loopUpper) || (myRop.loop_interval.size() < mySI.loopLower) ) { continue; }

            // Reporting section!!!     report_descr
            //reportStruct myRS = reportParse (mySI.report_descr);
            //reportSEQ myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);
            myRS = reportParse (mySI.report_descr);
            myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);
            myRSeq.h5 = h5loop_ORI + myRSeq.h5;
            //myRSeq.h3 = h3loop_ORI;
            // 
            if ((depth >= 3) && (mySI.pushNo == 0) &&
                (from_string<int>((v_descr[depth-2]).substr(0, 1), std::dec)==3 ||
                 from_string<int>((v_descr[depth-3]).substr(0, 1), std::dec)==3 ) && 
                 !v_seq_wait_combine.empty()                                    )
            {
                myRSeq.loop = myRSeq.loop + myRSeq.h3;
                myRSeq.loop = myRSeq.loop + v_seq_wait_combine.back();
                v_seq_wait_combine.pop_back();
                while ((depth == global_depth) && (!v_seq_wait_combine.empty())) {
                    myRSeq.loop = myRSeq.loop + v_seq_wait_combine.back();
                    v_seq_wait_combine.pop_back();
                }
            } else {
                myRSeq.loop = myRSeq.loop + myRSeq.h3;
            }

            // add struct
            for (int i=0; i<myRop.sec_struct.size(); i++) {
                sec_struct[backOfSS + startP +i] = myRop.sec_struct[i];
            }

// FOR DEBUG
/*
report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);
report_checking (myRS);    // test report value and print the values
seq_checking (myRSeq);
dump_vector (sec_struct);
cout << complete_seq << "\ndepth: " << depth << "\tglobal: " << global_depth << endl;
*/


            // setting push No.
            if (mySI.pushNo == 1) {
                if(myRS.oMinus == -1) {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] +1 +startP));
                } else {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP));
                }
            }
//cout << v_seq_left.back() << " " << sec_struct[sec_struct.size()-1] << endl;

            // check next calling method
            if (depth < global_depth) {
                //int gtype;// = atoi(((v_descr[depth]).substr(0, 1)).c_str());
                int gtype = from_string<int>((v_descr[depth]).substr(0, 1), std::dec);
                //cout << "-NEXT TYPE: " << gtype << endl;

                switch (gtype)
                {
                    case 2:
//for(int i=0; i<v_seq_left.size(); i++){ cout << v_seq_left[i] << endl; }
                        searchStem_gt2(v_seq_left, v_seq_wait_combine, myRSeq.h5 + myRSeq.loop, h3loop_ORI, sec_struct, depth+1);
                        break;
                    case 3:
                        //searchStem_gt3(seqStr, h5loop, h3loop, sec_structure, depth+1);
//if(depth==5) cout << v_seq_left_ORI.back() << endl;
                        searchStem_gt3(v_seq_left, v_seq_wait_combine, myRSeq.h5 + myRSeq.loop, h3loop_ORI, sec_struct, depth+1);
                        break;
                }
            } else {    // report the final result
//if(depth==5) cout << "SEQ LEFT: " << v_seq_left_ORI.back() << endl;
                struct_counter++;   // Found 1 more structrue
                startPstruct++;

                //cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                if (seq_size == 0) {
                    cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                } else {
                    cout << ">" << fasta << "=anti=" << seq_size - globalStartPoint << "s" << struct_counter << endl;
                }
                cout << myRSeq.h5 + myRSeq.loop + h3loop_ORI << endl;
                cout << "# ";
                dump_vector (sec_struct);
                cout << "# " << complete_seq << endl;
            }
            if (mySI.pushNo == 1) {
                sec_struct.pop_back();
            }

            operate.clear();
        }
    }
}



//void searchStem_gt3(const vector <string> v_seq_left_ORI, const vector <string> v_seq_wait_combine_ORI, const string h5loop_ORI, const string h3loop_ORI, vector <int> sec_struct_ORI, const int depth) {
void searchStem_gt3(vector <string> const &v_seq_left_ORI, vector <string> const &v_seq_wait_combine_ORI, string const &h5loop_ORI, string const &h3loop_ORI, vector <int> const &sec_struct_ORI, int const &depth) {

//if (struct_counter >= 1000 ) { return; }
if (startPstruct >= STRUCT_LIMIT ) { return; }

    shape_inf mySI = tokenSplit (string(v_descr[depth-1]));
    
    string operate;

vector <string> v_seq_left;
vector <string> v_seq_wait_combine;
vector <int> sec_struct;

// move forward
    string seq_left = "";
    int backOfSS = 0;
    string seq_before_op = "";
    string seq_after_op = "";
    reportStruct myRS;
    reportSEQ myRSeq;
// move forward

    for (int slideW=mySI.swU; slideW>=mySI.swL; slideW--) {
        for (int startP=mySI.wts_min; startP<=mySI.wts_max; startP++) {

//if (depth==6) cout << v_seq_left_ORI.back() << endl << v_seq_left_ORI[0] << endl << "\t" << endl;
//            vector <string> v_seq_left = v_seq_left_ORI;
//            vector <string> v_seq_wait_combine = v_seq_wait_combine_ORI;

            // clear vector contain
            vector<string>().swap(v_seq_left);
            vector<string>().swap(v_seq_wait_combine);
            vector<int>().swap(sec_struct);

            v_seq_left = v_seq_left_ORI;
            v_seq_wait_combine = v_seq_wait_combine_ORI;

            //string seq_left = v_seq_left.back();
            seq_left = v_seq_left.back();
            v_seq_left.pop_back();
//            vector <int> sec_struct = sec_struct_ORI;
            sec_struct = sec_struct_ORI;

            //int backOfSS = sec_struct.back();
            backOfSS = sec_struct.back();
            if (seq_left.size() < startP) { continue; }
            operate = seq_left.substr(startP, slideW);
            if (slideW > seq_left.size()) { continue; }
            sec_struct.pop_back();

            //string seq_before_op = seq_left.substr(0, startP);
            //string seq_after_op = seq_left.substr(startP + operate.size());
            seq_before_op = seq_left.substr(0, startP);
            seq_after_op = seq_left.substr(startP + operate.size());

            resultOfPair myRop = getRnP (operate, mySI, startP, seq_before_op, seq_after_op);
            if (myRop.conti) { continue; }
            if (myRop.arrayPaired[0][0]!=0) { continue; }
            if ( (myRop.loop_interval.size() > mySI.loopUpper) || (myRop.loop_interval.size() < mySI.loopLower) ) { continue; }

            // Reporting section!!!     report_descr
            //reportStruct myRS = reportParse (mySI.report_descr);
            //reportSEQ myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);
            myRS = reportParse (mySI.report_descr);
            myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);
            myRSeq.h5 = h5loop_ORI + myRSeq.h5;
            //myRSeq.h3 = h3loop_ORI;
            // move sequence wait to combine into vector
            v_seq_wait_combine.push_back(myRSeq.h3);


            // add struct
            for (int i=0; i<myRop.sec_struct.size(); i++) {
                sec_struct[backOfSS + startP +i] = myRop.sec_struct[i];
            }

// FOR DEBUG
/*
report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);
report_checking (myRS);    // test report value and print the values
seq_checking (myRSeq);
dump_vector (sec_struct);
cout << complete_seq << endl;
*/



            // setting push No.
            if (mySI.pushNo == 2) {
                // sequence after the branch
                if(myRS.oMinus == -1) {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] +1 +startP));
                } else {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP));
                }
                // sequence inside the branch
                sec_struct.push_back(backOfSS + myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1 +startP);
                v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1 +startP, myRop.loop_interval.size()));
            }

            // setting push No.
            if (mySI.pushNo == 1) {
/*                if(myRS.oMinus == -1) {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] +1 +startP));
                } else {
                    sec_struct.push_back(backOfSS + myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP);
                    v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[0][2] + myRS.oMinus +1 +startP));
                }
*/
                sec_struct.push_back(backOfSS + myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1 +startP);
                v_seq_left.push_back(seq_left.substr(myRop.arrayPaired[myRop.arrayPaired.size()-1][0] +1 +startP, myRop.loop_interval.size()));
            }

/*cout << v_seq_left[0] << " " << sec_struct[sec_struct.size()-2] << endl;
cout << v_seq_left[1] << " " << sec_struct[sec_struct.size()-1] << endl;
cout << v_seq_left.back() << " " << sec_struct[sec_struct.size()-1] << endl;*/

            // check next calling method
            if (depth < global_depth) {
                //int gtype;// = atoi(((v_descr[depth]).substr(0, 1)).c_str());
                int gtype = from_string<int>((v_descr[depth]).substr(0, 1), std::dec);
                //cout << "-NEXT TYPE: " << gtype << endl;

                switch (gtype)
                {
                    case 2:
                        searchStem_gt2(v_seq_left, v_seq_wait_combine, myRSeq.h5 + myRSeq.loop, h3loop_ORI, sec_struct, depth+1);
                        break;
                    case 3:
                        //searchStem_gt3(seqStr, h5loop, h3loop, sec_structure, depth+1);
                        searchStem_gt3(v_seq_left, v_seq_wait_combine, myRSeq.h5 + myRSeq.loop, h3loop_ORI, sec_struct, depth+1);
                        break;
                }
            } else {    // report the final result
                struct_counter++;   // Found 1 more structrue
                startPstruct++;

                //cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                if (seq_size == 0) {
                    cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                } else {
                    cout << ">" << fasta << "=anti=" << seq_size - globalStartPoint << "s" << struct_counter << endl;
                }
                cout << myRSeq.h5 + myRSeq.loop + h3loop_ORI << endl;
                cout << "# ";
                dump_vector (sec_struct);
                cout << "# " << complete_seq << endl;
            }
            if (mySI.pushNo == 1) {
                sec_struct.pop_back();
            }

            operate.clear();
        }
    }
}


void searchStem_gt4 (const string seqStr, const int depth) {

    // Remove the structure type identifier
    // throw description sequence to 'tokenSplit' process
    shape_inf mySI = tokenSplit (string(v_descr[depth-1]));

    string operate;
    vector <string> v_seq_left;
vector <string> v_seq_wait_combine;

// move forward
    int outer_stem1 = 0;
    string after_stem = "";
    reportStruct myRS;
    resultOfPair myRop;
    vector <int> sec_struct;
// move forward

    struct_counter = 0;

    for (int slideW=mySI.swL; slideW<=mySI.swU; slideW+=mySI.swStep) {
// Every start point, get only 1000 structures
//startPstruct = 0;
        for (int startP=0; startP<=(int)seqStr.length()-slideW; startP+=mySI.spStep) {

startPstruct = 0;
// Every start point, get only 1000 structures
//startPstruct = 0;

            // SET THE START POINT FOR BACKTRACK
            globalStartPoint = startP + baseOfGSP;  // Global variable
            //v_seq_left.clear();
            vector<string>().swap(v_seq_left);

            operate = seqStr.substr(startP, slideW);
            complete_seq = operate;

            // Reporting section!!!     report_descr
            //reportStruct myRS = reportParse (mySI.report_descr);
            myRS = reportParse (mySI.report_descr);
            //int outer_stem1 = 0;
            //string after_stem = "";
            outer_stem1 = 0;
            after_stem = "";
/*
            if ((myRS.oMinus != -1) && (startP+slideW < seqStr.size())) {
                // for riboswitches like PreQ1 and SAM_alpha,
                //  those have outer sequence after stem 1
                outer_stem1 = myRS.oMinus;
                after_stem += seqStr.substr(startP+slideW, myRS.oMinus);
                complete_seq += seqStr.substr(startP+slideW, outer_stem1);
            }
*/

            //resultOfPair myRop;
            //vector <int> sec_struct(operate.size(), '.');
            sec_struct.assign(operate.size(), '.');
            myRop.sec_struct = sec_struct;  // Recording 2nd structure
            //sec_struct.~vector();
            myRop.h5_pairing_seq = "";
            myRop.h3_pairing_seq = "";
            myRop.loop_interval = operate;
            myRop.outer_h5 = "";
            myRop.outer_h3 = "";
//report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);

            // Checking if sequence pattern is corresponding provided restrict
            if ( (mySI.seq_info_descr != "") && !(report_seq_info(myRop, mySI.seq_info_descr)) ) { //cout << myRop.loop_interval << endl << mySI.seq_info_descr << endl;
                continue; }


/*
            // Reporting section!!!     report_descr
            if ((myRS.oMinus != -1) && (startP+slideW < seqStr.size())) {
                // for riboswitches like PreQ1 and SAM_alpha,
                //  those have outer sequence after stem 1
                for (int i=0, full_size=myRop.sec_struct.size(); i<complete_seq.size()-full_size; i++) {
                    myRop.sec_struct.push_back('.');
                }
            }
*/

//            reportSEQ myRSeq = reporter (myRS, myRop.h5_pairing_seq, myRop.h3_pairing_seq, myRop.loop_interval, myRop.outer_h5, myRop.outer_h3);

// FOR DEBUG
/*
report_pairing_result4debug (operate, myRop.arrayPaired, myRop.sec_struct);
report_checking (myRS);    // test report value and print the values
seq_checking (myRSeq);
dump_vector (myRop.sec_struct);
cout << complete_seq << "\ndepth: " << depth << "\tglobal: " << global_depth << endl;
*/

            // setting push No.
            if (mySI.pushNo == 1) {
                if (myRS.lPlus == -1 && myRS.lMinus == -1) {    // no need loop
                    myRop.sec_struct.push_back(0);
                    v_seq_left.push_back(operate);
                } else if (myRS.lMinus == -1){                  // need loop plus
                    myRop.sec_struct.push_back(myRS.lPlus);
                    v_seq_left.push_back(operate.substr(myRS.lPlus, operate.size() -myRS.lPlus));
                } else if (myRS.lPlus == -1){                   // need loop minus
                    myRop.sec_struct.push_back(0);
                    v_seq_left.push_back(operate.substr(0, operate.size() -myRS.lMinus));
                } else {                                        // need loop plus and minus
                    myRop.sec_struct.push_back(myRS.lPlus);
                    v_seq_left.push_back(operate.substr(myRS.lPlus, operate.size() - myRS.lPlus - myRS.lMinus));
                }
            }

            // check next calling method
            if (depth < global_depth) {
                int gtype = from_string<int>((v_descr[depth]).substr(0, 1), std::dec);

                switch (gtype)
                {
                    case 2:
                        //searchStem_gt2(v_seq_left, v_seq_wait_combine, myRSeq.h5, myRSeq.h3, myRop.sec_struct, depth+1);
                        searchStem_gt2(v_seq_left, v_seq_wait_combine, "", "", myRop.sec_struct, depth+1);
                        break;
                    case 3:
                        //searchStem_gt3(v_seq_left, v_seq_wait_combine, myRSeq.h5, myRSeq.h3, myRop.sec_struct, depth+1);
                        searchStem_gt3(v_seq_left, v_seq_wait_combine, "", "", myRop.sec_struct, depth+1);
                        break;
                }
            } else {    // report the final result
                struct_counter++;   // Found 1 more structrue
                startPstruct++;

                if (seq_size == 0) {
                    cout << ">" << fasta << "=" << globalStartPoint << "s" << struct_counter << endl;
                } else {
                    cout << ">" << fasta << "=anti=" << seq_size - globalStartPoint << "s" << struct_counter << endl;
                }
                //cout << myRSeq.h5 + myRSeq.loop + myRSeq.h3 << endl;
                cout << "# ";
                dump_vector (myRop.sec_struct);
                cout << "# " << complete_seq << endl;
            }
            if (mySI.pushNo == 1) {
                myRop.sec_struct.pop_back();
            }

            operate.clear();
        }
    }
}

