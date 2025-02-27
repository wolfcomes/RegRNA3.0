#include <vector>
//#include "include.h"

using namespace std;

inline bool match (char nuc1, char nuc2);

vector< vector <int> > searchStemPair (string h5, string h3, int stem_upper, int stem_lower, int op_length, const int BULGE) {
    /*
    cout << "\033[33mStart to search stem pairing!!" << endl
         << "Stem H5 is " << h5 << endl
         << "Stem H3 is " << h3 << "\033[m" << endl;
    */

    int h5_nuc, h3_nuc, missH5, missH3, pairs, bulge;
    h5_nuc = h3_nuc = missH5 = missH3 = pairs = bulge = 0;
    int h5p, h3p;
    h5p = h3p = -1;

    vector< vector <int> > arrayPaired;

    while ( (h5_nuc < stem_upper) && (h3_nuc < stem_upper) && (h5_nuc + h3_nuc < op_length) ) {
        ////cout << stem_upper << "\t" << h5_nuc << "\t" << h3_nuc << endl;
        if(match ((char)h5[h5_nuc], (char)h3[h3_nuc])) {
            h5p = h5_nuc;
            h3p = h3_nuc;
            arrayPaired.push_back(vector <int>());
            arrayPaired[pairs].push_back((int)h5p);
            arrayPaired[pairs].push_back((int)h5[h5p]);
            arrayPaired[pairs].push_back((int)(op_length-h3p-1));
            arrayPaired[pairs].push_back((int)h3[h3p]);
            /*
            cout << "\033[35mGOT PAIRING\033[m" << endl
                << h5[h5_nuc] <<" and " << h3[h3_nuc] << endl;
            cout << "pairs = " << pairs+1 << "\th5p = " << h5p << "\th3p = " << h3p << endl;
            cout << "position " << h5p << " [" << h5[h5p] << "] and \tposition " << -h3p << " [" << h3[h3p] << "] pairs" << endl;
            */
            h5_nuc++; h3_nuc++; pairs++; missH5=0; missH3=0; bulge=0;
        } else if (missH3 < BULGE && h3_nuc<stem_upper-1) {
            h3_nuc++;
            missH3++;
        } else if (missH5 < BULGE) {
            h5_nuc++;
            missH5++;
            bulge++;
            h3_nuc = h3p + bulge;
            missH3 = 0;
        } else {
            break;
        }
    }

    vector< vector <int> > arrayPaired2;
    h5_nuc = h3_nuc = missH5 = missH3 = pairs = bulge = 0;
    h5p = h3p = -1;

    while ( (h5_nuc < stem_upper) && (h3_nuc < stem_upper) && (h5_nuc + h3_nuc < op_length) ) {
        ////cout << stem_upper << "\t" << h5_nuc << "\t" << h3_nuc << endl;
        if(match ((char)h5[h5_nuc], (char)h3[h3_nuc])) {
            h5p = h5_nuc;
            h3p = h3_nuc;
            arrayPaired2.push_back(vector <int>());
            arrayPaired2[pairs].push_back((int)h5p);
            arrayPaired2[pairs].push_back((int)h5[h5p]);
            arrayPaired2[pairs].push_back((int)(op_length-h3p-1));
            arrayPaired2[pairs].push_back((int)h3[h3p]);
            /*
            cout << "\033[35mGOT PAIRING\033[m" << endl
                << h5[h5_nuc] <<" and " << h3[h3_nuc] << endl;
            cout << "pairs = " << pairs+1 << "\th5p = " << h5p << "\th3p = " << h3p << endl;
            cout << "position " << h5p << " [" << h5[h5p] << "] and \tposition " << -h3p << " [" << h3[h3p] << "] pairs" << endl;
            */
            h5_nuc++; h3_nuc++; pairs++; missH5=0; missH3=0; bulge=0;
        } else if (missH3 < BULGE && h3_nuc<stem_upper-1 && missH3<=missH5) {
            h3_nuc++;
            missH3++;
        } else if (missH5 < BULGE) {
            h5_nuc++;
            missH5++;
            bulge++;
            h3_nuc = h3p + bulge;
            missH3 = 0;
        } else {
            break;
        }
    }


    // FOLLOW INTERNAL LOOP
    vector< vector <int> > arrayPaired3;
    h5_nuc = h3_nuc = missH5 = missH3 = pairs = bulge = 0;
    h5p = h3p = -1;

    while ( (h5_nuc < stem_upper) && (h3_nuc < stem_upper) && (h5_nuc + h3_nuc < op_length) ) {
        ////cout << stem_upper << "\t" << h5_nuc << "\t" << h3_nuc << endl;
        if(match ((char)h5[h5_nuc], (char)h3[h3_nuc])) {
            h5p = h5_nuc;
            h3p = h3_nuc;
            arrayPaired3.push_back(vector <int>());
            arrayPaired3[pairs].push_back((int)h5p);
            arrayPaired3[pairs].push_back((int)h5[h5p]);
            arrayPaired3[pairs].push_back((int)(op_length-h3p-1));
            arrayPaired3[pairs].push_back((int)h3[h3p]);
            /*
            cout << "\033[35mGOT PAIRING\033[m" << endl
                << h5[h5_nuc] <<" and " << h3[h3_nuc] << endl;
            cout << "pairs = " << pairs+1 << "\th5p = " << h5p << "\th3p = " << h3p << endl;
            cout << "position " << h5p << " [" << h5[h5p] << "] and \tposition " << -h3p << " [" << h3[h3p] << "] pairs" << endl;
            */
            h5_nuc++; h3_nuc++; pairs++; missH5=0; missH3=0; bulge=0;
        } else if (missH3 < BULGE && h3_nuc<stem_upper-1 && missH3<=missH5) {
            h3_nuc++;
            missH3++;
            h5_nuc++;
            missH5++;
        } else if (missH5 < BULGE) {
            h5_nuc++;
            missH5++;
            //bulge++;
            //h3_nuc = h3p + bulge;
            //missH3 = 0;
        } else {
            break;
        }
    }

    /*/ report result
    if ( pairs >= stem_lower )
    for ( int ii=0; ii<arrayPaired.size(); ii++ ) {
        for ( int jj=0; jj<arrayPaired[ii].size(); jj++ ) {
            if (jj%2==0) cout << (int)arrayPaired[ii][jj] << " ";
            if (jj%2==1) cout << arrayPaired[ii][jj] << " ";
        }
        cout << endl;
    }*/
    // records how many nucleotides paired
    //arrayPaired.push_back(vector <char>(1, (char)pairs));
    ////cout << "Total " << (int)arrayPaired[arrayPaired.size()-1][0] << " nucleotides paired." << endl;

    if ( arrayPaired.size() > arrayPaired2.size() && arrayPaired2.size() > arrayPaired3.size()) {
        return arrayPaired;
    } else if (arrayPaired2.size() > arrayPaired3.size()){
        return arrayPaired2;
    } else {
        return arrayPaired3;
    }
}

inline bool match (char nuc1, char nuc2) {
    if
    (
        (nuc1=='a' && nuc2=='u') ||
        (nuc1=='u' && nuc2=='a') ||
        (nuc1=='c' && nuc2=='g') ||
        (nuc1=='g' && nuc2=='c') ||
        (nuc1=='g' && nuc2=='u') ||
        (nuc1=='u' && nuc2=='g')
    )   // check if (i,j) is pairing
    {
        return true;
    }
    else {
        return false;
    }
}
