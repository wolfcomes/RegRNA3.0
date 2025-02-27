bool report_seq_info (resultOfPair myRop, string seq_info_descr) {

    vector <string> v_seq_info_descr;
    char * pch;
    char * cstr4strtok;
    cstr4strtok = strdup(seq_info_descr.c_str());
    pch = strtok(cstr4strtok, ",");
    if (pch == NULL) { return true; }

    while (pch != NULL) {
        v_seq_info_descr.push_back(string(pch));
        pch = strtok(NULL, ",");
    }

    while (v_seq_info_descr.size() != 0) {
        string si = v_seq_info_descr.back();

        if (si.find("l+")!=string::npos) {
            if ( ! pattern_match(si.substr(si.find("l+")+2), myRop.loop_interval) ) { return false; }
        } else if (si.find("l-")!=string::npos) {
            if ( ! pattern_match(si.substr(si.find("l-")+2), myRop.loop_interval) ) { return false; }
        } else if (si.find("l=")!=string::npos) {
        } else if (si.find("h+")!=string::npos) {
            if ( ! pattern_match(si.substr(si.find("h+")+2), myRop.h5_pairing_seq) ) { return false; }
        } else if (si.find("h-")!=string::npos) {
            if ( ! pattern_match(si.substr(si.find("h-")+2), myRop.h3_pairing_seq) ) { return false; }
        } else if (si.find("h=")!=string::npos) {
        } else if (si.find("o+")!=string::npos) {
            if ( ! pattern_match(si.substr(si.find("o+")+2), myRop.outer_h5) ) { return false; }
        } else if (si.find("o-")!=string::npos) {
//cout << si.substr(si.find("o-")+2) << endl << myRop.outer_h3 << endl;
            if ( ! pattern_match(si.substr(si.find("o-")+2), myRop.outer_h3) ) { return false; }
        } else if (si.find("o=")!=string::npos) {
        }

//cout << v_seq_info_descr.back() << endl;
        v_seq_info_descr.pop_back();
    }

    free(cstr4strtok);

    return true;
}
