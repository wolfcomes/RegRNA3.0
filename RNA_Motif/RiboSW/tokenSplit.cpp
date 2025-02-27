shape_inf tokenSplit (const string v_descr_ORI) {
    shape_inf mySI;

    // Remove the structure type identifier
    string v_descr = v_descr_ORI.substr(v_descr_ORI.find(";")+1);
    
    string shape_descr;
    string report_descr;
    string seq_info_descr = "";
    
    // Token split, shape and report section
    char * pch;
    char * cstr4strtok;
    cstr4strtok = strdup(v_descr.c_str());
    pch = strtok(cstr4strtok, "'");
    if (pch == NULL) { cout << "No description!!" << endl; } 
    shape_descr = pch;
    pch = strtok(NULL, "'");
    report_descr = pch;
    pch = strtok(NULL, "'");
    if (pch != NULL) {
        seq_info_descr = pch;
    }

    // Token split, get shape information
    /*
        0: swL
        1: swU
        2: swStep
        3: spStep
        4: stemLower
        5: stemUpper
        6: loopLower
        7: loopUpper
        8: bulge
        9: pushNo
        [optional]
            10: where to start (min) - wts_min
            11: where to start (max) - wts_max
            12: where to end (min) - wte_min
            13: where to end (max) - wte_max
     */
    vector <int> v_shape_descr;
    cstr4strtok = strdup(shape_descr.c_str());
    pch = strtok(cstr4strtok, ",");
    while (pch != NULL) { 
        v_shape_descr.push_back(atoi(pch));
        pch = strtok(NULL, ",");
    }
    free(cstr4strtok);

    // for DEBUG
    //for ( int i=0; i<v_shape_descr.size(); i++) { cout << v_shape_descr[i] << " "; } cout << endl;

    mySI.swL            = v_shape_descr[0];
    mySI.swU            = v_shape_descr[1];
    mySI.swStep         = v_shape_descr[2];
    mySI.spStep         = v_shape_descr[3];
    mySI.stemLower      = v_shape_descr[4];
    mySI.stemUpper      = v_shape_descr[5];
    mySI.loopLower      = v_shape_descr[6];
    mySI.loopUpper      = v_shape_descr[7];
    mySI.bulge          = v_shape_descr[8];
    mySI.pushNo         = v_shape_descr[9];
    mySI.report_descr   = report_descr;
    mySI.seq_info_descr = seq_info_descr;
    if (v_shape_descr.size()>=12) {
        mySI.wts_min    = v_shape_descr[10];
        mySI.wts_max    = v_shape_descr[11];
    } else {
        mySI.wts_min    = -1;
        mySI.wts_max    = -1;
    }
    if (v_shape_descr.size()==14) {
        mySI.wte_min    = v_shape_descr[12];
        mySI.wte_max    = v_shape_descr[13];
    } else {
        mySI.wte_min    = -1;
        mySI.wte_max    = -1;
    }

    return mySI;
}
