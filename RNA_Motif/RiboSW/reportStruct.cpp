using namespace std;

reportStruct reportParse (string report_descr) {
    reportStruct myRS;
    memset(&myRS, -1, sizeof(myRS));

    vector <string> v_reportValues;
    char * pch;
    char * cstr4strtok;
    cstr4strtok = strdup(report_descr.c_str());
//cout << cstr4strtok << endl;
    pch = strtok(cstr4strtok, ",");
    while (pch != NULL) {
        v_reportValues.push_back(pch);
        pch = strtok(NULL, ",");
    }
    
    while (v_reportValues.size() != 0) {
        string pv = v_reportValues.back();

        if (pv.find("l+")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT loop plus, value= " << pv.substr(pv.find("l+")+2) << endl; }
            myRS.lPlus  = from_string<int>(pv.substr(pv.find("l+")+2), std::dec);
        } else if (pv.find("l-")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT loop minus, value= " << pv.substr(pv.find("l-")+2) << endl; }
            myRS.lMinus = from_string<int>(pv.substr(pv.find("l-")+2), std::dec);
        } else if (pv.find("l=")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT loop full, value= " << pv.substr(pv.find("l=")+2) << endl; }
            myRS.lFull  = from_string<int>(pv.substr(pv.find("l=")+2), std::dec);
        } else if (pv.find("h+")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT helix plus, value= " << pv.substr(pv.find("h+")+2) << endl; }
            myRS.hPlus  = from_string<int>(pv.substr(pv.find("h+")+2), std::dec);
        } else if (pv.find("h-")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT helix minus, value= " << pv.substr(pv.find("h-")+2) << endl; }
            myRS.hMinus = from_string<int>(pv.substr(pv.find("h-")+2), std::dec);
        } else if (pv.find("h=")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT helix full, value= " << pv.substr(pv.find("h=")+2) << endl; }
            myRS.hFull  = from_string<int>(pv.substr(pv.find("h=")+2), std::dec);
        } else if (pv.find("o+")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT outer plus, value= " << pv.substr(pv.find("o+")+2) << endl; }
            myRS.oPlus  = from_string<int>(pv.substr(pv.find("o+")+2), std::dec);
        } else if (pv.find("o-")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT outer minus, value= " << pv.substr(pv.find("o-")+2) << endl; }
            myRS.oMinus = from_string<int>(pv.substr(pv.find("o-")+2), std::dec);
        } else if (pv.find("o=")!=string::npos) {
            if(DEBUG_struct) { cout << "GOT outer full, value= " << pv.substr(pv.find("o=")+2) << endl; }
            myRS.oFull  = from_string<int>(pv.substr(pv.find("o=")+2), std::dec);
        }

        v_reportValues.pop_back();
    }

    free(cstr4strtok);

    return myRS;
}


reportSEQ reporter (const reportStruct myRS, const string h5ps, const string h3ps, const string loopi, const string oh5, const string oh3) {
    reportSEQ myReportSeq;
    //memset(&myReportSeq, 0, sizeof(myReportSeq));
    //*
    myReportSeq.h5 = "";
    myReportSeq.h3 = "";
    myReportSeq.loop = "";
    //*/
    
    string tempStr = "";
    
    if (myRS.oPlus != -1) {
        string oh5_copy = oh5;
        reverse(oh5_copy.begin(), oh5_copy.end());
//cout << "outer h5:\t" << oh5 << endl;
        myReportSeq.h5 += oh5_copy.substr(0, myRS.oPlus);
        reverse(myReportSeq.h5.begin(), myReportSeq.h5.end());
        if (myRS.oPlus == 0) {
            myReportSeq.h5 += oh5;
        }
//cout << "outer my h5:\t" << myReportSeq.h5 << endl;
    } if (myRS.oMinus != -1) {
        myReportSeq.h3 += oh3.substr(0, myRS.oMinus);
        if (myRS.oMinus == 0) {
            myReportSeq.h3 += oh3;
        }
    }
    if (myRS.hFull != -1) {     // need helix full sequence
        tempStr = h3ps;
        myReportSeq.h5 += h5ps;
        reverse(tempStr.begin(), tempStr.end());
        myReportSeq.h3 = tempStr + myReportSeq.h3;
    } else {
        if (myRS.hPlus != -1) {     // need 5' terminal seq. of helix
            myReportSeq.h5 += h5ps.substr(0, myRS.hPlus);
            tempStr = h3ps.substr(0, myRS.hPlus);
            reverse(tempStr.begin(), tempStr.end());
            myReportSeq.h3 = tempStr + myReportSeq.h3;
        }
        if (myRS.hMinus != -1 ) {   // need 3' terminal seq. of helix
            //myReportSeq.h5 += h5ps.substr(h5ps.size()-myRS.hMinus, myRS.hMinus);      // maybe error, cuz start pos may negative
            string tempHxps = h5ps;
            reverse(tempHxps.begin(), tempHxps.end());
            tempStr = tempHxps.substr(0, myRS.hMinus);
            reverse(tempStr.begin(), tempStr.end());
            myReportSeq.h5 += tempStr;
            //tempStr = h3ps.substr(h3ps.size()-myRS.hMinus, myRS.hMinus);
            tempHxps = h3ps;
            reverse(tempHxps.begin(), tempHxps.end());
            tempStr = tempHxps.substr(0, myRS.hMinus);
            myReportSeq.h3 = tempStr + myReportSeq.h3;
        }
    }
    if (myRS.lFull != -1) {     // need loop full sequence
        myReportSeq.loop = loopi;
    } else {
        if (myRS.lPlus != -1) {     // need 5' terminal seq. of loop region
            myReportSeq.h5 += loopi.substr(0, myRS.lPlus);
        }
        if (myRS.lMinus != -1) {    // need 3' terminal seq. of loop region
            //myReportSeq.h3 = loopi.substr(loopi.size()-myRS.lMinus, myRS.lMinus);
            string loopi_copy = loopi;
            reverse(loopi_copy.begin(), loopi_copy.end());
            tempStr = loopi_copy.substr(0, myRS.lMinus);
            reverse(tempStr.begin(), tempStr.end());
            myReportSeq.h3 = tempStr + myReportSeq.h3;
        }
    }
    //cout << "h5: '" << myReportSeq.h5 << "'\th3: '" << myReportSeq.h3 << "'\tloop: '" << myReportSeq.loop << "'\t after report" << endl;

    return myReportSeq;
}
