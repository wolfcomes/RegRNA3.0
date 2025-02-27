using namespace pcrecpp;

inline bool pattern_match (string pattern, string src) {
    if(RE(pattern).PartialMatch(src)) {
        return true;
    }
    else {
        return false;
    }
}
