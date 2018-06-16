#include <iostream>
#include <limits>
#include <vector>

#include "global_alignment.h"


using namespace std;

vector<pair<string,string>> results;

bool choose_maximum_size(const std::pair<string, string> &lhs, const std::pair<string,string> &rhs) {
    return lhs.first.length() < rhs.first.length() ;
}

int main()
{
    vector<string> dna_s = {"ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT", "ACTGACC"};
    int maximum = numeric_limits<int>::max() * -1;
    string center;
    int position_center;
    for(int i=0; i<dna_s.size(); i++){
        string actual_center = dna_s[i];
        int acum_score = 0;
        for(int j=0; j<dna_s.size(); j++){
            if(i!=j) acum_score += global_neddleman::get_maximum_score(actual_center, dna_s[j]);
        }
        if(acum_score > maximum){
            maximum = acum_score;
            center = actual_center;
            position_center = i;
        }
    }
    dna_s.erase(dna_s.begin() + position_center);


    for(int i=0; i<dna_s.size(); i++){
        pair<string,string> actual_alignment = global_neddleman::get_global_alignment(center, dna_s[i]);
        results.push_back(actual_alignment);
    }
    auto maximum_alignment = max_element(results.begin(), results.end(), choose_maximum_size);
    int maximum_size_alignment = maximum_alignment->first.length();

    for(int i=0; i<results.size(); i++){
        int size_a = results[i].first.length();
        results[i].first += string(maximum_size_alignment-size_a, '-');
        results[i].second += string(maximum_size_alignment-size_a, '-');
        cout << results[i].first << endl;
        cout << results[i].second << endl;

    }

    return 0;
}
