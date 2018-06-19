#include <iostream>
#include <limits>
#include <vector>
#include <fstream>

#include "global_alignment.h"


using namespace std;

vector<pair<string,string>> results;

bool choose_maximum_size(const std::pair<string, string> &lhs, const std::pair<string,string> &rhs) {
    return lhs.first.length() < rhs.first.length() ;
}

int main(int argc, char **argv)
{
    vector<string> dna_s;

    if(argc<2){
        perror("No dataset file specified");
        return 0;
    }

    std::ifstream file(argv[1]);
    std::string str;
    while (std::getline(file, str))
    {
        dna_s.push_back(str);
    }

    int maximum = numeric_limits<int>::max() * -1;
    string center;
    int position_center;
    for(int i=0; i<dna_s.size(); i++){
        string actual_center = dna_s[i];
        int acum_score = 0;
        for(int j=0; j<dna_s.size(); j++){
            if(i!=j){
                acum_score += global_neddleman::get_maximum_score(actual_center, dna_s[j], true);
                global_neddleman::old_delete_pointers(actual_center.length()+1, true);
                //cout << "deleted" << endl;
            }
        }
        if(acum_score > maximum){
            maximum = acum_score;
            center = actual_center;
            position_center = i;
        }
        //cout << i << endl;
    }

    dna_s.erase(dna_s.begin() + position_center);

//    cout << "finish center find " << endl;
//    cout << endl;

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
