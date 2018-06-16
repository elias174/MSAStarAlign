#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>


#define LEFT -1
#define UP 1
#define DIAG 0

using namespace std;

namespace global_neddleman {
    int **_similarity_matrix;
    vector<pair<int,int>> **_position_matrix;
    vector<pair<string, string>> results;

    int score_beetween_char(char a, char b){
        if(a!=b)
            return -1;
        return 1;
    }

    void print_matrix(int **_similarity_matrix, int m, int n){
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << _similarity_matrix[i][j] << '\t';
            }
            cout << endl;
        }
    }

    void print_vector_pairs(vector<pair<int,int>> pairs){
        for(int i=0; i<pairs.size(); i++){
            cout << "(" << pairs[i].first << "," << pairs[i].second << "),";
        }
    }

    void print_matrix_positions(vector<pair<int,int>> **_position_matrix, int m, int n){
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                print_vector_pairs(_position_matrix[i][j]);
                cout << '\t' << '\t';
            }
            cout << endl;
        }
    }

    void write_to_file(string result_a, string result_b){
        ofstream m_file;
        m_file.open("results.txt", std::ofstream::out | std::ofstream::app);
        m_file << result_a.c_str() << endl;
        m_file << result_b.c_str() << endl << endl << endl;
    }


    void traceback(int i, int j, string dna_a, string dna_b, string tmp_result, string tmp_result_b){
        if( i < 1 && j < 1){
            std::reverse(tmp_result.begin(), tmp_result.end());
            std::reverse(tmp_result_b.begin(), tmp_result_b.end());
            //        cout << tmp_result << '\t';
            //        cout << tmp_result_b << endl;
            results.push_back(make_pair(tmp_result, tmp_result_b));
            return;
        }
        for(int k=0; k<_position_matrix[i][j].size(); ++k){
            pair<int, int> current_pair = _position_matrix[i][j][k];
            if( i > current_pair.first && j > current_pair.second){
                tmp_result += dna_a[i-1];
                tmp_result_b += dna_b[j-1];
            }
            else if( i > current_pair.first && j == current_pair.second){
                tmp_result += dna_a[i-1];
                tmp_result_b += '-';

            }
            else if( i == current_pair.first && j > current_pair.second){
                tmp_result += '-';
                tmp_result_b += dna_b[j-1];
            }
            traceback(current_pair.first, current_pair.second, dna_a, dna_b, tmp_result, tmp_result_b);
            tmp_result = tmp_result.substr(0, tmp_result.size()-1);
            tmp_result_b = tmp_result_b.substr(0, tmp_result_b.size()-1);
        }

    }

    string read_sequence_from_file(const char *file_name){
        ifstream ifs(file_name);
        std::string content((std::istreambuf_iterator<char>(ifs)),(std::istreambuf_iterator<char>()));
        return content;
    }



    int get_maximum_score(string dna_a, string dna_b){
        int _length_of_a = dna_a.length();
        int _length_of_b = dna_b.length();

        int diagonal_value, up_value, left_value, selected;

        _similarity_matrix = new int*[_length_of_a + 1];
        _position_matrix = new vector<pair<int,int>>*[_length_of_a + 1];


        //Algorithm init
        for (int i = 0; i < (_length_of_a + 1); i++)
        {
            _similarity_matrix[i] = new int[_length_of_b + 1];
            _position_matrix[i] = new vector<pair<int,int>>[_length_of_b + 1];
        }

        for (int i = 0; i < _length_of_a; i++)
        {
            for (int j = 0; j < _length_of_b; j++)
            {
                _similarity_matrix[i][j] = 0;
            }
        }

        vector<pair<int,int>> init_vector;
        init_vector.push_back(make_pair(-1, -1));
        _position_matrix[0][0] = init_vector;

        int gap_penalty = -2;
        for (int i = 1; i < (_length_of_a+1); i++)
        {
            _similarity_matrix[i][0] = i * gap_penalty;
            vector<pair<int,int>> to_matrix;
            to_matrix.push_back(make_pair(i-1, 0));
            _position_matrix[i][0] = to_matrix;
        }

        for (int j = 1; j < (_length_of_b+1); j++)
        {
            _similarity_matrix[0][j] = j * gap_penalty;
            vector<pair<int,int>> to_matrix;
            to_matrix.push_back(make_pair(0, j-1));
            _position_matrix[0][j] = to_matrix;
        }

        int counter_matrix_score = 1;

        for (int i = 1; i < _length_of_a + 1; i++)
        {
            for (int j = 1; j < _length_of_b + 1; j++)
            {
                vector<pair<int,int>> positions;
                pair<int, int> diagonal = make_pair(i-1, j-1);
                pair<int, int> up = make_pair(i-1, j);
                pair<int, int> left = make_pair(i, j-1);

                diagonal_value = _similarity_matrix[i - 1][j - 1] + score_beetween_char(dna_a[i-1], dna_b[j-1]) ;
                up_value = _similarity_matrix[i - 1][j] + gap_penalty;
                left_value = _similarity_matrix[i][j - 1] + gap_penalty;

                selected = std::max(diagonal_value, up_value);
                selected = std::max(selected, left_value);

                if(selected == diagonal_value) positions.push_back(diagonal);
                if(selected == up_value) positions.push_back(up);
                if(selected == left_value) positions.push_back(left);

                _similarity_matrix[i][j] = selected;
                _position_matrix[i][j] = positions;
            }
        }

        return _similarity_matrix[_length_of_a][_length_of_b];
    }

    pair<string, string> get_global_alignment(string dna_a, string dna_b){
        results.clear();
        int _length_of_a = dna_a.length();
        int _length_of_b = dna_b.length();
        int score = get_maximum_score(dna_a, dna_b);
        string tmp_result = "";
        string tmp_result_2 = "";
        traceback(_length_of_a, _length_of_b, dna_a, dna_b, tmp_result, tmp_result_2);
        return results[0];
    }

}

#endif // GLOBAL_ALIGNMENT_H
