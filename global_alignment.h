#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>

#define LEFT 3
#define UP 1
#define DIAGONAL 2
#define END 0

using namespace std;

namespace global_neddleman {
int **_similarity_matrix;
vector<pair<int,int>> **_position_matrix;
unsigned char **_position_matrix_short;

int m,n;


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

void print_short_matrix_positions(unsigned char **_position_matrix_arg, int m, int n){
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << (int) _position_matrix_arg[i][j] << '\t';
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


void traceback(int i, int j, string dna_a, string dna_b, string tmp_result, string tmp_result_b, pair<string,string> &ret_pair){
    if( i < 1 && j < 1){
        std::reverse(tmp_result.begin(), tmp_result.end());
        std::reverse(tmp_result_b.begin(), tmp_result_b.end());
        ret_pair = make_pair(tmp_result, tmp_result_b);
        return;
    }

    pair<int, int> current_pair;

    switch (_position_matrix_short[i][j]) {
        case UP:
            current_pair = make_pair(i-1,j);
            break;
        case LEFT:
            current_pair = make_pair(i, j-1);
            break;
        case DIAGONAL:
            current_pair = make_pair(i-1,j-1);
            break;
        default:
            break;
    }

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
    traceback(current_pair.first, current_pair.second, dna_a, dna_b, tmp_result, tmp_result_b, ret_pair);
    tmp_result = tmp_result.substr(0, tmp_result.size()-1);
    tmp_result_b = tmp_result_b.substr(0, tmp_result_b.size()-1);


}

string read_sequence_from_file(const char *file_name){
    ifstream ifs(file_name);
    std::string content((std::istreambuf_iterator<char>(ifs)),(std::istreambuf_iterator<char>()));
    return content;
}

int get_maximum_score(string dna_a, string dna_b, bool only_score=false){
    int _length_of_a = dna_a.length();
    int _length_of_b = dna_b.length();

    m = _length_of_a+1;
    n = _length_of_b+1;

    int diagonal_value, up_value, left_value, selected;

    _similarity_matrix = new int*[_length_of_a + 1];

    if(!only_score){
        _position_matrix_short = new unsigned char*[_length_of_a + 1];
    }


    //Algorithm init
    for (int i = 0; i < (_length_of_a + 1); i++)
    {
        _similarity_matrix[i] = new int[_length_of_b + 1];
        if(!only_score) _position_matrix_short[i] = new unsigned char[_length_of_b + 1];
    }

    for (int i = 0; i < _length_of_a; i++)
    {
        for (int j = 0; j < _length_of_b; j++)
        {
            _similarity_matrix[i][j] = 0;
        }
    }

    if(!only_score){
        _position_matrix_short[0][0] = END;
    }

    int gap_penalty = -2;
    for (int i = 1; i < (_length_of_a+1); i++)
    {
        _similarity_matrix[i][0] = i * gap_penalty;
        if(!only_score){
            _position_matrix_short[i][0] = UP;
        }
    }

    for (int j = 1; j < (_length_of_b+1); j++)
    {
        _similarity_matrix[0][j] = j * gap_penalty;
        if(!only_score){
            _position_matrix_short[0][j] = LEFT;
        }
    }

    for (int i = 1; i < _length_of_a + 1; i++)
    {
        for (int j = 1; j < _length_of_b + 1; j++)
        {
            char target_position;

            diagonal_value = _similarity_matrix[i - 1][j - 1] + score_beetween_char(dna_a[i-1], dna_b[j-1]) ;
            up_value = _similarity_matrix[i - 1][j] + gap_penalty;
            left_value = _similarity_matrix[i][j - 1] + gap_penalty;

            selected = std::max(diagonal_value, up_value);
            selected = std::max(selected, left_value);

            if(selected == diagonal_value) target_position = DIAGONAL;
            if(selected == up_value) target_position = UP;
            if(selected == left_value) target_position = LEFT;

            _similarity_matrix[i][j] = selected;
            if(!only_score) _position_matrix_short[i][j] = target_position;
        }
    }
    return _similarity_matrix[_length_of_a][_length_of_b];
}

int get_maximum_score_optimized(string dna_a, string dna_b, bool only_score=false){
    int _length_of_a = dna_a.length();
    int _length_of_b = dna_b.length();

    int diagonal_value, up_value, left_value, selected;

    _similarity_matrix = new int*[2];

    _similarity_matrix[0] = new int[_length_of_b + 1];
    _similarity_matrix[1] = new int[_length_of_b + 1];

    //Algorithm init
    if(!only_score){
        _position_matrix_short = new unsigned char*[_length_of_a + 1];
        for (int i = 0; i < (_length_of_a + 1); i++)
        {
            _position_matrix_short[i] = new unsigned char[_length_of_b + 1];
        }
        _position_matrix_short[0][0] = END;
    }

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < _length_of_b; j++)
        {
            _similarity_matrix[i][j] = 0;
        }
    }



    int gap_penalty = -2;
    for (int i = 1; i < (_length_of_a+1); i++)
    {
        if(!only_score) _position_matrix_short[i][0] = UP;
    }

    for (int j = 1; j < (_length_of_b+1); j++)
    {
        _similarity_matrix[0][j] = j * gap_penalty;
        if(!only_score) _position_matrix_short[0][j] = LEFT;
    }

    int current_row, previous_row;
    int i;
    for (i = 1; i < _length_of_a+1; i++)
    {
        current_row = i%2;
        previous_row = !current_row;
        _similarity_matrix[current_row][0] =  _similarity_matrix[previous_row][0] - 2;
        for (int j = 1; j < _length_of_b + 1; j++)
        {
            char target_position;

            diagonal_value = _similarity_matrix[previous_row][j - 1] + score_beetween_char(dna_a[i-1], dna_b[j-1]) ;
            up_value = _similarity_matrix[previous_row][j] + gap_penalty;
            left_value = _similarity_matrix[current_row][j - 1] + gap_penalty;

            selected = std::max(diagonal_value, up_value);
            selected = std::max(selected, left_value);

            if(selected == diagonal_value) target_position = DIAGONAL;
            if(selected == up_value) target_position = UP;
            if(selected == left_value) target_position = LEFT;

            _similarity_matrix[current_row][j] = selected;
            if(!only_score) _position_matrix_short[i][j] = target_position;
        }
    }

    return _similarity_matrix[current_row][_length_of_b];
}

void delete_pointers(int m, bool only_score=false){
    for(int i = 0; i < 2; i++) {
        delete _similarity_matrix[i];
    }
    delete[] _similarity_matrix;
    if(!only_score){
        for(int i=0; i < m; i++){
            delete _position_matrix_short[i];
        }
        delete[] _position_matrix_short;
    }
}

void old_delete_pointers(int m, bool only_score=false){
    for(int i = 0; i < m; i++) {
        delete _similarity_matrix[i];
        if(!only_score) delete _position_matrix_short[i];
    }
    if(!only_score) delete[] _position_matrix_short;
    delete[] _similarity_matrix;
}

pair<string, string> get_global_alignment(string dna_a, string dna_b){
    pair<string, string> ret;
    int _length_of_a = dna_a.length();
    int _length_of_b = dna_b.length();
    int score = get_maximum_score(dna_a, dna_b);
    string tmp_result = "";
    string tmp_result_2 = "";
    traceback(_length_of_a, _length_of_b, dna_a, dna_b, tmp_result, tmp_result_2, ret);
    old_delete_pointers(_length_of_a+1, false);
    return ret;
}

pair<string, string> get_global_alignment_optimized(string dna_a, string dna_b){
    pair<string, string> ret;
    int _length_of_a = dna_a.length();
    int _length_of_b = dna_b.length();
    int score = get_maximum_score_optimized(dna_a, dna_b);
    string tmp_result = "";
    string tmp_result_2 = "";
    traceback(_length_of_a, _length_of_b, dna_a, dna_b, tmp_result, tmp_result_2, ret);
    delete_pointers(_length_of_a+1, false);
    return ret;
}

}

#endif // GLOBAL_ALIGNMENT_H
