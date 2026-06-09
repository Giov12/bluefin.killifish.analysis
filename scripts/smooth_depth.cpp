#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <algorithm>
#include <math.h>

using std::cout;
using std::vector;
using std::string;
using std::fstream;
using std::stringstream;
using std::getline;
using std::stoi;
using std::stof;


// compile with -lz


/*
    smooth each depth position for the pool-seq results
    standard using: 17
*/

// make global objs
vector<float> weights;
int sigma = 150000;


struct Chrom{
    vector<int> depths;
    vector<int> positions;
    string chrom;
    void clear(){
        depths.clear();
        positions.clear();
    }
    // construct
    Chrom(){};
    // destruct
    ~Chrom(){};
    // funcs
    void add_stat(int depth, int pos){
        depths.push_back(depth);
        positions.push_back(pos);
    }
    void add_chrom(string c){
        chrom = c;
    }
};

void calc_weights(int &window_size, int &sigma){

    int window_span = (3 * sigma) + 1;

    // resize
    weights.resize(window_span);

    // fill weights once
    for (int i = 0; i < window_span; i++){
        weights[i] = (float)exp((-1 * pow(i, 2)) / (2 * pow(sigma, 2)));
    }
}

void write_smooth_array(Chrom &chr, fstream &outfh){

    int start_pos = 0, end_pos = 0;
    int window_cen, window_start, window_end, n_window_sites, distance, pos_num;
    float window_avg = 0.0, window_weight = 0.0, weight, pos_val;
    int n_sites = chr.positions.size();
    int window_span = 3 * sigma;
    vector<float> smooth_container(n_sites);
    vector<int> snp_cnts(n_sites);

    for (int site = 0; site < n_sites; site++){
        window_cen = chr.positions[site]; // go to center
        window_start = window_cen - window_span; // go to start
        if (window_start < 0){ 
            window_start = 1; //keep within bound
        }
        window_end = window_cen + window_span; // go to end
        while (chr.positions[start_pos] < window_start){ // start pointer in positions
            start_pos++;
        }
        while (end_pos < n_sites && chr.positions[end_pos] < window_end){ // end pointer in positions
            end_pos++;
        }
        n_window_sites = 0; // window size
        window_weight = 0.0;
        for (int pos = start_pos; pos < end_pos; pos++){
            n_window_sites++;
            pos_val = chr.depths[pos];
            pos_num = chr.positions[pos]; // get current position
            distance = (int)abs(pos_num - window_cen); // relative to center
            weight = weights[distance]; // get weight at position
            window_avg += (pos_val * weight); // accumulate sum
            window_weight += weight; // accumulate weight
        }
        if (n_window_sites == 0){
            continue;
        }
        if (window_avg != 0.0){
            window_avg = window_avg / window_weight;
        }
        if (isnan(window_avg)){
            continue;
        }
        smooth_container[site] = window_avg;
        snp_cnts[site] = n_window_sites;
    }

    for (int i = 0; i < n_sites; i++){
        outfh << chr.chrom     << '\t' << chr.positions[i] << '\t'
              << chr.depths[i] << '\t' <<  snp_cnts[i]     << '\t'
              << smooth_container[i]   << '\n';
    }
}

bool check_for_uncompressed(string file_name){

    bool Compressed;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // get file extension
    std::string file_extension =
        file_basename.substr(file_basename.find_last_of('.') + 1);
    
    return (file_extension == "gz");
}


string rename_file(string &file_name){

    string outname;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add suffix
    outname = "Smoothed_" + file_basename.substr(0, file_basename.find_last_of('.'));
    return outname;
}

void create_smth_file(const char *fst_file, string &output){

    int buffer_line_size = 100, position, depth; // read in size
    struct gzFile_s *gzf;      
    char line_buffer[buffer_line_size]; // char array to hold line
    gzf = gzopen(fst_file, "r");  
    bool eol = false;                        
    memset(line_buffer, '\0', buffer_line_size); 
    string line, chrom, tmp, cur_chrom;
    vector<int> positions;
    vector<float> depths;
    fstream outfh;
    outfh.open(output, std::ios::out);
    Chrom chr = Chrom();
    do
    {
        if (gzgets(gzf, line_buffer, buffer_line_size) == NULL)
        {
            break;
        }
        if (line_buffer[0] == '\n')
        {
            continue;
        }
        else
        {
            // convert char array to string
            line = std::string(line_buffer);
            // now to get the fields
            stringstream ss(line);
            getline(ss, chrom, '\t');
            getline(ss, tmp, '\t');
            position = stoi(tmp);
            getline(ss, tmp, '\n'); // todo : do a for loop in the future
            depth = stoi(tmp);
            if (cur_chrom == ""){
                chr.add_stat(depth, position);
                chr.add_chrom(chrom);
                cur_chrom = chrom;
            }
            else if (cur_chrom == chrom){
                chr.add_stat(depth, position);
            }
            else{
                write_smooth_array(chr, outfh);
                chr.clear();
                chr.add_stat(depth, position);
                chr.add_chrom(chrom);
                cur_chrom = chrom;
            }
        }
     } while (!gzeof(gzf) && !eol);

     if (!chr.depths.empty()){
        write_smooth_array(chr, outfh);
        chr.clear();
     }

     outfh.close();
}


int main(int argc, char *args[])
{
    // params
    string depth_file, output_file, message, file_basename;
    bool compressed;

    // in case of help or error
    message = "./depth_file -d file.depth.gz -s int\n";

    // read args
    string param;

    for (int i = 0; i < argc; i++)
    {
        param = string(args[i]);

        if (param == "-h" || param == "--help")
        {
            cout << message;
            return 0;
        }
        else if (param == "-d" || param == "--depth_file")
        {
            depth_file = string(args[i + 1]);
        }
        else if (param == "-s" || param == "--sigma")
        {
            sigma = stoi(args[i + 1]);
        }
    }

    if (depth_file == "")
    {
        cout << "depth file is required\n";
        cout << message;
        exit(EXIT_FAILURE);
    }
    if (sigma <= 0){
        string e = "Sigma must be greater than 0\n";
        cout << e;
        exit(EXIT_FAILURE);
    }

    
    // create weights
    int window_span = (3 * sigma) + 1;
    calc_weights(window_span, sigma);
    
    
    // start processing
    if (depth_file != ""){
        // check if gzipped b/c that's how i wrote the program
        compressed = check_for_uncompressed(depth_file);
        if (compressed){
                const char *ffile = depth_file.c_str();
                if (output_file == ""){
                    output_file = rename_file(depth_file);
                 }
                create_smth_file(ffile, output_file);
                exit(EXIT_SUCCESS);  
            }
        else{
            std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
            exit(EXIT_FAILURE);
        }
    }
    
    return 0;
}