#include <iostream>
#include <cstring>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <math.h>

using std::string;
using std::vector;
using std::cout;
using std::fstream;

/*
Had to compile using the following flags -lz
*/

// sigma & weights will be global
const int SIGMA = 150000;
const int window_span = (3 * SIGMA) + 1;
const int buffer_line_size = 1024; 
double weights[window_span];

void calc_weights(double (&weights)[window_span]){

    // fill weights once
    for (int i = 0; i < window_span; i++){
        weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(SIGMA, 2)));
    }
}

struct Chrom {

    string name = "";
    vector<int>    positions;
    vector<double> qvals;
    vector<double> smoothed;

};

struct Entry {

    string chrom;
    int    position;
    double qval;

};

int
parse_line(Entry &e, const string &line){

    vector<string> parts;
    parts.reserve(10);

    uint i = 0;
    uint j = 0;
    uint dist = 0;
    const uint len = line.size();

    while (i < len && line[i] != '\n'){
        while (j < len && line[j] != '\t' && line[j] != '\n'){
            j++;
        }
        dist = j - i;
        parts.push_back(line.substr(i, dist));
        j++;
        i = j;
    }

    e.chrom    = parts[0];
    e.position = std::stoi(parts[1]);

    // get around NA's
    try {
        e.qval = std::stod(parts[9]);
    }
    catch (const std::invalid_argument &ex) {
        e.qval = 0.0;
    }
    

    return 0;
}


void 
smooth_cmh(Chrom &current_chrom, fstream &filestream){
    
    /*
    positions = bp positions
    fsts = fsts at each position
    weights = precomputed kernel smoothed weights
    start_pos = beginning of window on chrom
    end_pos = end of window on chrom
    window_cen = current position on chrom
    window_start = used to find start_pos by adjusting from center
    window_end = used to find end_pos by adjusting to window_start
    n_window_sites = number of positions
    distnace = distance from window center
    window_avg = will hold weighted average
    window_weight = will accumulate weights to use for weighing average
    weight = kernel weight at position
    */

    int start_pos = 0, end_pos = 0;
    int window_cen, window_start, window_end, n_window_sites, distance, pos_num;
    double qval;
    double avg;
    double window_weight, weight;
    const int n_sites = current_chrom.positions.size();

    if (n_sites == 0){
        return;
    }

    for (int site = 0; site < n_sites; site++){
        window_cen   = current_chrom.positions[site]; // go to center
        window_start = window_cen - (window_span / 2); // go to start
        if (window_start < 0){ 
            window_start = 1; //keep within bound
        }
        window_end = window_start + window_span; // go to end
        while (start_pos < n_sites && current_chrom.positions[start_pos] < window_start){ // start pointer in positions
            start_pos++;
        }
        while (end_pos < n_sites && current_chrom.positions[end_pos] < window_end){ // end pointer in positions
            end_pos++;
        }
        n_window_sites = 0;   // window size
        avg            = 0.0; // reset average
        window_weight  = 0.0;
        for (int pos = start_pos; pos < end_pos; pos++){
            n_window_sites++;
            qval = current_chrom.qvals[pos];           // qvalue
            pos_num  = current_chrom.positions[pos];   // get current position
            distance = (int)abs(pos_num - window_cen); // relative to center
            weight   = weights[distance];              // get weight at position
            avg      += (qval * weight);               // accumulate sum
            window_weight += weight;                   // accumulate weight
        }
        if (n_window_sites == 0){
            exit(EXIT_FAILURE);
        }
        current_chrom.smoothed.push_back((avg/window_weight));
    }

    // now to write to output
    for (int i = 0; i < n_sites; i++){
        filestream << current_chrom.name        << '\t' << current_chrom.positions[i] << '\t' 
                   << current_chrom.smoothed[i] << '\n'; 
    }

}

bool check_for_uncompressed(const string &file_name){

    // get basename of file
    string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // get file extension
    string file_extension =
        file_basename.substr(file_basename.find_last_of('.') + 1);

    return file_extension == "gz";
}

string 
make_output_name(const char *in_file){

    string outname;
    string file_name = in_file;

    // get basename of file
    string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add suffix
    outname = "Smoothed_" + file_basename.substr(0, file_basename.find_last_of('.'));
    return outname;
}
int
update_chrom(Chrom &current_chrom, Entry &e){

    // clear all the vectors

    current_chrom.positions.clear();
    current_chrom.qvals.clear();
    current_chrom.smoothed.clear();

    // update
    current_chrom.name = e.chrom;
    current_chrom.positions.push_back(e.position);
    current_chrom.qvals.push_back(e.qval);

    return 0;
}

void 
parse_compressed(const char *fst_file){

    struct gzFile_s *gzf;                        // pointer to file
    char line_buffer[buffer_line_size];          // char array to hold line
    gzf = gzopen(fst_file, "r");               // file handle
    bool eol = false;                            // end of file
    string line;                              // char arr -> string
    memset(line_buffer, '\0', buffer_line_size); // fill in nulls
    string outfile = make_output_name(fst_file);
    fstream file_stream;
    file_stream.open(outfile, std::ios::out);
    Chrom current_chrom;
    Entry e;

    // fill weights
    calc_weights(weights);

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
            line = string(line_buffer);
            parse_line(e, line);

            if (current_chrom.name == ""){
                current_chrom.name = e.chrom;
            }
            if (current_chrom.name == e.chrom){
                current_chrom.positions.push_back(e.position);
                current_chrom.qvals.push_back(e.qval);
            }
            else {
                smooth_cmh(current_chrom, file_stream);
                update_chrom(current_chrom, e);
            }
    
        }
     } while (!gzeof(gzf) && !eol);

    //
    // get last chrom
    //
    smooth_cmh(current_chrom, file_stream);

    file_stream.close();
}

int main(int argc, char *args[])
{
    // params
    string fst_file;
    string output_file;
    string message;
    string file_basename;
    bool compressed;

    // in case of help or error
    message = "./kernel_smooth_cmh -f input_file [required]\n";

    // read args
    std::string param;

    for (int i = 0; i < argc; i++)
    {
        param = std::string(args[i]);

        if (param == "-h" || param == "--help")
        {
            cout << message;
            return 0;
        }
        else if (param == "-f" || param == "--file")
        {
            fst_file = std::string(args[i + 1]);
        }
    }

    if (fst_file == "")
    {
        cout << "No file supplied.\n" << message;
        exit(EXIT_FAILURE);
    }

    // check if gzipped b/c that's how i wrote the program
    compressed = check_for_uncompressed(fst_file);

    // start processing
    if (compressed){
        const char *ffile = fst_file.c_str();
        parse_compressed(ffile);
    }
    else{
        std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(EXIT_FAILURE);
    }
    
    return 0;
}