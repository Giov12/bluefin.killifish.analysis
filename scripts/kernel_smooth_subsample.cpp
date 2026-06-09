#include <iostream>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <math.h>

/*
Had to compile using the following flags -lz
*/

// sigma & weights will be global
const int SIGMA = 150000;
const int window_span = (3 * SIGMA) + 1;
double weights[window_span];

void calc_weights(double (&weights)[window_span]){

    // fill weights once
    for (int i = 0; i < window_span; i++){
        weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(SIGMA, 2))); }
}

void create_smooth_array(std::vector<int> &positions, std::vector<int> &depths, const int n_sites, double *smth_ar, std::string chrom, std::string outfile_n){
    
    /*
    positions = bp positions
    depths = depths at each position
    weights = precomputed kernel smoothed weights
    start_pos = beginning of window on chrom
    end_pos = end of window on chrom
    window_cen = current position on chrom
    window_start = used to find start_pos by adjusting from center
    window_end = used to find end_pos by adjusting to window_start
    n_window_sites = number of positions
    distnace = distance from window center
    pos_val = depth at position
    pos_num = position in bp
    window_avg = will hold weighted average
    window_weight = will accumulate weights to use for weighing average
    weight = kernel weight at position
    */

    int start_pos = 0, end_pos = 0;
    int window_cen, window_start, window_end, n_window_sites, distance, pos_val, pos_num;
    double window_avg, window_weight, weight;

    for (int site = 0; site < n_sites; site++){
        window_cen = positions[site]; // go to center
        window_start = window_cen - window_span; // go to start
        if (window_start < 0){ 
            window_start = 1; //keep within bound
        }
        window_end = window_start + window_span; // go to end
        while (positions[start_pos] < window_start){ // start pointer in positions
            start_pos++;
        }
        while (end_pos < n_sites && positions[end_pos] < window_end){ // end pointer in positions
            end_pos++;
        }
        n_window_sites = 0; // window size
        window_avg = 0.0; // reset averages
        window_weight = 0.0;
        for (int pos = start_pos; pos < end_pos; pos++){
            n_window_sites++;
            pos_val = depths[pos]; // get depth
            pos_num = positions[pos]; // get current position
            distance = (int)abs(pos_num - window_cen); // relative to center
            weight = weights[distance]; // get weight at position
            window_avg += (pos_val * weight); // accumulate sum
            window_weight += weight; // accumulate weight
        }
        if (n_window_sites == 0){ // do I really need this?
            exit(EXIT_FAILURE); }
        window_avg = window_avg/window_weight; // adjust
        smth_ar[site] = window_avg; // set
    }

    // now to write to output
    std::fstream outfile;
    outfile.open(outfile_n, std::ios_base::app); // append
    for (int i = 0; i < positions.size(); i++){
        outfile << chrom << '\t' << positions[i] << '\t' << smth_ar[i] << '\n';
    }
    outfile.close();

}

bool check_for_uncompressed(std::string file_name){

    bool Compressed;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // get file extension
    std::string file_extension =
        file_basename.substr(file_basename.find_last_of('.') + 1);
    
    if (file_extension == "gz") {Compressed = true;}
    else { Compressed = false; }

    return Compressed;

}

std::string rename_file(const char *depth_file, int subsample){

    std::string outname;
    std::string file_name = depth_file;
    std::string subsamp = std::to_string(subsample);

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add suffix
    outname = "SS_" + subsamp + "_" + file_basename.substr(0, file_basename.find_last_of('.'));

    return outname;
}

int median(std::vector<int> &temp){
    
    int median;
    int arraysize = temp.size();

    //sort vector
    std::sort(temp.begin(), temp.end());

    if (arraysize % 2 == 0){ // even
        median = int((temp[arraysize/2] + temp[(arraysize/2) - 1])/2); }
    else{ // odd
        median = temp[arraysize/2];}

    return median;

}

void parse_compressed(const char *depth_file, int subsample){

    int buffer_line_size = 80; // should be enough         
    struct gzFile_s *gzf;                        
    char line_buffer[buffer_line_size]; 
    gzf = gzopen(depth_file, "r"); // why does r have to be a string ?
    bool eol = false;               
    memset(line_buffer, '\0', buffer_line_size); // read up to NULL char
    std::string chrom, depth_s, position_s, cur_chrom, line;
    int depth, position, median_depth;
    std::vector<int> positions;
    std::vector<int> depths;
    // TODO I am 100% sure there's a better way to do this
    std::vector<int> temp(subsample); // hold temp depths
    std::vector<int> temp2(subsample); // hold temp pos
    int temp_size = 0; // should be faster than computing vector.size(), but is it worth it?
    std::string outfile = rename_file(depth_file, subsample);  

    // fill weights
    calc_weights(weights);

    do
    {
        if (gzgets(gzf, line_buffer, buffer_line_size) == NULL) { break;}
        if (line_buffer[0] == '\n'){ continue; }
        else{
            // convert char array to string
            line = std::string(line_buffer);
            // now to get the fields
            std::stringstream ss(line);
            std::getline(ss, chrom, '\t');
            std::getline(ss, position_s, '\t');
            std::getline(ss, depth_s, '\n');
            depth = std::stoi(depth_s);
            position = std::stoi(position_s);
            if (cur_chrom == ""){
                cur_chrom = chrom;
                temp.push_back(depth);
                temp2.push_back(position);
                temp_size++;
            }
            // doing .size() takes time so using a counter
            else if (chrom == cur_chrom && temp_size < subsample){
                temp.push_back(depth);
                temp2.push_back(position);
                temp_size++;
            }
            else if (chrom == cur_chrom && temp_size == subsample){
                median_depth = median(temp);
                depths.push_back(median_depth);
                positions.push_back(position);
                temp.clear();
                temp2.clear();
                temp.push_back(depth);
                temp2.push_back(position);
                temp_size = 1; // restart 
            }
            else if (positions.size() > 0){
                    const unsigned int n_sites = positions.size(); //unsigned to keep +value
                    double * smoothed_vals = new double[n_sites]; // passing pointer
                    create_smooth_array(positions, depths, n_sites, smoothed_vals, cur_chrom, outfile);
                    depths.clear(); // clear objects
                    positions.clear();
                    delete []smoothed_vals;
                    temp_size = 0;
                    temp.clear();
                    cur_chrom = chrom;
                    temp[temp_size] = depth;
                    temp2[temp_size] = position;
                    temp_size++;
                    }     
            else { continue; }  
            } 
        } while (!gzeof(gzf) && !eol);

    // if left over
    if (positions.size() > 0 && temp.size() > 0){
        // add any left over elements from the temp arrays
        median_depth = median(temp);
        depths.push_back(median_depth);
        positions.push_back(position); // last called position
        const unsigned int n_sites = positions.size();
        double * smoothed_vals = new double[n_sites];
        create_smooth_array(positions, depths, n_sites, smoothed_vals, cur_chrom, outfile);
        depths.clear();
        positions.clear();
        delete []smoothed_vals; }
}

int main(int argc, char *args[])
{
    // params
    std::string depth_file;
    std::string output_file;
    std::string message;
    std::string file_basename;
    bool compressed;
    int subsample = 10000; // default value

    // in case of help or error
    message = "./kernel_smooth -d depth_file [required] -s Int [default = " + std::to_string(subsample) + "]\n";

    // read args
    std::string param;

    for (int i = 0; i < argc; i++) {
        param = std::string(args[i]);

        if (param == "-h" || param == "--help")
        {
            std::cout << message;
            return 0;
        }
        else if (param == "-d" || param == "--depth_file"){ 
            depth_file = std::string(args[i + 1]); }
        else if (param == "-s" || param == "--subsample"){
            subsample = std::stoi(args[i + 1]); }
    }

    if (depth_file == ""){
        std::cout << "No depth file supplied.\n";
        std::cout << message;
        exit(EXIT_FAILURE); }

    // check if gzipped b/c that's how i wrote the program
    compressed = check_for_uncompressed(depth_file);

    // start processing
    if (compressed){
        const char *dfile = depth_file.c_str();
        parse_compressed(dfile, subsample);
    }
    else{
        std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(EXIT_FAILURE);
    }
    
    return 0;
}