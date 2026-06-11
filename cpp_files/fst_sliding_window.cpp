#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>
#include <omp.h>

using std::cout;
using std::vector;
using std::string;
using std::fstream;
using std::stringstream;
using std::getline;
using std::stoi;
using std::stof;
using std::random_device;
using std::uniform_int_distribution;

// compile with -lz & -fopenmp


/*
    fst sliding window from fst files generated from popoolation2
    standard using: 17
*/

// make global objs
vector<float> weights;
int sigma = 150000;


struct Chrom{
    vector<float> fsts;
    vector<int> positions;
    string chrom;
    void clear(){
        fsts.clear();
        positions.clear();
    }
    // construct
    Chrom(){};
    // destruct
    ~Chrom(){};
    // funcs
    void add_stat(float fst, int pos){
        fsts.push_back(fst);
        positions.push_back(pos);
    }
    void add_chrom(string c){
        chrom = c;
    }
};

struct Window{

    float fst;
    float smth_fst;
    int num_snps; // number of snps in window

    Window(float f, float s, int n){
        fst = f;
        smth_fst = s;
        num_snps = n;
    }

    ~Window(){};

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
            pos_val = chr.fsts[pos];
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
        outfh << chr.chrom << '\t' << chr.positions[i] << '\t'
              << chr.fsts[i] << '\t' <<  snp_cnts[i] << '\t'
              << smooth_container[i] << '\n';
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
    
    if (file_extension == "gz")
    {
        Compressed = true;
    }
    else{
        Compressed = false;
    }

    return Compressed;
}

int load_fsts(string &smth_file, vector<float> &fst_pool){
    
    // objects
    fstream fh;
    string line, tmp;
    float fst;
    int snp_cnt = 0;

    // now do the reading
    fh.open(smth_file, std::ios::in);
    if (fh.is_open()){
        while (getline(fh, line)){
            snp_cnt++;
            stringstream ss(line);
            for (int i = 0; i < 5; i++){
                getline(ss, tmp, '\t');
            }
            fst = stof(tmp);
            fst_pool.push_back(fst);
        }
    }

    fh.close();

    return snp_cnt;
}


float gen_mean(vector<float> &fst_pool, random_device &RD, uniform_int_distribution<int> &UD, int &num_samples){

    vector<float> sampled_fsts(num_samples);
    double val = 0.0;
    int rand_pos;
    float fmean;

    for (int i = 0; i < num_samples; i++){
        rand_pos = UD(RD);
        sampled_fsts[i] = fst_pool[rand_pos];
    }
    
    for (float &f : sampled_fsts){
        val += f;
    }

    if (val != 0.0){
        fmean = static_cast<float>(val) / static_cast<float>(num_samples);
    }
    else{
        fmean = 0.0;
    }
    
    return fmean;
}


Window ParseLine(string &line){

    string tmp;
    stringstream ss(line);
    float fst, smth_fst;
    int num_snps;
    
    for (int i = 0; i < 4; i++){
        getline(ss, tmp, '\t');
        if (i == 2){
            fst = stof(tmp);
        }
        else if (i == 3){
            num_snps = stoi(tmp);
        }
    }
    getline(ss, tmp, '\n'); // go to last column
    smth_fst = stof(tmp);

    Window W = Window(fst, smth_fst, num_snps);

    return W;
}

float estimate_pval(float &fst, int &bootsamples, vector<float> &means){

    float pval;
    int rank = 0;

    // sort than rank
    std::sort(means.begin(), means.end(), [](float a, float b) {return a < b;});

    float mymin = means[0];
        
    for (int i = 0; i < bootsamples; i++){
        if (means[i] < fst){
            rank++;
        }
        else{
            break;
        }
    }
    
    float percent = static_cast<float>(rank) / static_cast<float>(bootsamples);
    
    pval = 1.0 - percent;

    return pval;
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

    int buffer_line_size = 100, position; // read in size
    struct gzFile_s *gzf;      
    char line_buffer[buffer_line_size]; // char array to hold line
    gzf = gzopen(fst_file, "r");  
    bool eol = false;                        
    memset(line_buffer, '\0', buffer_line_size); 
    string line, chrom, tmp, position_s, cur_chrom;
    float fst;
    vector<int> positions;
    vector<float> fsts;
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
            getline(ss, position_s, '\t');
            getline(ss, tmp, '='); // todo : do a for loop in the future
            getline(ss, tmp, 'n');
            fst = stof(tmp);  
            position = stoi(position_s);
            if (cur_chrom == ""){
                chr.add_stat(fst, position);
                chr.add_chrom(chrom);
                cur_chrom = chrom;
            }
            else if (cur_chrom == chrom){
                chr.add_stat(fst, position);
            }
            else{
                write_smooth_array(chr, outfh);
                chr.clear();
                chr.add_stat(fst, position);
                chr.add_chrom(chrom);
                cur_chrom = chrom;
            }
        }
     } while (!gzeof(gzf) && !eol);

     if (!chr.fsts.empty()){
        write_smooth_array(chr, outfh);
        chr.clear();
     }
}



void Bootstrap(int &bootstraps, string &smth_file, vector<float> &fst_pools, int &snp_count, int &num_threads){

    // generate random generators
    std::random_device RD;
    unsigned seed = RD();
    std::uniform_int_distribution<int> UD(0, snp_count);

    // objects to parse each line
    fstream fh;
    fh.open(smth_file, std::ios::in);
    string line;

    // outfile
    string outname = smth_file.substr(smth_file.find_last_of("/\\") + 1);
    outname = "Bootstrapped_" + outname;
    fstream ofh;
    ofh.open(outname, std::ios::out);

    // threading
    int tcount = 0;
    vector<string> lines(num_threads);
    vector<float> results(num_threads);
    omp_set_num_threads(num_threads);

    if (fh.is_open()){
        while (getline(fh, line)){
            lines.push_back(line);
            tcount++;
            if (tcount < num_threads){
                continue;
            }
            else{
                #pragma omp parallel for{
                for (int t = 0; t < num_threads; t++){
                    Window W = ParseLine(lines[t]);
                    vector<float> means(bootstraps);
                    for (int i = 0; i < bootstraps; i++){
                        means[i] = gen_mean(fst_pools, RD, UD, W.num_snps);
                    }
                    results[t] = estimate_pval(W.fst, bootstraps, means);
                    }     
                };
                for (int k = 0; k < num_threads; k++){
                    ofh << lines[k] << '\t' << std::setprecision(5) << results[k] << '\n';
                }
                tcount = 0;
                results.clear();
                lines.clear();
            }
            }
        if (!lines.empty()){
            #pragma omp parallel for{
            for (int t = 0; t < num_threads; t++){
                Window W = ParseLine(lines[t]);
                vector<float> means(bootstraps);
                for (int i = 0; i < bootstraps; i++){
                    means[i] = gen_mean(fst_pools, RD, UD, W.num_snps);
                }
                results[t] = estimate_pval(W.fst, bootstraps, means);
            }
            };
            for (int k = 0; k < num_threads; k++){
                ofh << lines[k] << '\t' << std::setprecision(5) << results[k] << '\n';
            }
            tcount = 0;
            results.clear();
            lines.clear();
        }
        
        }

    fh.close();
    ofh.close();   
}


int main(int argc, char *args[])
{
    // params
    string fst_file, output_file, message, file_basename, smth_file;
    bool compressed;
    int bootstraps = 1000, num_threads = 1;

    // in case of help or error
    message = "./fst_sliding_window -f fst_file -S Smoothed_Res -b BOOTNUM\n";

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
        else if (param == "-f" || param == "--fst_file")
        {
            fst_file = string(args[i + 1]);
        }
        else if (param == "-b" || param == "--bootstraps")
        {
            bootstraps = stoi(args[i + 1]);
        }
        else if (param == "-s" || param == "--sigma")
        {
            sigma = stoi(args[i + 1]);
        }
        else if (param == "-S" || param == "--smth_file")
        {
            smth_file = string(args[i + 1]);
        }
        else if (param == "-t" || param == "--threads")
        {
            num_threads = stoi(args[i + 1]);
        }
    }

    if (fst_file == "" && smth_file == "")
    {
        cout << "Neither fst or smoothed fst file supplied.\n";
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
    if (fst_file != ""){
        // check if gzipped b/c that's how i wrote the program
        compressed = check_for_uncompressed(fst_file);
        if (compressed){
            if (smth_file == ""){
                const char *ffile = fst_file.c_str();
                if (output_file == ""){
                    output_file = rename_file(fst_file);
                 }
                create_smth_file(ffile, output_file);
                cout << "Created " << output_file << "\nRerun with new output"
                    << " instead of " << fst_file << " using -S or --smth_file\n";
                exit(EXIT_SUCCESS);  
            }
        }
        else{
            std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
            exit(EXIT_FAILURE);
        }
    }
    else if (smth_file != ""){
        vector<float> fst_pool;
        int snp_count = load_fsts(smth_file, fst_pool);
        Bootstrap(bootstraps, smth_file, fst_pool, snp_count, num_threads);
        
    }
    
    return 0;
}