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

struct Record{
    string chrom;
    int pos;
    float fst;
    float zscore;

    Record(string c, int p, float f){
        chrom = c;
        pos = p;
        fst = f;
    }

    ~Record(){};
};


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


Record ParseLine(string &line){
    
    int pos;
    string tmp, chrom;
    float fst;
    stringstream ss(line);
    
    getline(ss, tmp, '\t');
    chrom = tmp;
    getline(ss, tmp, '\t');
    pos = stoi(tmp);
    getline(ss, tmp, '='); // todo : do a for loop in the future
    getline(ss, tmp, 'n');
    fst = stof(tmp);

    Record R = Record(chrom, pos, fst);

    return R;
}

string rename_file(string &file_name, char &c){

    string outname;
    string prefix;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add prefix
    prefix = "Top_";
    prefix = prefix + c;
    prefix = prefix + "_Per_SNPs_";

    outname = prefix + file_basename.substr(0, file_basename.find_last_of('.'));
    return outname;
}

float get_mean(const vector<Record> &records){

    float sum = 0.0;

    for (auto r: records){
        sum += r.fst;
    }

    return (sum / static_cast<float>(records.size()));
}

float get_SD(const vector<Record> &records, float &mean){

    float varianceSum = 0.0, sd;

    for (auto &r: records){
        varianceSum += (r.fst - mean) * (r.fst - mean);
    }

    float variance = varianceSum / static_cast<float>(records.size());

    sd = sqrt(variance);

    return sd;

}

float get_ZScore(float &fst, float &SD, float &mean){
    return (fst - mean) / SD;
}


void Score_Fsts(vector<Record> &records, float &SD, float &mean, char percent, const char *fst_file){

    float percentage;

    // impute zscores
    for (auto &r: records){
        r.zscore = get_ZScore(r.fst, SD, mean);
    }

    // sort
    std::sort(records.begin(), records.end(), [](const auto &lhs, const auto &rhs) {return lhs.zscore > rhs.zscore;});

    // get size of the array to estimate percentile
    int num_recs = records.size();

    if (percent == '1'){
        percentage = 0.01;
    }
    else if (percent == '5'){
        percentage = 0.05;
    }

    int count = (int)(num_recs * percentage);

    fstream outfh;
    string f = string(fst_file);
    string output = rename_file(f, percent);
    outfh.open(output, std::ios::out);

    for (int i = 0; i < count; i++){
        outfh << records[i].chrom << '\t' << records[i].pos << '\t' << records[i].fst << '\n';
    }

    outfh.close();

}

void Scoring_Pipeline(const char *fst_file){

    int buffer_line_size = 100, position; // read in size
    struct gzFile_s *gzf;      
    char line_buffer[buffer_line_size]; // char array to hold line
    gzf = gzopen(fst_file, "r");  
    bool eol = false;                        
    memset(line_buffer, '\0', buffer_line_size); 
    string line, chrom, tmp, position_s;
    float mean, sd;
    vector<Record> records;

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
            // now to get the fields
            if (line == ""){
                continue;
            }
            Record R = ParseLine(line);
            records.push_back(R);
        }
     } while (!gzeof(gzf) && !eol);

    // get mean and standard deviation
    mean = get_mean(records);
    sd = get_SD(records, mean);

    // write out rankings
    Score_Fsts(records, sd, mean, '1', fst_file);
    Score_Fsts(records, sd, mean, '5', fst_file);

}

int main(int argc, char *args[])
{
    // params
    string fst_file, message, file_basename;
    bool compressed = false;

    // in case of help or error
    message = "./Rank_Fsts -f fst_file [required]\n";

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
    }

    if (fst_file == "")
    {
        cout << "No fst file supplied.\n";
        cout << message;
        exit(EXIT_FAILURE);
    }

    // start processing
    // check if gzipped b/c that's how i wrote the program
    compressed = check_for_uncompressed(fst_file);
    if (compressed){
            // prep input
            const char *ffile = fst_file.c_str();
            Scoring_Pipeline(ffile);
        }
    else{
        std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(EXIT_FAILURE);
        }    
    return 0;
}

