#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <iomanip>

/*
Had to compile using the following flags -lz
*/

bool check_for_uncompressed(std::string &file_name){

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

std::string rename_file(const char *infile){

    std::string outname;
    std::string file_name = infile;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add suffix
    outname = "Dxy_" + file_basename.substr(0, file_basename.find_last_of('.'));
    return outname;
}

/*
    takes a column such as A:T:C:G:N:D
    and returns a vector of size 4 containing the counts

    [A count, T count, C count, G count]
*/
std::vector<int> parse_column(std::string &allele_cnt){

    int counter = 0;
    std::string num = "";
    std::vector<int> genotypes{6};
    int val;

    for (auto &c : allele_cnt){
        if (c == ':'){
            genotypes[counter] = std::stoi(num);
            num = "";
            counter++;
        }
        else{
            num += c;
        }
    }

    // add last val
    genotypes[5] = std::stoi(num);

    return genotypes;
}

/*
    compares two vectors of size 6 containing allele counts
                    [A, T, C, G, N, INDEL]
    and returns whether the site is a possible candidate of interest
*/

bool is_biallelic(std::vector<int> &pop1, std::vector<int> &pop2){

    // if either are not diploid, return false
    int p1_cnt;
    int p2_cnt;
    int g1_cnt = 0;
    int g2_cnt = 0;
    bool biallelic;

    for (int i = 0; i < 6; i++){
        p1_cnt = pop1[i];
        p2_cnt = pop2[i];
        if (p1_cnt != 0){
            g1_cnt++;
        }
        if (p2_cnt != 0){
            g2_cnt++;
        }
    }

    if (g1_cnt > 2 || g2_cnt > 2){
        biallelic = false;
    }
    // saw a few cases where this is true
    else if (g1_cnt == 0 || g2_cnt == 0){
        biallelic = false;
    }
    else if (g1_cnt == 1 && g2_cnt == 1){
        biallelic = false; // may lose fixed sites, but this should be few
    }
    else{
        biallelic = true;
    }

    return biallelic;
 
}

// formula from https://onlinelibrary.wiley.com/doi/pdf/10.1111/eva.13488

float calc_dxy(std::vector<int> &pop1_genos, std::vector<int> &pop2_genos){
    
    float dxy;
    float p1_cnt = 0.0;
    float p2_cnt = 0.0;
    float p1_tot = 0.0;
    float p2_tot = 0.0;
    float p1, p2, q1, q2;
  
    for (int i = 0; i < 6; i++){
        p1_tot += static_cast<float>(pop1_genos[i]);
        p2_tot += static_cast<float>(pop2_genos[i]);
    }
    for (int i = 0; i < 6; i++){
        p1_cnt = static_cast<float>(pop1_genos[i]);
        p2_cnt = static_cast<float>(pop2_genos[i]);
        if (p1_cnt > 0 || p2_cnt > 0){
            if (p1_cnt != 0){
                p1 = p1_cnt/p1_tot;
            }
            if (p2_cnt != 0){
                p2 = p2_cnt/p2_tot;
            }
            q1 = 1.0 - p1;
            q2 = 1.0 - p2;
            dxy = ((p1 * q2) + (q1 * p2));
            return dxy;
        }
    }
    // just in case
    return dxy;
}

void parse_compressed(const char *infile, std::string &outname){

    // general objects for streaming compressed data
    int buffer_line_size = 200;                   
    struct gzFile_s *gzf;                      
    char line_buffer[buffer_line_size];          
    gzf = gzopen(infile, "r");              
    bool eol = false;                           
    memset(line_buffer, '\0', buffer_line_size);
    // objects for this analysis
    std::string line, temp, pop1_geno, pop2_geno, chrom;  
    std::vector<int> p1_geno; // will hold columns A:T:C:G:N:INDEL
    std::vector<int> p2_geno;
    bool biallelic; // true or false
    int pos; // for location
    float dxy;
    // for output
    std::fstream outfile;
    outfile.open(outname, std::ios::out);

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
            std::stringstream ss(line);
            // parse out files
            std::getline(ss, chrom, '\t');
            std::getline(ss, temp, '\t');
            pos = std::stoi(temp);
            std::getline(ss, temp, '\t');
            std::getline(ss, pop1_geno, '\t');
            std::getline(ss, pop2_geno, '\n');

            // get arrays
            p1_geno = parse_column(pop1_geno);
            p2_geno = parse_column(pop2_geno);
            }
            biallelic = is_biallelic(p1_geno, p2_geno);
            if (biallelic){
                dxy = calc_dxy(p1_geno, p2_geno);
                outfile << std::fixed << chrom << '\t' << pos << '\t'
                        << pop1_geno << '\t' << pop2_geno
                        << '\t' << std::setprecision(5) << dxy << '\n';
            }

        } while (!gzeof(gzf) && !eol);
    // the end
    outfile.close();
}

int main(int argc, char *args[])
{
    // params
    std::string sync_file;
    std::string output_file;
    std::string message;
    std::string file_basename;
    bool compressed = false;

    // in case of help or error
    message = "./CalcDxySync -s sync_file [required]\n";

    // read args
    std::string param;

    for (int i = 0; i < argc; i++)
    {
        param = std::string(args[i]);

        if (param == "-h" || param == "--help")
        {
            std::cout << message;
            return 0;
        }
        else if (param == "-s" || param == "--sync_file")
        {
            sync_file = std::string(args[i + 1]);
        }
        else if (param == "-o" || param == "--output")
        {
            output_file = std::string(args[i + 1]);
        }
    }

    if (sync_file == "")
    {
        std::cout << "No sync file supplied.\n";
        std::cout << message;
        exit(EXIT_FAILURE);
    }
    
    // check if gzipped b/c that's how i wrote the program
    compressed = check_for_uncompressed(sync_file);

    // start processing
    if (compressed){
        const char *infile = sync_file.c_str();
        if (output_file == ""){
            output_file = rename_file(infile);
        }
        parse_compressed(infile, output_file);
    }
    else{
        std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(EXIT_FAILURE);
    }
    
    return 0;
}