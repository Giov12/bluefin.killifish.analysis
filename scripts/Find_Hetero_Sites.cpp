#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>

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
    outname = "Hetero_Sites_" + file_basename.substr(0, file_basename.find_last_of('.'));
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
    std::vector<int> genotypes(4);
    int val;

    for (auto &c : allele_cnt){
        if (c == ':'){
            genotypes[counter] = std::stoi(num);
            num = "";
            counter++;
            if (counter == 4){
                break;
            }
        }
        else{
            num += c;
        }
    }

    return genotypes;
}

/*
    compares two vectors of size 4 containing allele counts
                    [A, T, C, G]
    and returns whether the site is a possible candidate of interest
*/

bool is_differ(std::vector<int> &pop1, std::vector<int> &pop2, int &min_cov, bool diploid){

    // pop 1 should be homozygous while pop 2 should be hetero or homo alternative
    int pop1_cnt, pop2_cnt;
    int evidence = 0;
    bool diff = false;
    int num_geno1 = 0;
    int num_geno2 = 0;

    for (int i = 0; i < 4; i++){
        pop1_cnt = pop1[i];
        pop2_cnt = pop2[i];
        if (pop1_cnt == 0 && pop2_cnt != pop1_cnt){
            if (pop2_cnt >= min_cov){
                diff = true;
            }
            num_geno2++;
        }
        else if (pop2_cnt == 0 && pop2_cnt != pop1_cnt){
            if (pop1_cnt >= min_cov){
                diff = true;
            }
            num_geno1++;
        }
        else if (pop2_cnt > 0 && pop1_cnt > 0){
            num_geno1++;
            num_geno2++;
        }
        if (pop1_cnt > pop2_cnt && pop1_cnt/2 > pop2_cnt){
            if (pop1_cnt >= min_cov){
                evidence++;
            }
        }
        else if (pop2_cnt > pop1_cnt && pop1_cnt == 0){
            if (pop2_cnt >= min_cov){
                evidence++;
            }
        }
    }

    if (evidence >= 2 && diff){
        if (!diploid){
            return true;
        }
        else if (num_geno1 > 2 || num_geno2 > 2){
            return false;
        }
        else {
            return true;
        }
    }
    else{
        return false;
    }
}

void parse_compressed(const char *infile, int min_cov, std::string &outname, bool diploid){

    int buffer_line_size = 200;                   // read in size
    struct gzFile_s *gzf;                        // pointer to file
    char line_buffer[buffer_line_size];          // char array to hold line
    gzf = gzopen(infile, "r");               // file handle
    bool eol = false;                            // end of file
    std::string line, temp, pop1_geno, pop2_geno;  
    memset(line_buffer, '\0', buffer_line_size); // fill in nulls
    std::vector<int> p1_geno;
    std::vector<int> p2_geno;
    bool difference;
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
            // skip over to cols of interest
            for (int j = 0; j < 3; j++){
                std::getline(ss, temp, '\t');
            }
            // now put cols of interest
            std::getline(ss, pop1_geno, '\t');
            std::getline(ss, pop2_geno, '\n');

            // get arrays
            p1_geno = parse_column(pop1_geno);
            p2_geno = parse_column(pop2_geno);
            }
            difference = is_differ(p1_geno, p2_geno, min_cov, diploid);
            if (difference){
                outfile << line;
            }


        } while (!gzeof(gzf) && !eol);
}

int main(int argc, char *args[])
{
    // params
    std::string sync_file;
    std::string output_file;
    std::string message;
    std::string file_basename;
    int min_cov = 1;
    bool compressed = false;
    bool diploid = false;

    // in case of help or error
    message = "./Find_Hetero_Sites -s sync_file [required]\n";

    // read args
    std::string param;

    for (int i = 0; i < argc; i++)
    {
        param = std::string(args[i]);

        if (param == "-h" || param == "--help")
        {
            std::cout << "This?\n";
            std::cout << message;
            return 0;
        }
        else if (param == "-s" || param == "--sync_file")
        {
            sync_file = std::string(args[i + 1]);
        }
        else if (param == "-c" || param == "--cov")
        {
            min_cov = std::stoi(std::string(args[i + 1]));
        }
        else if (param == "-o" || param == "--output")
        {
            output_file = std::string(args[i + 1]);
        }
        else if (param == "-d" || param == "--diploid")
        {
            diploid = true;
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
        parse_compressed(infile, min_cov, output_file, diploid);
    }
    else{
        std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
        exit(EXIT_FAILURE);
    }
    
    return 0;
}