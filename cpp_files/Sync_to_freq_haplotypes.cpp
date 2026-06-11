#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <filesystem>

using std::cout;
using std::vector;
using std::string;
using std::fstream;
using std::stringstream;
using std::getline;
using std::stoi;
using std::stof;

/*
    program to estimate the major and minor haplotypes based on the frequencies of the alleles

    NOTE: make sure to compile with -lz
*/

class Allele_Info{

    public:
        vector<char> alleles;
        vector<int> counts;
        char major, minor;
        bool monomorphic = false; // assume variant
        bool missing_data = true; // until proven otherwise

        Allele_Info(){}; // initialize as empty
        ~Allele_Info(){};

        void add_allele(char allele, int count){
            alleles.push_back(allele);
            counts.push_back(count);
        }

        // this is where we determine the major and minor alleles

        void determine_major(){

            if (missing_data){
                major = 'N';
                minor = 'N';
            }
            else{
                int highest_count = counts[0], second_highest = 0; // set as first allele
                major = alleles[0];

            if (alleles.size() == 1){
                monomorphic = true;
                minor = major; // to avoid an extra conditional
                }
            else{
                // used to only report the second most freq
                int second_highest = counts[1]; //
                minor = counts[1];
                }
            if (!monomorphic){
                for (int i = 1; i < alleles.size(); i++){
                    if (counts[i] > highest_count){
                        if (highest_count > second_highest){
                            // swap
                            second_highest = highest_count;
                            minor = alleles[i];
                        }
                        highest_count = counts[i];
                        major = alleles[i];
                        }
                    else if (counts[i] > second_highest){
                        second_highest = counts[i];
                        minor = alleles[i];
                        }
                    }
                }    
            }
        }
};

bool check_for_uncompressed(char *file_name){

    bool Compressed;
    int str_len = std::strlen(file_name);

    if (file_name[str_len - 1 - 2] == '.' && file_name[str_len - 1 - 1] == 'g' && file_name[str_len - 1] == 'z'){
        Compressed = true;
    }
    else{
        Compressed = false;
    }
    
    return Compressed;
}

string rename_file(char *file_name){

    string outname;

    // get basename of file
    string file_basename = std::filesystem::path(file_name).filename();

    // renaming here
    outname = "Freq_Based_Haplotypes_" + file_basename.substr(0, file_basename.find_last_of('.')) + ".fa";
    return outname;
}

Allele_Info get_allele_info(string &sample_alleles){

    int count = 0;
    string tmp = "";
    char nuc;
    Allele_Info AI = Allele_Info(); // start off empty

    for (int i = 0; i < sample_alleles.length(); i++){
        if (sample_alleles[i] != ':'){
            tmp += sample_alleles[i];
        }
        else {
            count++;
            if (tmp != "0"){
                switch (count)
                {
                case 1:
                    nuc = 'A';
                    break;
                case 2:
                    nuc = 'T';
                    break;
                case 3:
                    nuc = 'C';
                    break;
                case 4:
                    nuc = 'G';
                    break;
                case 5:
                    nuc = 'N';
                    break;
                case 6:
                    nuc = 'D';
                    break;
                } // end of switch case
                int cnt = stoi(tmp);
                AI.add_allele(nuc, cnt);
                AI.missing_data = false; // we have data
            }
            tmp = ""; // reset
        }
    }

    // check last count
    if (tmp != "0"){
        int cnt = stoi(tmp);
        AI.add_allele('D', cnt);
    }

    // determine major/minor
    AI.determine_major();

    return AI;

}

int count_num_samples(string &line){
    /* count number of columns based on tabs
        forumula
        num columns = total tabs - 3 + 1
    */
   int tab_count = 0;
   for (int i = 0; i < line.length(); i++){
        if (line[i] == '\t'){
            tab_count++;
        }
   }

   return tab_count - 3 + 1; // add last sample
}

void create_haplotypes(char *sync_file, string output){

    // params for reading in data
    int buffer_line_size = 100, start_position, end_position, num_samples; // read in size
    struct gzFile_s *gzf;      
    char line_buffer[buffer_line_size]; // char array to hold line
    gzf = gzopen(sync_file, "r");  
    bool eol = false, first = true;                        
    memset(line_buffer, '\0', buffer_line_size); 
    string line, chrom, pos, tmp, cur_chrom;
    // params for the output
    vector<string> haplotypes; // will be of size 2n, where n == number of columns
    fstream outfh;
    outfh.open(output, std::ios::out);
    string major_haplo, minor_haplo;

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
            if (first){
                // get num of samples
                num_samples = count_num_samples(line);
                for (int i = 0; i < num_samples * 2; i++){
                    haplotypes.push_back(""); // fill with empty strings
                }
            }
            // now to get the fields
            stringstream ss(line);
            getline(ss, chrom, '\t');
            getline(ss, tmp, '\t');
            if (first){
                start_position = stoi(tmp);
                // now we can set first to false
                first = false;
            }
            end_position = stoi(tmp); // eventually will become the final end
            getline(ss, tmp, '\t'); // this is the reference allele
            // now loop through
            for (int j = 0; j < num_samples - 1; j++){
                getline(ss, tmp, '\t');
                Allele_Info AI = get_allele_info(tmp);
                haplotypes[j * 2] += AI.major; // major
                haplotypes[j * 2 + 1] += AI.minor;
            }
            // then get last one
            getline(ss, tmp, '\n');
            Allele_Info AI = get_allele_info(tmp);
            haplotypes[(num_samples * 2) - 2] += AI.major; // major
            haplotypes[(num_samples * 2) - 1] += AI.minor;
        }
     } while (!gzeof(gzf) && !eol);

     // now to write the output

    for (int i = 0; i < num_samples; i++){
        outfh << ">Sample" << i + 1 << "_Major_" << chrom << '_' << start_position << '-' << end_position << '\n';
        outfh << haplotypes[i * 2] << '\n';
        outfh << ">Sample" << i + 1 << "_Minor_" << chrom << '_' << start_position << '-' << end_position << '\n';
        outfh << haplotypes[i * 2 + 1] << '\n';
    }

     outfh.close();
}

int main(int argc, char *args[])
{
    // params
    string file_basename;
    char *sync_file, *output_file; // some empty char arrays
    bool compressed;

    // in case of help or error
    const char *message = "./Sync_to_freq_haplotypes -s sync_file.gz [required] -o output_name [optional]\n";

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
        else if (param == "-s")
        {
            sync_file = args[i + 1];
        }
        else if (param == "-o")
        {
            output_file = (args[i + 1]);
        }
    }

    if (sync_file[0] == '\0')
    {
        cout << "sync file is required\n";
        cout << message;
        exit(EXIT_FAILURE);
    }

    // a safety check for now
    compressed = check_for_uncompressed(sync_file);
    if (compressed != true){
        cout << "sync file must be compressed\n";
        cout << message;
        exit(EXIT_FAILURE);
    }
    
    // start processing
    if (output_file[0] == '\0'){
        string output_file = rename_file(sync_file);
    }
    create_haplotypes(sync_file, output_file);

    return 0;
}