#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>

using std::cout;
using std::string;
using std::fstream;
using std::stringstream;
using std::getline;
using std::stoi;
using std::stof;


// compile with -lz

/*
    find windows where the coverage is 0
    standard using: 17
*/

string rename_file(string &file_name){

    string outname;

    // get basename of file
    std::string file_basename =
        file_name.substr(file_name.find_last_of("/\\") + 1);

    // add suffix
    outname = "No_Coverage_" + file_basename.substr(0, file_basename.find_last_of('.'));
    return outname;
}

void create_smth_file(const char *d_file, string &output){

    int buffer_line_size = 100, position, depth; // read in size
    struct gzFile_s *gzf;      
    char line_buffer[buffer_line_size]; // char array to hold line
    gzf = gzopen(d_file, "r");  
    bool eol = false;                        
    memset(line_buffer, '\0', buffer_line_size); 
    string line, chrom, tmp, cur_chrom;
    fstream outfh;
    int start = -1, end = -1;
    bool need_start = true; // will keep track of whether we are expanding or adding window
    outfh.open(output, std::ios::out);

    // let's put a header
    outfh << "Sequence\tStart_Pos\tEnd_Pos\tSize\n";
    
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
            if (cur_chrom != chrom){
                if (!need_start){ // close current window when moving to chromosomes
                    outfh << cur_chrom << '\t' << start << '\t' << end << '\t' << end - start << '\n';
                    start = -1;
                    end = -1;
                    need_start = true;
                }
                cur_chrom = chrom;
            }
            if (depth == 0){
                if (need_start){
                    start = position;
                    end = position; // size of 1 for now
                    need_start = false;
                }
                else{
                    end = position;
                }
            }
            else{
                if (!need_start){
                    // again, close window
                    outfh << chrom << '\t' << start << '\t' << end << '\t' << end - start << '\n';
                    start = -1;
                    end = -1;
                }
                need_start = true; // just to be safe
            }
        }
     } while (!gzeof(gzf) && !eol);


     outfh.close();
}


int main(int argc, char *args[])
{
    // params
    string depth_file, output_file, message, file_basename;
    bool compressed;

    // in case of help or error
    message = "./depth_file -d file.depth.gz\n";

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

    }

    if (depth_file == "")
    {
        cout << "depth file is required\n";
        cout << message;
        exit(EXIT_FAILURE);
    }

    
    // start processing
    if (depth_file != ""){
        // check if gzipped b/c that's how i wrote the program
        output_file = rename_file(depth_file);
        const char *dfile = depth_file.c_str();
        create_smth_file(dfile, output_file);
        exit(EXIT_SUCCESS);  
            }
        else{
            std::cout << "Need to write an overload function to work with uncompressed files. Sorry\n";
            exit(EXIT_FAILURE);
        }
    return 0;
}