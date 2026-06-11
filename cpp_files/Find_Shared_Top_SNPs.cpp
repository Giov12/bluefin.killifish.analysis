#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <vector>
#include <map>

using std::cout;
using std::vector;
using std::string;
using std::fstream;
using std::stringstream;
using std::getline;
using std::stoi;
using std::stof;
using std::map;

struct Record{
    string chrom;
    int pos;
    float fst;

    Record(string c, int p, float f){
        chrom = c;
        pos = p;
        fst = f;
    }

    ~Record(){};
};


Record ParseLine(string &line){
    
    int pos;
    string tmp, chrom;
    float fst;
    stringstream ss(line);
    
    getline(ss, tmp, '\t');
    chrom = tmp;
    getline(ss, tmp, '\t');
    pos = stoi(tmp);
    getline(ss, tmp, '\n'); // todo : do a for loop in the future
    fst = stof(tmp);

    Record R = Record(chrom, pos, fst);

    return R;
}



void Find_Shared(vector<string> fst_files){

    // loop params
    string line;
    map<string, map<int, int>> snp_map;
    map<int, int> def_map;

    for (string &fstfile: fst_files){
        
        fstream fh;
        fh.open(fstfile, std::ios::in);
        
        if (fh.is_open()){
            while (getline(fh, line)){
                if (line == ""){
                    continue;
                }
                Record R = ParseLine(line);
                if (snp_map.count(R.chrom)){
                    if (snp_map[R.chrom].count(R.pos)){
                        snp_map[R.chrom][R.pos] += 1;
                    }
                    else{
                        snp_map[R.chrom][R.pos] = 1;
                    }
                }
                else{
                    snp_map[R.chrom] = def_map;
                    snp_map[R.chrom][R.pos] = 1;
                }
            }
            fh.close();
        }
    }

    // now to write to a file
    string outname = "Shared_Top_SNPs.tsv";
    fstream ofh;
    ofh.open(outname, std::ios::out);

    for (auto const &cmap: snp_map){
        for (auto &pos : cmap.second){
            if (pos.second == 4){
                ofh << cmap.first << '\t' << pos.first << '\n';
            }
        }
    }

    ofh.close();
}


int main(int argc, char *args[])
{
    // params
    string f1, f2, f3, f4, message;

    // in case of help or error
    message = "./Find_Shared_Top_SNPs -f1 fst_file1 -f2 fst_file2 -f3 fst_file3 -f4 fst_file4 [ALL required]\n";

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
        else if (param == "-f1")
        {
            f1 = string(args[i + 1]);
        }
        else if (param == "-f2")
        {
            f2 = string(args[i + 1]);
        }
        else if (param == "-f3")
        {
            f3 = string(args[i + 1]);
        }
        else if (param == "-f4")
        {
            f4 = string(args[i + 1]);
        }
    }

    if (f1 == "" || f2 == "" || f3 == "" || f4 == "")
    {
        cout << "Please provide all 4 files.\n";
        cout << message;
        exit(EXIT_FAILURE);
    }

    vector<string> snp_files = {f1, f2, f3, f4};

    Find_Shared(snp_files);

    return 0;
}

