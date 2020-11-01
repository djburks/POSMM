# include <fstream>
# include <string>
# include <sstream>
# include <iostream>
# include <unordered_map>
# include <vector>
# include <math.h>
# include <algorithm>
# include <complex>

// Globals 

std::vector<std::vector<double *>> fmasta;
std::vector<std::vector<double *>> rmasta;
std::vector<int> mflen;  
int klim = (int)pow(4,12) + 1; 
std::vector<std::vector<double>> smmodel(5,std::vector<double> (klim,1));


// Transition nucleotide index.

std::unordered_map<char,int> nucloc = {
    {'A',1},
    {'T',2},
    {'G',3},
    {'C',4}
};

std::unordered_map<char,int> rnucloc = {
    {'A',2},
    {'T',1},
    {'G',4},
    {'C',3}
};

// Fasta import and filter function.
// Skips id lines, returns only standard ATGC nucleotides in uppercase.

std::string fasta(std::string infile) {
    std::ifstream fna(infile);
    std::string line;
    std::string genome;
    while(std::getline(fna,line)) {
        if (line[0] == '>') {
            continue;
        }
        for(char const &c: line) {
            if (c == 'A' || c == 'a') {genome+=('A');}
            if (c == 'T' || c == 't') {genome+=('T');}
            if (c == 'G' || c == 'g') {genome+=('G');}
            if (c == 'C' || c == 'c') {genome+=('C');}
        }
        
        
    }
    return genome;
}

// Kmer index functions.

int index(std::string k) {
    int idex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        idex = idex + ((pow(4,len-i-1))*nucloc[k[i]]);
    }
    return idex;
}

// Build multidimensional array model for storage of counts.
// 1-4 for transitional specific counts, 0 for total.

void smm(std::string genomefile, int order, std::vector<std::vector<double>> &model) {
    
    std::string line;
    std::string genome = "";
    std::ifstream infile(genomefile);
    while (std::getline(infile,line)) {
        if (line[0] == '>') {
            continue;
        }
        for(char const &c: line) {
            if (c == 'A' || c == 'a') {genome+=('A');}
            if (c == 'T' || c == 't') {genome+=('T');}
            if (c == 'G' || c == 'g') {genome+=('G');}
            if (c == 'C' || c == 'c') {genome+=('C');}
        }
    }

    // Iterate through genome, build counts, and generate probabilities.
    int giter = genome.length() - order;
    std::string r(order,'A');
    std::string e(order,'C');
    int reducer = index(r);
    int ender = index(e);
    int totalc = giter + (pow(4,order)*4);

    
    for(int i=0; i<giter; i++) {
        std::string&& kmer =  genome.substr(i,order);
        char tnuc = genome[i+order];
        int idex = (index(kmer) - reducer);
        model[nucloc[tnuc]][idex]++;
        model[0][idex]++;
    }
    

    for(int i=0; i<=(ender-reducer); i++) {
        model[4][i] = log10((model[4][i])/(model[0][i] + 3));
        model[3][i] = log10((model[3][i])/(model[0][i] + 3));
        model[2][i] = log10((model[2][i])/(model[0][i] + 3));
        model[1][i] = log10((model[1][i])/(model[0][i] + 3));
        model[0][i] = log10((model[0][i] + 3)/totalc);
    }

}

// Read in metafasta file.
// Store filtered reads in vector.

void readlibrary(std::string infile, int order) {
    //std::cout << "Indexing FASTA..." << std::endl;
    std::ifstream readfile(infile);
    std::string readline = "";
    std::string filterline = "";
    std::string fread;
    std::string rread;
    std::string r(order,'A');
    std::string e(order,'C');
    int reducer = index(r);
    int ender = index(e);
    int fidex = 0;
    int ridex =0;
    while (std::getline(readfile,readline)) {
        if (readline[0] == '>') {
            if (filterline != "") {
            	std::vector<double *> tfmasta,trmasta;
            	fread = filterline;
            	rread = filterline;
            	std::reverse(rread.begin(),rread.end());
            	std::string&& fchunk = fread.substr(0,order);
            	std::string&& rchunk = rread.substr(0,order);
            	fidex = 0;
            	ridex = 0;
            	for(int j=0; j<order; j++) {
                	fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                	ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
        	}
            	fidex = fidex - reducer;
        	ridex = ridex - reducer;
        	// Add to *vec
        	double* fob = &smmodel[0][fidex];
        	double* rob = &smmodel[0][ridex];
        	tfmasta.push_back(fob);
        	trmasta.push_back(rob);
        	for(int i=0;i<fread.length()-order;i++) {
            		std::string&& fchunk = fread.substr(i,order);
            		std::string&& rchunk = rread.substr(i,order);
        
           	fidex = 0;
            	ridex = 0;
            	for(int j=0; j<order; j++) {
               	fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                	ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
            	}
            	fidex = fidex - reducer;
            	ridex = ridex - reducer;
            
        	double* fob = &smmodel[nucloc[fread[i+order]]][fidex];
        	double* rob = &smmodel[rnucloc[rread[i+order]]][ridex];
        	tfmasta.push_back(fob);
        	trmasta.push_back(rob);
        	}
        	fmasta.push_back(tfmasta);
        	rmasta.push_back(trmasta);
        	mflen.push_back(fread.length());
        	filterline = "";
            }
        }
        else {
            for(char const &c: readline) {
                if (c == 'A' || c == 'a') {filterline+=('A');}
                if (c == 'T' || c == 't') {filterline+=('T');}
                if (c == 'G' || c == 'g') {filterline+=('G');}
                if (c == 'C' || c == 'c') {filterline+=('C');}
            }
        }
    }
    if (filterline != "") {
            	std::vector<double *> tfmasta,trmasta;
            	fread = filterline;
            	rread = filterline;
            	std::reverse(rread.begin(),rread.end());
            	std::string&& fchunk = fread.substr(0,order);
            	std::string&& rchunk = rread.substr(0,order);
            	int fidex = 0;
            	int ridex = 0;
            	for(int j=0; j<order; j++) {
                	fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                	ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
        	}
            	fidex = fidex - reducer;
        	ridex = ridex - reducer;
        	// Add to *vec
        	double* fob = &smmodel[0][fidex];
        	double* rob = &smmodel[0][ridex];
        	tfmasta.push_back(fob);
        	trmasta.push_back(rob);
        	for(int i=0;i<fread.length()-order;i++) {
            		std::string&& fchunk = fread.substr(i,order);
            		std::string&& rchunk = rread.substr(i,order);
        
           	int fidex = 0;
            	int ridex = 0;
            	for(int j=0; j<order; j++) {
               	fidex = fidex + (pow(4,order-j-1)*nucloc[fchunk[j]]);
                	ridex = ridex + (pow(4,order-j-1)*rnucloc[rchunk[j]]);
            	}
            	fidex = fidex - reducer;
            	ridex = ridex - reducer;
            
        	double* fob = &smmodel[nucloc[fread[i+order]]][fidex];
        	double* rob = &smmodel[rnucloc[rread[i+order]]][ridex];
        	tfmasta.push_back(fob);
        	trmasta.push_back(rob);
        	}
        	fmasta.push_back(tfmasta);
        	rmasta.push_back(trmasta);
        	mflen.push_back(fread.length());
        	filterline = "";
            }
    readfile.close();
    //std::cout << fmasta.size() << std::endl;
    return;
}



int rindex(std::string k) {
    int ridex = 0;
    int len = k.length();
    for(int i=0; i<=len; i++) {
        ridex = ridex + (pow(4,len-i-1)*rnucloc[k[i]]);
    }
    return ridex;
}

// Normalization for average 12th order scores based on read length across all models.

double normScore(float rawscore,int readlen) {
    float denom;
    denom = (-0.605185747917117*readlen) + 0.228066296087579;
    return rawscore/denom;
}

int main(int argc, char *argv[]) {
    
    std::string readfile = argv[2];
    int order = std::stoi(argv[3]);
    std::string scoretype = argv[4];
    std::string outfilenameprefix = argv[5];
    readlibrary(readfile,order);
    std::string r(order,'A');
    std::string e(order,'C');
    int reducer = index(r);
    int ender = index(e);
    
    
    
// For every genome in the linked genome fasta list, cycle through and build model.
// Use fasta lookups to generate score.
    std::vector<std::string> genomes;
    std::string genomeline;
    std::ifstream genomefastafile(argv[1]);
    while(std::getline(genomefastafile,genomeline)) {
        genomes.push_back(genomeline);
    }
    

    int genomenum = 0;
    std::vector<double> probline;
    std::vector<int> idline;
    int gid = -1;
    float cur_best;
    for(int i2=0; i2 < fmasta.size(); i2++) {
          probline.push_back(-999.9);
          idline.push_back(gid);
     }
     for (auto& genome: genomes) {
     	//Loop reset the SMM interiors.
     	for (auto &n : smmodel) {
     		std::fill(n.begin(),n.end(),1);
     	}
          gid++;
	  //std::cout << gid << std::endl;
          smm(genome,order,smmodel);
          for(int i2=0; i2 < fmasta.size(); i2++) {
                    double forprob = 0;
                    double revprob = 0;
                    for(int i3=0; i3 < fmasta[i2].size(); i3++) {
                    	forprob += *fmasta[i2][i3];
                    	revprob += *rmasta[i2][i3];
                    }
                    cur_best = std::max(forprob,revprob);
                    if (cur_best > probline[i2]) {
                    	probline[i2] = cur_best;
                    	idline[i2] = gid;
                    }
                }
      }
    std::ofstream outfile;
    outfile.open (outfilenameprefix + ".markov", std::ofstream::out);
    for(int i4=0;i4 < probline.size();i4++) {
    	outfile << probline[i4];
    	outfile << "\t";
    	outfile << genomes[idline[i4]];
    	outfile << "\n";
    }
    
    outfile.close();
    
    return 0;
}
